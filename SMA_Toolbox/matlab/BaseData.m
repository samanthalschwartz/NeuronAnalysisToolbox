

classdef BaseData < Pickle & GUIBuilder
    properties
        globalTBounds; %1-based [minT maxT] no useful data outside these bounds.  This allows non-useful frames to be fully removed from the data presentation.

        %Data region selections

        ROI={}; % cell array of saved ROIs (1-based) [xmin xmax ymin ymax Lmin Lmax tmin tmax]
        ROIname={}; %Name for regions of interest
    end

    %These properties are not saved
    properties (Access=protected, Transient=true)
        dirtyData=false; %Marks fact that data still needs to be written to rawDataPath.

        %Lazily loaded properties
        raw_frames_; 
        frames_;
        %Lazily loaded property flags
        raw_frames_loaded=false; 
        frames_loaded=false;
    end
   
    methods
        function save(obj)
            obj.assertInitialized();
            obj.updateWaitbar(0,'Saving');
            if isempty(obj.globalTBounds) && isfield(obj.preservedProperties,'globalTBounds');
                obj.globalTBounds = obj.preservedProperties.globalTBounds;
                obj.globalTBounds(1)=min(obj.globalTBounds(1), obj.nFrames); %Make sure globalTBounds stays valid               
                obj.globalTBounds(2)=min(obj.globalTBounds(2), obj.nFrames); %Make sure globalTBounds stays valid               
            end
            save@Pickle(obj);
            obj.updateWaitbar(1);
        end

        function resetObject(obj)
            %Prevent clearing of obj.guiFig and obj.waitbarH
            guiFig = obj.guiFig;
            waitbarH = obj.waitbarH;
            resetObject@Pickle(obj);
            obj.guiFig = guiFig;
            obj.waitbarH = waitbarH;
        end

         
        function setPreservedProperties(obj, defaultProperties)
            setPreservedProperties@Pickle(obj, defaultProperties);
            if obj.initialized && isempty(obj.ROI)
                obj.filterValidROI();
            end
        end

        %% ROI Management
        function addROI(obj, roi_in, roi_name)
            % Add a new ROI to the list of current ROI's and check that it is
            % valid
            %
            % IN:
            %   roi - [xmin xmax ymin ymax Lmin Lmax tmin tmax] (1-based indexs (not real unit values))
            if nargin==2
                roi_name=GUIBuilder.nextUnusedName(obj.ROIname, 'ROI%i', length(obj.ROI)+1);
            end
            roi=obj.checkROI(roi_in, roi_name);
            %check for duplicate
            for n=1:length(obj.ROI)
                if all(obj.ROI{n}==roi) || strcmp(obj.ROIname{n}, roi_name)
                    error('SPData:addROI',['Already have nameed roi:', roi_name]);
                end
            end
            %append to list
            obj.ROI=[obj.ROI roi];
            obj.ROIname=[obj.ROIname roi_name];
            obj.dirty=true;
        end

        function deleteROI(obj, roi_idx)
            if roi_idx>0 && roi_idx<=length(obj.ROI)
                obj.ROI(roi_idx)=[];
                obj.ROIname(roi_idx)=[];
                obj.dirty=true;
            end
        end
        
        function clearROI(obj)
            obj.ROI={};
            obj.ROIname={};
        end
        
        function modifyROI(obj, roi_idx, new_roi, new_roi_name)
            if nargin<4
                new_roi_name=obj.ROIname{roi_idx};
            end
            if roi_idx>0 && roi_idx<=length(obj.ROI)
                if strcmp(new_roi_name, obj.ROIname{roi_idx})
                    if all(new_roi==obj.ROI{roi_idx})
                        return %No changes
                    end
                    new_roi=obj.checkROI(new_roi);
                else
                    new_roi=obj.checkROI(new_roi,new_roi_name);
                end
                obj.ROI{roi_idx}=new_roi;
                obj.ROIname{roi_idx}=new_roi_name;
                obj.dirty=true;
            end
        end



        %% File path helper functions
        

        function [file_path, file_name, file_ext] = getROIFileNameParts(obj, roi_in, class_name)
            % This defines how the subsidiary classes are named for a particular ROI relative to this
            % SPData file.  This makes a common notaion allowing us to look for all relavent files for
            % that ROI.  This returns the same format that fileparts does. [path, name, ext]
            if isscalar(roi_in) && roi_in>=1 && roi_in<=length(obj.ROI)
                roi_name = obj.ROIname{roi_in};
            elseif ischar(roi_in)
                roi_name = roi_in;
            end
            file_path = obj.getFilePath(class_name);
            file_name = sprintf('%s_%s',obj.saveFileBaseName, roi_name);
            file_ext = eval([class_name '.saveFileExt']); %Have to do this since classes are not first class objects in matlab            
        end

        function file_paths=getROIFiles(obj,roi_in, class_name, selectSingle)
            %Get a cellarray of paths to files
            % [in]
            % selectSingle - boolean [defualt:false].  Force slection of just a single filename using
            %                           user input dialog if necessary
            [~,roi_idx] = obj.getROI(roi_in);
            roi_name = obj.ROIname{roi_idx};
            [class_path, roi_file_name, file_ext] = obj.getROIFileNameParts(roi_idx, class_name);
            pattern = sprintf('%s*%s',roi_file_name, file_ext);
            existing_paths = Pickle.listExistingFileNames(class_path,pattern);
            if nargin==3 || ~selectSingle || isempty(existing_paths) % return all files found
                file_paths = existing_paths;
            elseif numel(existing_paths)==1
                file_paths = existing_paths{1};
            else % Force user to choose just one
                dataformats = eval(sprintf('%s.SaveableDataFormats',class_name));
                title = sprintf('Load an Existing %s Filename for ROI:%s',class_name, roi_name);
                file_paths = Pickle.selectExistingFileName(class_path, pattern, dataformats, title);
            end
        end
        
        function clearROIFiles(obj, roi_in, class_name, confirm)
            % clear a class type of Pickle files associated with this .spdata for a given ROI.
            % [IN]
            %  roi_name - a name or an index
            %  class_name - char - a string giving
            %  confirm - [optional] boolean if we should warn.  Default: true
            if nargin==3
                confirm = true;
            end
            file_paths=obj.getROIFiles(roi_in, class_name);
            roi_name = obj.ROIname{roi_in};
            if isempty(file_paths) && obj.inGui
                warndlg(sprintf('No %s files for "%s" in directory: "%s"',class_name,roi_name,obj.getFilePath(class_name)),...
                                'Clear Files Warning');
                return
            end
            if confirm 
                msg = sprintf('Clear %i %s files saved for %s?',length(file_paths),class_name,roi_name);
                heading = sprintf('Clear %s File',class_name);
                switch questdlg(msg,heading,'Yes','No','Cancel','Cancel');
                    case {'No','Cancel'}
                        return
                end
            end
            delete(file_paths{:});
        end
    end


    methods (Access=protected)
        function updateFromPreservedProperties(obj)
            % Update persistant properties using the obj.preservedProperties value as the default if no other value has been set.
            % This is called on a save of a partially filled in object.
            updateFromPreservedProperties@Pickle(obj);
        end

        function filterValidROI(obj)
            %Try to preserve ROIs from previous file if usefull
            if isfield(obj.preservedProperties,'ROI')
                ROIs=obj.preservedProperties.ROI;            
                for n=1:length(ROIs)
                    roi=ROIs{n};
                    roi_name=obj.preservedProperties.ROIname{n};
                    if roi(1)<obj.sizeX && roi(3)<obj.sizeY
                        obj.addROI(roi(1:4),roi_name);
                    end
                end
            end
        end
       

    end

    methods(Abstract=true)
        roi_size=getROIsize(obj, roi_in);
        
        roi=getROI(obj, roi_in);
    end %Public abstract methods
    
    methods(Abstract=true, Access=protected)
        roi=checkROI(obj, roi_in, roi_name);
        
%         saveRawData(obj);
    end
end %classdef
