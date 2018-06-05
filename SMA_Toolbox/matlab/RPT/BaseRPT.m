% BaseRPT.m - The base class for 2D and 3D RPT tracking
% 05/2015
% Mark J. Olah (mjo@cs.unm.edu)
%

classdef BaseRPT < Pickle & GUIBuilder
    %The base class for 2D and hyperspectral tracking
    properties (Constant=true, Hidden=true)
        %The obj.phaseIdx property determines the highest phase of the tracking
        % that has been _completed_.  These are the values obj.phaseIdx can refer to.
        % The public dependant property obj.phase gives the current phase as a string from this cellarray
        N_PHASES=7;
        PHASES={ 'Invalid',...           % [phaseIdx=1] No data or ROI
                 'Initialized',...       % [phaseIdx=2] Data and ROI set and calibrated frames ready to track
                 'MaximaFound',...       % [phaseIdx=3] Data filtered and maxima identified
                 'MaximaFiltered',...    % [phaseIdx=4] The maxima have been thresholded
                 'EmittersLocalized',... % [phaseIdx=5] The candidate emitters ROI's have been localizaed
                 'EmittersFiltered',...  % [phaseIdx=6] The localizations have been filtered for wuality
                 'Tracked'...            % [phaseIdx=7] The tracks have been formed
                };
        PHASE_TITLES={ '',...
                        '',...
                        'Find Maxima',...
                        'Filter Maxima',...
                        'Localize Emitters',...
                        'Filter Emitters',...
                        'Track'...
                      };
        EmitterTableTitle='Internal Emitters Format';
        EmitterTableDimensionNames={'Internal Emitter Index', 'Parameter'}; 
                  
        LocalizationTableTitle='Filtered Localizations';
        LocalizationTableDimensionNames={'Localization', 'Parameter'}; 

        TrackTableTitle='Tracks';
        TrackTableDimensionNames={'Localization', 'Parameter'}; 

   end

   properties (Constant=true, Hidden=true, Abstract=true)
        DataClass; % A string with the name of the data class, e.g., 'SPData' 'HSData'
        DataFileExt; % The extsion of the save file type of the data class, e.g., '.spdata', .'hsdata'

        %The internal emitter format description
        NEmitterColumns;
        EmitterColumnNames;
        EmitterColumnUnits;
        % This defines a common localization format for other processes that need access to the localization list
        % we keep units in pixels as these units are most natural for the localization phase        
        NLocalizationColumns;
        LocalizationColumnNames;
        LocalizationColumnUnits;
        LocalizationColumnDescriptions;
        
        NTrackColumns;
        TrackColumnNames;
        TrackColumnUnits;
        TrackColumnDescriptions;
    end

    properties (Abstract=true)
        ParamsFindMaxima;       %[phaseIdx==3]
        ParamsFilterMaxima;     %[phaseIdx==4]
        ParamsLocalizeEmitters; %[phaseIdx==5]
        ParamsFilterEmitters;   %[phaseIdx==6]
        ParamsTrack;            %[phaseIdx==7]
    end

    properties
        %These properties are common the the 2D and HS tracking

        % ROI - 2D Format:[min_x, max_x, min_y, max_y, min_t, max_t] - Same format as from SPData
        %     - HS Format:[min_x, max_x, min_y, max_y, min_L, max_L, min_t, max_t] - Same format as from HSData
        ROI; 
        
        % A name for this ROI.  This forms part of the filename.  To distinguish tracked ROIs within a dataset
        ROIname;
 
        % A 2D image of the data as combination of maximum and sum images
        sumImage; 

        %Paths - The relevent file paths
        %All Paths are relative paths from obj.workingDir
        Paths=struct(... 
             'saveFile',[],... %filename(w/extension) of the .rpt file where this will save to
             'data',[],... % base name for raw image file (w/ extension).
             'figures','Figures'...  %relative path from workingDir where figures should be saved
             );

      
        % The results structures are populated by each of the 5 computation phases 3-7.  The names and
        % are common between 2D and HS tracking, but the value types are often different to accound for
        % the added dimension in the HS data
        ResultsFindMaxima=struct(...
            'rawMaxima',[],... %Maxima from filtered image. Size 2D:[3,n] or HS:[4,n]  Coords [X Y T] or [X Y L T]??
            'rawMaximaVals',[]... %Value of filtered image at each maxima
            );

        ResultsFilterMaxima=struct(...
            'maxima',[],... %Maxima from filtered image. Size 2D:[3,n] or HS:[4,n]  Coords [X Y T] or [X Y L T]??
            'maximaVals',[],... %Value of filtered image at each maxima
            'filter',[],... %boolean if a raw maxima was kept or not
            'maximaThreshold',[],...
            'boxCoords',[],... % a BoxCoords object that has the coords for the localizationROIs
            'emitterImages',[],... %These are the sub-images to fit - later replace with data structure to handle different sizes
            'emitterFrames',[]...
            );
        
        ResultsLocalizeEmitters=struct(...
            'thetaInit',[],...    % Initial guess for theta value.
            'rawTheta',[],...     % Direct thetas from the fit. In local box x/y coords.
            'rawEmitters',[]...  % Using column format given by obj.EmitterColumnNames
            );
        
        ResultsFilterEmitters=struct(...
            'theta',[],... %Direct thetas from the fit after filtering for quality. In local box x/y coords.
            'filter',[],...%flag if a localization is filtered out.  0=keep 
            'filterDescriptions',[],...% Cell array of descriptions for filtering phases that each give a tag
            'rawIdx',[],... % size=[1,nLocalizations] - index of each localization back to the corresponding raw loaclization
            'localizations',[]... % The filterd localizations.  Using column format given by obj.LocalizationColumns
            );
        
        ResultsTrack=struct(...
            'tracks',[],...
            'trackLengths',[]...
            );

        times; % a structre for storing timing info
    end

    properties (Abstract=true, Dependent=true)       
        % frameSize - [Req: phase>=2] - The size (2D or 3D) of an individual frame from the global ROI.
        % Note this is potentially smaller than the actual frame size.  This is the size of a single 
        % frame from this obj.ROI region of the raw data.
        frameSize;
        % nFrames - [Req: phase>=2] The number of frames in the global ROI again may be smaller than the actual number of frames in the raw dataset
        nFrames;
    end

    properties (Dependent=true)       
        phase;   % The current entry from obj.PHASES corresponding to internal index
        nScales; % The number of filtering scales >=1
        nMaxima; % [Req: phase>=4] The number of maxima identified after background filtering. 
        nLocalizations; % [Req: phase>=6] The number of localizations after filtering
        nTracks; % [Req: phase==7] The number of tracks in the final dataset
        
        emitterBoxCoords; % [Req: phase>=4] The box coords object with all boxes from phase=4.
        filteredEmitterBoxCoords; % [Req: phase>=6] The remaining box coords after filtering
    end

    properties (Access=protected)
        phaseIdx=1; %Index into global PHASES constant; 1='Invalid'
    end

    properties (Access=protected, Transient=true)
        frames_;
        frames_loaded=false;
        filtered_frames_;
        filtered_frames_loaded=false;
    end
    
    properties (Transient=true)
        data;
        boxxer;
        emitterModel;
        tracker;
        srimage; % A super-res rendering of the data if phaseIdx>=6
    end

    methods
        function load(obj, varargin)
            %Args same as for constructor, except they cannot be empty
            if isempty(varargin)
                error('RPT:load','No arguments given');
            end
            if ischar(varargin{1})
                [pathn, ~, ext] = fileparts(varargin{1});
                DataConstructor = str2func(obj.DataClass);
                if strcmp(ext, obj.saveFileExt)
                    % Option (2)
                    % Load from saved file
                    obj.resetObject();
                    s = load(varargin{1},'-mat');
                    obj.copyobj(s.obj); %copy object into this one.
                    obj.workingDir = pathn;
                    obj.data=DataConstructor(obj.getFilePath('data'));
                    obj.getFilePath('data');
                    obj.initialized = true;
                    obj.dirty = false;
                elseif strcmp(ext, obj.DataFileExt)
                    % Option (3)
                    obj.loadData(DataConstructor(varargin{1}),varargin{2:end});
                else
                    error('RPT:load','Uknown file extension "%s"',ext)
                end
            elseif isa(varargin{1}, obj.DataClass)
                % Option (4)
                obj.loadData(varargin{:});
            end
        end

        function save(obj)
            obj.assertInitialized();
            obj.updateWaitbar(0,'Saving');
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

        function fs = getFrames(obj)
            %Lazy-loaded cached copy of data
            %Not a dependant property so it doesnt get default displayed with the variable viewer
            if obj.phaseIdx<2
                fs = []; % not initialized yet
                return
            end
            if ~obj.frames_loaded
                obj.frames_ = obj.data.getFrames(obj.ROI);
                obj.frames_loaded = true;
            end
            fs = obj.frames_;
        end

        function fs = getFilteredFrames(obj)
            %Lazy-loaded cached copy of filtered data
            %Not a dependant property so it doesnt get default displayed with the variable viewer
            if obj.phaseIdx<2
                fs = []; % not initialized yet
                return
            end
            if ~obj.filtered_frames_loaded
                obj.filterFrames();
            end
            fs = obj.filtered_frames_;
        end

        function p=getParams(obj, idx)
            switch idx
                case 3; p = obj.ParamsFindMaxima;
                case 4; p = obj.ParamsFilterMaxima;
                case 5; p = obj.ParamsLocalizeEmitters;
                case 6; p = obj.ParamsFilterEmitters;
                case 7; p = obj.ParamsTrack;
                otherwise; p = [];
            end
        end

        function setParams(obj, idx, new_params)
            params = obj.getParams(idx);
            if isempty(params); return; end
            nfs = fieldnames(new_params);
            for i = 1:length(nfs)
                params.(nfs{i}) = new_params.(nfs{i});
            end
            switch idx
                case 3; obj.ParamsFindMaxima = params;
                case 4; obj.ParamsFilterMaxima = params;
                case 5; obj.ParamsLocalizeEmitters = params;
                case 6; obj.ParamsFilterEmitters = params;
                case 7; obj.ParamsTrack = params;
            end
            obj.dirty=true;
        end

        function setDefaultParams(obj, propertiesStruct)
            % Reset object to phase 2 and set all Params structures
            if isa(propertiesStruct, 'BaseRPT');
                propertiesStruct = propertiesStruct.getParamStruct();
            end
            defaultableProperties={'ParamsFindMaxima', 'ParamsFilterMaxima', 'ParamsLocalizeEmitters', 'ParamsFilterEmitters', 'ParamsTrack'};
            obj.resetPhase(2);
            for p=defaultableProperties
                obj.(p{1}) = propertiesStruct.(p{1});
            end
            obj.phaseIdx = 2;
        end


        function r=getResults(obj, idx)
            switch idx
                case 3; r = obj.ResultsFindMaxima;
                case 4; r = obj.ResultsFilterMaxima;
                case 5; r = obj.ResultsLocalizeEmitters;
                case 6; r = obj.ResultsFilterEmitters;
                case 7; r = obj.ResultsTrack;
                otherwise; r = [];
            end
        end

        function autoTrack(obj, endPhase)
            %
            % automatically complete necessary tracking phases up to a
            % desired endPhase
            %
            % In:
            %  endPhase - [default=7] The desired phase to complete the
            %                         tracking to
            if nargin<2
                endPhase = length(obj.PHASES);
            end
            if obj.phaseIdx == 1
                error('RPT:autoTrack','Object is invalid (phaseIdx=1)');
            end
            methods = {[],[],@obj.findMaxima,@obj.filterMaxima,@obj.localizeEmitters,...
                             @obj.filterEmitters,@obj.trackEmitters};
            obj.keepWaitBar = true;
            try
                for i = 3:obj.N_PHASES
                    if obj.phaseIdx<i && endPhase>=i
                        methods{i}();
                    end
                end
            catch err
                %Close waitbars and report error to caller
                obj.keepWaitBar = false;
                close(obj.waitbarH);
                rethrow(err);
            end
            obj.keepWaitBar = false;
            close(obj.waitbarH);
            if obj.phaseIdx ~= endPhase
                error('RPT:autoTrack','Could not complete tracking to desired endPhase:%s.  Current phase: %s',obj.PHASES{endPhase},obj.phase);
            end             
        end
 
        function resetPhase(obj, phaseIdx)
            % Reset object so that phaseIdx is not yet complete
            if nargin == 1 || phaseIdx < 2
                phaseIdx = 2;
            elseif phaseIdx >= 7
                warning('RPT:resetPhase','Can only reset to phaseIdx < N_PHASES=%i',obj.N_PHASES);
                return;
            elseif phaseIdx > obj.phaseIdx-1
                warning('RPT:resetPhase','Current phase is %i no need to reset to phase%i',obj.phaseIdx,phaseIdx);
                return;
            end
            obj.setPhase(phaseIdx-1);
        end

        function rerunPhase(obj, phaseIdx)
            % Reset object phaseIdx-1 and re-run phase corresponding to phaseIdx.
            %
            % This is like a combination of a reset and auto track that ensures a particular phase
            % is explicitly run.  This is usefull for testing a phase parameters change, etc.
            if  phaseIdx < 3  || phaseIdx > 7
                warning('RPT:resetPhase','Can only rerun computational phases 3-7');
            elseif  phaseIdx-1 >= obj.phaseIdx
                obj.autoTrack(phaseIdx);
            else                
                obj.setPhase(phaseIdx-1);
                obj.autoTrack(phaseIdx);
            end
        end

        %% Data Accessing Methods
        function boxCoords = getBoxCoords(obj)
            % The boxCoords object holds all the information on box coordinate and indexs to look up boxes
            % by frame scale and size
            obj.checkPhase(4);
            boxCoords = obj.ResultsFilterMaxima.boxCoords;
        end

        
        function [E, col_names, col_units, col_desc] = getRawEmitters(obj)
            % Note: The rawEmitters output is mainly useful for internal plotting methods, and it includes all
            % emitter localization information in an internal pixel-based format and represents the raw unfiltered
            % data which in general should not be used unless there is a clear application for them.  Please use
            % the getLocalizations() and related methods if you want a list of localizations.
            %
            %The rawEmitters is a varibale produced by the LocalizeEmitters (Phase=5).  It is in internal
            %emitters format an is a matrix of [Nemitters x NEmitterColumns] where each column is a
            %property simillar to Localizations, but normally in pixels
            % [out]
            %    E - Matrix of raw (unfiltered) emitters.  Each emitter is a row and features are columns.
            %    col_names - [optional] the value of obj.EmitterColumnNames for convenience
            %    col_units - [optional] the value of obj.EmitterColumnUnits for convenience
            %    col_desc  - [optional] the value of obj.EmitterColumnDescriptions for convenience
            obj.checkPhase(5);
            E = obj.ResultsLocalizeEmitters.rawEmitters;
            if nargout>1 % Provide these for convenience only
                col_names = obj.EmitterColumnNames;
                col_units = obj.EmitterColumnUnits;
                col_desc = obj.EmitterColumnDescriptions;
            end
        end

        function tableE = convertEmittersTable(obj, Es)
            % Note: The emitters output is mainly useful for internal plotting methods.  In general the
            % getLocalizations() is the preferred way to get information of the indvidual localizations and the
            % localizations are avilible in structure and table format too.
            tableE = array2table(Es);
            tableE.Properties.Description = obj.EmitterTableTitle;
            %Add a Track index column to the table output.
            tableE.Properties.VariableNames = obj.EmitterColumnNames;
            tableE.Properties.VariableDescriptions = obj.EmitterColumnDescriptions;
            tableE.Properties.VariableUnits = obj.EmitterColumnUnits;
        end
        
        function tableE = getRawEmittersTable(obj)
            % Note: The emitters output is mainly useful for internal plotting methods.  In general the
            % getLocalizations() is the preferred way to get information of the indvidual localizations and the
            % localizations are avilible in structure and table format too.
            % [out]
            %    tableE - Matrix of raw (unfiltered) emitters.  Each emitter is a row and features are columns.
            tableE = obj.convertEmittersTable(obj.getRawEmitters);
        end
        
        function exportRawEmitters(obj, fileName)
            % Note: The emitters output is mainly useful for internal plotting methods.  In general the
            % getLocalizations() is the preferred way to get information of the indvidual localizations and the
            % localizations are avilible in structure and table format too.
            %
            % Saves the raw emitters to a .csv or .xls file
            % [in]
            %  filename - path to filename to save.  Should end in .csv or .xls
            writetable(obj.getRawEmittersTable(), fileName);
        end

        
        function [E, col_names, col_units, col_desc] = getEmitters(obj)
            % Note: The emitters output is mainly useful for internal plotting methods.  In general the
            % getLocalizations() is the preferred way to get information of the indvidual localizations and the
            % localizations are avilible in structure and table format too.
            %
            % The emitters is an internal varibale produced by the FiltersEmitters (Phase=6).  It is in internal
            % emitters format an is a matrix of [Nemitters x NEmitterColumns] where each column is a
            % property simillar to Localizations, but normally in pixels
            % [out]
            %    E - Matrix of emitters.  Each emitter is a row and features are columns.
            %    col_names - [optional] the value of obj.EmitterColumnNames for convenience
            %    col_units - [optional] the value of obj.EmitterColumnUnits for convenience
            %    col_desc  - [optional] the value of obj.EmitterColumnDescriptions for convenience
            obj.checkPhase(6);
            E = obj.ResultsFilterEmitters.emitters;
            col_names = obj.EmitterColumnNames;
            if nargout>1 % Provide these for convenience only
                col_names = obj.EmitterColumnNames;
                col_units = obj.EmitterColumnUnits;
                col_desc = obj.EmitterColumnDescriptions;
            end
        end
        
        function tableE = getEmittersTable(obj)
            % Note: The emitters output is mainly useful for internal plotting methods.  In general the
            % getLocalizations() is the preferred way to get information of the indvidual localizations and the
            % localizations are avilible in structure and table format too.
            % [out]
            %    tableE - Matrix of (filtered) emitters.  Each emitter is a row and features are columns.
            tableE = obj.convertEmittersTable(obj.getEmitters);
        end
        
        function exportEmitters(obj, fileName)
            % Note: The emitters output is mainly useful for internal plotting methods.  In general the
            % getLocalizations() is the preferred way to get information of the indvidual localizations and the
            % localizations are avilible in structure and table format too.
            %
            % Saves the emitters to a .csv or .xls file
            % [in]
            %  filename - path to filename to save.  Should end in .csv or .xls
            writetable(obj.getEmittersTable(), fileName);
        end
        
        function tableL = getLocalizationsTable(obj)
            % Localizations are sorted by frameIdx.  Localizations differ from Emitters in that the localizations
            % are reported in physical units and don't include internal indexing columns and other extraneous info.  
            % In general the localizations are a better format for later post-processing, and the Emitter format is
            % mainly used internally for plotting.  Note that the localizations variable is produced by the
            % getLocalizations method on demand from the Emitters internal format.
            % [out] 
            %   tableL - localizations in table format with nice column names and units specified   
            obj.checkPhase(6);
            L = obj.getLocalizations();
            tableL = array2table(L);
            tableL.Properties.DimensionNames = obj.LocalizationTableDimensionNames;
            tableL.Properties.Description = obj.LocalizationTableTitle;
            tableL.Properties.VariableNames = obj.LocalizationColumnNames;
            tableL.Properties.VariableDescriptions = obj.LocalizationColumnDescriptions;
            tableL.Properties.VariableUnits = obj.LocalizationColumnUnits;
        end

        function [structL, field_units, field_desc] = getLocalizationsStruct(obj)
            % Localizations are sorted by frameIdx.
            % [out] 
            %   structL - localizations in struct format with named parameters   
            %   field_units - [optional] the value of obj.TrackColumnUnits for convenience
            %   field_descriptions - [optional] the value of obj.TrackColumnDescriptions for convenience
            obj.checkPhase(6);
            tableL = obj.getLocalizationsTable();
            structL = table2struct(tableL,'ToScalar',true);
            if nargout>1 % Provide these for convenience only
                field_units = obj.TrackColumnUnits;
                field_desc = obj.TrackColumnDescriptions;
            end
        end

        function exportLocalizations(obj, fileName)            
            % Saves the localizations to a .csv or .xls file
            % [in]
            %  filename - path to filename to save.  Should end in .csv or .xls
            writetable(obj.getLocalizationsTable(), fileName);
        end
        
        function [Ts, col_names, col_units, col_desc] = getTracks(obj,trackIdxs)
            % The RPT track format is defined by the RPT constant properties:
            %  TrackColumnNames, TrackColumnUnits, and TrackColumnDescriptions
            %
            % This method returns tracks as a cell array of track matricies.  Each
            % track is an element in the cell array sorted by number of localizations 
            % then by first frame.  Each track matrix has columns defined by TrackColumnNames and
            % rows represent the sequence of localizations for the track.
            % [in]
            %     trackIdxs - [optional] An array of selected track indexes [default=all tracks]
            % [out] 
            %     Ts - A cell array of matricies.  One matrix for each track.  
            %               * rows = localization,
            %               * cols = named properties
            %     col_names - [optional] the value of obj.TrackColumnNames for convenience
            %     col_units - [optional] the value of obj.TrackColumnUnits for convenience
            %     col_desc  - [optional] the value of obj.TrackColumnDescriptions for convenience
            obj.checkPhase(7);
            Ts = obj.ResultsTrack.tracks;
            if nargin==2
                Ts = Ts(trackIdxs);%select tracks if requested
            end
            if nargout>1 % Provide these for convenience only
                col_names = obj.TrackColumnNames;
                col_units = obj.TrackColumnUnits;
                col_desc = obj.TrackColumnDescriptions;
            end
        end
        
        function tableT = getTracksTable(obj, varargin)
            % This is a reformatting of the RPT tracks format into a single table where all tracks have
            % been concatenated and distinguised by an extra trackID column.  This is slower an more inefficient
            % than the getTracks() method, but can be convienet as tables have the ability to store the names,
            % units, and descritptions of all columns and we take advantage of this.
            %
            % Note: if you have an RPT tracks cell array it can be converted to a table using the
            % static method BaseRPT.convertTracksTable().  This is how this method works.
            %
            % [in]
            %   trackIdxs - [optional] An array of selected track indexes [default=all tracks]
            % [out] 
            %   tableT - A table giving each tracks information using an inital column for trackID
            tableT = obj.convertTracksTable(obj.getTracks(varargin{:}),varargin{:});
        end

        function [structT, field_units, field_desc] = getTracksStruct(obj, varargin)
            % This is a reformatting of the RPT tracks format into a structure array where each element of the
            % array is a single track and the fields represent the columns in the normal RPT format.  This is
            % more inefficeint to use and is provided for convenience.
            %
            % Note: if you have an RPT tracks cell array it can be converted to a table using the
            % static method BaseRPT.convertTracksStruct().  This is how this method works.
            %
            % [in]
            %   trackIdxs - [optional] An array of selected track indexes [default=all tracks]
            % [out] 
            %   structT - tracks in structure array format.  Each element in the array is one track
            %                 the structure has named properties for each parameter
            %   field_units - [optional] the value of obj.TrackColumnUnits for convenience
            %   field_descriptions - [optional] the value of obj.TrackColumnDescriptions for convenience
            structT = obj.convertTracksStruct(obj.getTracks(varargin{:}));
            if nargout>1 % Provide these for convenience only
                field_units = obj.TrackColumnUnits;
                field_desc = obj.TrackColumnDescriptions;
            end
        end
        
        function exportTracks(obj, fileName, varargin)            
            % Saves the tracks to a .csv or .xls file.
            % This uses the same formate as getTracksTable() where an initial column gives the trackID
            % and all track matricies are then concatenated vertically.
            % [in]
            %  filename - path to filename to save.  Should end in .csv or .xls
            %  trackIdxs - [optional] An array of track indexes to save [default=all tracks]
            writetable(obj.getTracksTable(varargin{:}), fileName);
        end
    end % public methods

    methods %Dependant properties
        function ph = get.phase(obj)
            ph = obj.PHASES{obj.phaseIdx};
        end
        function N = get.nScales(obj)
            %Number of scales to filter on.
            if obj.phaseIdx < 2
                N = 0;
            else
                N = numel(obj.ParamsFindMaxima.filterSigmas);
            end
        end

        function N = get.nMaxima(obj)
            if obj.phaseIdx < 4
                N = 0;
            else
                N = size(obj.ResultsFilterMaxima.maxima,2);
            end
        end

        function N = get.nLocalizations(obj)
            if obj.phaseIdx < 6
                N = 0;
            else
                N = size(obj.ResultsFilterEmitters.emitters,1);
            end
        end

        function N = get.nTracks(obj)
            if obj.phaseIdx < 7
                N = 0;
            else
                N = length(obj.ResultsTrack.tracks);
            end
        end

        function bc = get.emitterBoxCoords(obj)
            if obj.phaseIdx < 4
                bc = [];
            else
                bc = obj.ResultsFilterMaxima.boxCoords;
            end
        end
        
        function bc = get.filteredEmitterBoxCoords(obj)
            if obj.phaseIdx < 6
                bc = [];
            else
                keepers = ~obj.ResultsFilterEmitters.filter;
                bc = BoxCoords.filter(obj.ResultsFilterMaxima.boxCoords, keepers);
            end
        end
    end % dependant properties

    methods (Abstract=true)
        % Localizations are sorted by frame
        % L - matrix format with column names given and each row a localization.  
        % col_names - Cell array of descriptions for each column of L.
        [L,col_names] = getLocalizations(obj)

        %% Computation functions
        % These are the actual computation functions where the real work is done
        findMaxima(obj);       %[phase 2->3]
        filterMaxima(obj);     %[phase 3->4]
        localizeEmitters(obj); %[phase 4->5]
        filterEmitters(obj);   %[phase 5->6]
        trackEmitters(obj);    %[phase 6->7]
    end %^public abstract methods

    methods (Abstract=true, Access=protected)
        loadData(obj, data, roi_in, roi_name);
        filterFrames(obj);
    end %prtoected abstract methods

    methods (Access=protected)
        function checkPhase(obj,minPhase)
            %Check we have met a minimum phase
            if obj.phaseIdx<minPhase
                error('RPT:InvalidPhase','RPT object must have phase "%s"(%i), but currently has phase "%s"(%i)',...
                        obj.PHASES{minPhase},minPhase,obj.PHASES{obj.phaseIdx},obj.phaseIdx);
            end
        end
        
        function setPhase(obj, newPhase)
            %Request the object's phase be reset to given value.
            if newPhase<1 || newPhase>length(obj.PHASES)
                error('RPT:setPhase','Invalid Phase:%i',newPhase)
            end
            if newPhase<3 && obj.phaseIdx>=3                
                obj.ResultsFindMaxima = structfun(@(~) [],obj.ResultsFindMaxima,'Uniform',0); %zero all fields
            end
            if newPhase<4 && obj.phaseIdx>=4                
                obj.ResultsFilterMaxima = structfun(@(~) [],obj.ResultsFilterMaxima,'Uniform',0);
            end
            if newPhase<4 && obj.phaseIdx>=4                
                obj.ResultsLocalizeEmitters = structfun(@(~) [],obj.ResultsLocalizeEmitters,'Uniform',0);
            end
            if newPhase<5 && obj.phaseIdx>=5                
                obj.ResultsFilterEmitters = structfun(@(~) [],obj.ResultsFilterEmitters,'Uniform',0);
            end
            if newPhase<6 && obj.phaseIdx>=6
                obj.srimage=[]; %zero out srimage
                obj.ResultsTrack = structfun(@(~) [],obj.ResultsTrack,'Uniform',0);
            end
            obj.phaseIdx=newPhase;
            obj.dirty=true;
        end

        function plotEventsPerFrame(obj, axH, frame_series, series_names)
            axes(axH);
            ts= 1:obj.nFrames;
            for i = 1:length(frame_series)
                [count, ~] = histc(frame_series{i},ts);
                plot(ts, count,'DisplayName',series_names{i});
                hold('on')
            end
            hold('off');
            legend('Location','East');
            yl=ylim();
            ylim([0 yl(2)]);
            xlabel('Frame index');
            ylabel('Count');
        end
    end % protected methods

    methods (Static = true)

        
        function tBounds = trackTimeBounds(Ts)
            % [in]
            %  Ts - A cell array of track matricies in RPT format.
            % [out]
            %  tBounds = [tmin, tmax] (sec)
            Ts = makecell(Ts);
            tBounds = [min(cellfun(@(T) T(1,1), Ts)), max(cellfun(@(T) T(end,1), Ts))];
        end

        function frameBounds = trackFrameBounds(Ts)
            % [in]
            %  Ts - A cell array of track matricies in RPT format.
            % [out]
            %  frameBounds = [frame_min, frame_max] (sec)
            Ts = makecell(Ts);
            frameBounds = [min(cellfun(@(T) T(1,end), Ts)), max(cellfun(@(T) T(end,end), Ts))];
        end

        function pBounds = trackPhysicalBounds(Ts, sigmaMultiple)
            % [in]
            %  Ts - A cell array of track matricies in RPT format.
            %  sigmaMultiple - How many sigmas out shoulf the bounds extend.  Default =3;
            % [out]
            %  pBounds = [minx, maxy, miny, maxy] (microns)
            if nargin<2
                sigmaMultiple = 3;
            end
            Ts = makecell(Ts);
            pBounds = [min(cellfun(@(T) min(T(:,2)-T(:,7)*sigmaMultiple), Ts)),...
                       max(cellfun(@(T) max(T(:,2)+T(:,7)*sigmaMultiple), Ts)),...
                       min(cellfun(@(T) min(T(:,3)-T(:,8)*sigmaMultiple), Ts)),...
                       max(cellfun(@(T) max(T(:,3)+T(:,8)*sigmaMultiple), Ts))];
        end
        %Static plotting methods for reuse by both 2D and HS data
        

        function S = makeLocalizationSurface(spheresize)
            % Make a surface object to show a localization as radially shaded circle at given pos and zdepth.
            % This returns all the coordinate matricies nessary for matlab surface() plotting
            % Those The sphere returned is a flat circle in the XY plane centered at [0,0,0] with unit radius and color 0
            % It the resulting arrays can then be shifted to show the localization in the appropriate location
            % [in]
            %   spheresize - the level of detail in the sphere rendering (size in polygons) [default=16]
            % [out]
            %   S - Struct with X,Y,Z,C,A component matricies all ready to be passed to sphere.
            if nargin<1
                spheresize=16;
            end
            [S.X,S.Y,~] = sphere(spheresize);
            S.X = S.X(spheresize/2+1:end,:); %cut out the bottom
            S.Y = S.Y(spheresize/2+1:end,:); %cut out the bottom
            S.Z = zeros(size(S.X)); % Make it flatand at given frame time
            S.C = zeros(size(S.X)); % Make space to use for colors
            S.A = 1-sqrt(S.X.^2 + S.Y.^2); %alpha based on distance from origin
        end

        function h = drawLocalizationSurface(pos, radius, zdepth, color, S)
            % draw a surface object to show a localization as radially shaded circle at given pos and zdepth.
            % draws into the current axes
            % [in]
            %   pos - the position [X Y]
            %   radius - the sphere radius [radX radY]
            %   zdepth - the zdepth of the sphere
            %   color - the numerical color value of the spere (for colormap mapping). This should be a scalar
            %   sphere - A struct of sphere coordinates as returned by makeLocalizationSurface.  This will allow us to reuse
            %            the already created datastructrue for efficiency.  If omitted we make a new localization
            %            surface sphere.
            % [out]
            %   h - Handle to surface object
            if nargin<5
                S = BaseRPT.makeLocalizationSurface();
            end
            h = surface('XData',S.X.*radius(1) + pos(1), 'YData',S.Y.*radius(2)+pos(2),...
                        'ZData',S.Z + zdepth,'CData',S.C+color,'FaceAlpha','texturemap',...
                        'EdgeColor','none','FaceColor','texturemap','CDataMapping','Direct', ...
                        'AlphaData',S.A,'AlphaDataMapping','scaled');
        end

        function h = drawTrackSurface(T, C, varargin)
            %
            % [in]
            %  T - track as a RPT track matric representation
            %  C - color values as array same length as track
            %  varargin - All other arguments passed to surface
            % [out]
            %  h - Handle to surface object
            if iscell(T) && iscell(C)
                h = cellfun(@(t,c) surface(t(:,[2,2]),t(:,[3,3]),t(:,[1,1]),[c,c],'LineWidth',2,'EdgeColor','interp','CDataMapping','Direct', varargin{:}), T, C);
            else
                h = surface(T(:,[2,2]),T(:,[3,3]),T(:,[1,1]),[C,C],'LineWidth',2,'EdgeColor','interp','CDataMapping','Direct', varargin{:});
            end
        end
        
        function h = configureTracksColorBar(nTrackColors, nImageColors, trackColorRange, label)
            % nTrackColors - Number of colors in the track part of the colormap
            % nImageColors - Number of colors in the image part of the colormap
            % trackColorRange - [minC, maxC] - min and max for the colors in all of the tracks'
            textColor = [1,1,1];
            h = colorbar();
            h.YLim = [0 nTrackColors]+1+nImageColors;
            h.Label.String = label;
            h.Label.Interpreter = 'tex';
            h.Label.FontSize = 12;
            h.Label.Color = textColor;
            h.Color = textColor;
            ticks = h.Ticks - nImageColors-1;
            if ~strcmp(label, 'Sequence')
                if ticks(1)~=0
                    ticks = [0 ticks];
                end
                if ticks(end)~= nTrackColors
                    ticks = [ticks nTrackColors];
                end
                ticks = ticks ./ nTrackColors;
                ticks = ticks * trackColorRange(2) + trackColorRange(1);
                ticks = round(ticks, 3, 'significant');
                h.Ticks = (ticks - trackColorRange(1)) .* (nTrackColors/trackColorRange(2)) + nImageColors +1;
                h.TickLabels = cellmap(@num2str,ticks);
            else
                ticks = unique(round(ticks));
                h.Ticks = (ticks - trackColorRange(1)) .* (nTrackColors/trackColorRange(2)) + nImageColors +1+0.5;
                h.TickLabels = cellmap(@num2str,ticks);        
            end
            
        end
        
        function fixROIAxesTicks(prop,minval,maxval)
            ticks=get(gca(),prop);
            sp=ticks(2)-ticks(1);
            if ticks(end)<maxval
                if maxval-ticks(end)>=0.75*sp;
                    ticks=[ticks,maxval];
                else
                    ticks(end)=maxval;
                end                    
            end
            if ticks(1)>minval
                if ticks(1)-minval >=0.75*sp;
                    ticks=[minval,ticks];
                else
                    ticks(1)=minval;
                end
            end
            set(gca,prop,ticks);
        end
        
        function nim = normalizeIndividualFrames(frames, nImageColors)
            %Determine color range
            nim = zeros(size(frames));
            for i=1:size(frames,3)                
                nim(:,:,i) = cosmicNorm(frames(:,:,i))*(nImageColors-1)+1; %Normalize each frame to [1,ImColRng+1)
            end        
        end
        
        function nim = normalizeIndividualFramesGlobal(frames, nImageColors)
            %Determine color range
            nim =  cosmicNorm(frames)*(nImageColors-1)+1; %Normalize each frame to [1,ImColRng+1)
        end
        
        function [normCs, minC, maxC] = normalizeTrackColors(Cs, nTrackColors, nImageColors)
            maxC = max(cellfun(@max, Cs));
            minC = min(cellfun(@min, Cs));
            if minC==maxC
                normCs = cellmap(@(c) ones(size(c))+nImageColors, Cs);
            else
                scaleF = nTrackColors/(maxC-minC);
                normCs = cellmap(@(c) (c-minC)*scaleF+nImageColors+1, Cs);
            end
        end
        
        function plotEventsPerScale(axH, scale_series, series_names)
            axes(axH);
            scales= 1:max(cellfun(@max, scale_series));
            for i = 1:length(scale_series)
                [count, ~] = histc(scale_series{i},scales);
                bar(scales, count,'DisplayName',series_names{i});
                hold('on')
            end
            hold('off');
            legend('Location','East');
            yl=ylim();
            ylim([0 yl(2)]);
            xlabel('Scale index');
            ylabel('Count');
        end

        function plot_imagesc(axH,im,origin)
            % This function plots a 2D image with imagesc
            axes(axH);
            imagesc([.5,size(im,2)-.5]+origin(1),[.5,size(im,1)-.5]+origin(2),im);
            xlabel('X (px)');
            axH.TickDir = 'out';
            axH.Box = 'on';
            axH.BoxStyle = 'full';
            ylabel('Y (px)');
            axis('tight');
        end

        function plot_surface(axH,im,origin)
            %plot image as a matlab surface
            axes(axH);
            sz = size(im);
            [X,Y] = meshgrid(double(0:sz(2))+origin(1),double(0:sz(1))+origin(2));
            % Expand image to be size+1 because of how surf plots faces
            sim = zeros(size(im)+1);
            sim(1:end-1,1:end-1) = im;
            surf('XData',X,'YData',Y,'ZData',sim,'CData',sim,'EdgeColor','none')
            axH.YDir = 'reverse';
            axH.TickDir = 'out';
            axH.Box = 'on';
            axH.BoxStyle = 'full';
            axis('tight');
            ratio = max(im(:))/max(sz);
            daspect([1 1 5*ratio]);
            xlabel('X (px)');
            ylabel('Y (px)');
            rotate3d('on');
            view(0,90);
        end

        function Hs = plot_distributions(axH,dists, names, normalized)
            % This function plots a 2D image with imagesc
            if nargin<4
                normalized = false;
            end
            axes(axH);
            minvals = cellfun(@(v) min(v(:)),dists);
            maxvals = cellfun(@(v) max(v(:)),dists);
            Npts = numel(dists{1});
            nbins = max(15,min(100,Npts^0.65));
            edges = linspace(min(minvals),max(maxvals),nbins);
            N = numel(dists);
            Hs = cell(N,1);
            for i = 1:N
                [count, ~] = histc(dists{i},edges);
                if normalized
                    count = count / sum(count);
                end
                Hs{i} =plot(edges,count,'DisplayName',names{i},'LineWidth',2);
                hold('on')
            end
            hold('off');
            legend('Location','NorthEast');
            xlim([min(minvals) max(maxvals)]);
            if normalized
                ylabel('Normalized Probability');
            else
                ylabel('Count');
            end
            Hs = [Hs{:}];
        end
        function Hs = plot_logDistributions(axH,dists, names, normalized)
            % This function plots a 2D image with imagesc
            if nargin<4
                normalized = false;
            end
            axes(axH);
            minvals = cellfun(@(v) min(v(:)),dists);
            maxvals = cellfun(@(v) max(v(:)),dists);
            Npts = numel(dists{1});
            nbins = max(15,min(100,Npts^0.65));
            edges = logspace(log10(min(minvals)),log10(max(maxvals)),nbins);
            N = numel(dists);
            Hs = cell(N,1);
            for i = 1:N
                [count, ~] = histc(dists{i},edges);
                if normalized
                    count = count / sum(count);
                end
                Hs{i} =plot(edges,count,'DisplayName',names{i},'LineWidth',2);
                hold('on')
            end
            hold('off');
            legend('Location','NorthEast');
            axH.XScale = 'log';
            xlim([min(minvals) max(maxvals)]);
            if normalized
                ylabel('Normalized Probability');
            else
                ylabel('Count');
            end
            Hs = [Hs{:}];
        end
    end % public static methods
end

