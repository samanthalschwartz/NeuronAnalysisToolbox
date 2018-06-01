
%still do to:------------
% clearing things from the lists
% button functionality for Hist and CumProbDist check
% Save file name 
%--------------

classdef HistPlotter < handle  
    
    %class properties - access is private so nothing else can access these
    %variables. Useful in different sitionations
    properties (Access = private)
        
%         GroupName = ''; % don't actually need this. just get group name
%         from gui as needed --- delete this
        GroupNames = {}; % Group names as defined by user: cell index should match with Groupings
        Groupings = {}; % n x 1 list of full file paths and names, n is number of groups, each index corresponds to a cell containing all files in group. index of list should match with GroupNames
        GroupTrackLengths; % list of tracklengths. this gets populated after Analyze is run. Indices are same as for GroupNames and Groupings. Needs to be cleared and loaded properly  
        defaultpathname; % this is just so that the uigetfile function defaults to starting with a useful directory
        savedir; % directory where all the analysis objects are to be saved. this also becomes the default for the plot groups uigetfile
         plot_title;
        gui_h;
    end
    
    properties (Access = public)
        fighist;
        cpdhist;
    end
    
    methods    
        %function - class constructor - creates and init's the gui
        function obj = HistPlotter()
            
            %make the gui handle and store it locally
            obj.gui_h = guihandles(PlotTrackLengthHists_analyzeplot()); % corresponding gui code for hist
            
            %set the Plot Button to be defaulted to off- this is so that
            %you can't try and plot anything unless there is something
            %loaded
            set(obj.gui_h.Plotbutton10,'Enable','off')
            
            
            %set the callback for re-setting everything
            set(obj.gui_h.ResetAllbutton8, 'callback', @(src, event) Reset_callback(obj, src, event)); % re-set all inputs to the gui
            
            %set the callback functions to the File Selection Functions box's
%             set(obj.gui_h.Filelistbox1, 'callback', @(src, event) Edit_Files(obj, src, event)); % List of Files: either to be part of a group or from a loaded group
            set(obj.gui_h.selectFilesbutton1, 'callback', @(src, event) obj.SelectFiles_callback(src, event)); % to select the files to go into Filelistbox1
            set(obj.gui_h.clearFilesbutton6, 'callback', @(src, event) obj.ClearFiles_callback(src, event)) % selected file(s) should be cleared
            
            %set the callback functions for the Grouping Functions
            set(obj.gui_h.GroupFilesbutton3, 'callback', @(src, event) GroupFiles_callback(obj, src, event)); % button to put all files into one group - uses GroupName and adds to GroupstoPlotlistbox3
            set(obj.gui_h.AddFilesIndividuallybutton5, 'callback', @(src, event) AddIndividualFiles_callback(obj, src, event)); % button to add individual files instead to GroupstoPlotlistbox3. uses filename
%             set(obj.gui_h.GroupNameedit1 , 'callback', @(src, event) %Edit_GroupName(obj, src, event)); %Edit the particular groupname maybe don't need this callback
%             set(obj.gui_h.GroupstoPlotlistbox3 , 'callback', @(src, event) Edit_Groupings(obj, src, event)); %List of Groups:
            set(obj.gui_h.saveGroupDiredit3 , 'callback', @(src, event) Edit_SaveDir(obj, src, event)); %Directory for saving things: default should be TrackLength subfolder in Results Directory
            set(obj.gui_h.selectSaveDirbutton13, 'callback', @(src, event) GetSaveDir_callback(obj, src, event)); % button to select a different save dir
       
            %set the callback function to run the analysis
            set(obj.gui_h.Analyzebutton4, 'callback', @(src, event) Analyze_callback(obj, src, event)); % uses groups in GroupstoPlotlistbox3 and runs analysis. Also saves groups in savedir folder.
                            % groups in GroupstoPlotlistbox3 should be transferred to GroupstoPlotlistbox3_CreateFcn after running 
            
            %set the callback functions for loading groups for plotting
            set(obj.gui_h.LoadGroupsbutton9, 'callback', @(src, event) LoadGroups_callback(obj, src, event)); % button to select different groups to use for plotting analysis
%             set(obj.gui_h.GroupstoPlotlistbox3, 'callback', @(src, event) Edit_Groups(obj, src, event)); % List of groups to be used for plotting comparison
            set(obj.gui_h.clearGroupsbutton11, 'callback', @(src, event) ClearGroup_callback(obj, src, event)); % button to clear the selected group out of GroupstoPlotlistbox3
            set(obj.gui_h.ShowFilesButton12, 'callback', @(src, event) ShowFiles_callback(obj, src, event)); % button to display the files that are part of the group in the Filelistbox1
            set(obj.gui_h.PlotTitleedit2 , 'callback', @(src, event) Edit_PlotTitle(obj, src, event)); % title of comparison plot, also figure saved with this name
           
            %set the callback function for plotting 
            set(obj.gui_h.Plotbutton10, 'callback', @(src, event) PlotGroups_callback(obj, src, event)); % button to clear the selected group out of GroupstoPlotlistbox3
         
            %Set the selection change Fcn for the radio button box. This
            %function will be called when the selection changes within the
            %box
%             set(obj.gui_h.unitgroup, 'selectionchangefcn', @(src, event) Ui_callback      (obj, src, event));
            
            %sets the figure close function. This lets the class know that
            %the figure wants to close and thus the class should cleanup in
            %memory as well
            set(obj.gui_h.figure1,  'closerequestfcn', @(src,event) Close_fcn(obj, src, event));
            %reset the gui (not needed, but  obj is used to duplicate the
            %functionality of the matlab example)
            obj = Reset(obj);
            
            
        end
        
    end
    
    
    %Private Class Methods - these functions can only be access by the
    %class itself.
    methods (Access = private)
        
        %class deconstructor - handles the cleaning up of the class &
        %figure. Either the class or the figure can initiate the closing
        %condition, obj function makes sure both are cleaned up
        function delete(obj)
            %remove the closerequestfcn from the figure, obj prevents an
            %infitie loop with the following delete command
            set(obj.gui_h.figure1,  'closerequestfcn', '');
            %delete the figure
            delete(obj.gui_h.figure1);
            %clear out the pointer to the figure - prevents memory leaks
            obj.gui_h = [];
        end
        
        %function - Close_fcn
        %
        %obj is the closerequestfcn of the figure. All it does here is
        %call the class delete function (presented above)
        function obj = Close_fcn(obj, src, event)
            delete(obj);
        end
        
        %function - Reset
        %
        %resets the gui to initial values. Called from the Reset_btn
        %callback and when the gui init's.
        %This function is mainly kept to mirror the functionality of the
        %MATLAB guide example
        %--------- Could just get rid of obj function or need to add buttons ------------------
        function obj = Reset(obj)    
            obj.GroupNames = {}; % Group names as defined by user: cell index should match with Groupings
            obj.Groupings = {}; % list of filenames, cell index of list should match with GroupNames
            obj.GroupTrackLengths = {};
            obj.defaultpathname='';
            obj.savedir='';
            set(obj.gui_h.Filelistbox1, 'String', '');
            set(obj.gui_h.GroupstoPlotlistbox3, 'String', obj.GroupNames);
            set(obj.gui_h.GroupNameedit1, 'String','GroupName');
            set(obj.gui_h.GroupstoPlotlistbox3,'String',obj.GroupNames); 
            set(obj.gui_h.saveGroupDiredit3,'String','');
            set(obj.gui_h.Plotbutton10,'Enable','off')
        end 
        
        %function - Reset_callback
        %
        %the callback function for the reset button. This simply calls
        %Reset function directly
        function obj = Reset_callback(obj, src, event)
            obj = Reset(obj); 
        end
        % --------------------------------                  -----------------
        
        
        % callback functions for selecting files and grouping 
%         function obj=Edit_Files(obj, src, event) %callback that does nothing currently- setup so matlab doesn't complain?
%             
%         end
        
        function obj = SelectFiles_callback(obj, src, event)
    if isempty(obj.defaultpathname)
        [filename, pathname, ~] = uigetfile('*.spt', 'MultiSelect', 'on');
        if ~pathname % if user didn't select a pathname then just return without doing anything
            return;
        end
        obj.defaultpathname=pathname;
    else
        [filename, pathname, ~] = uigetfile(fullfile(obj.defaultpathname,'*.spt'), 'MultiSelect', 'on');
        if ~pathname % if user didn't select a pathname then just return without doing anything
            return;
        end
        obj.defaultpathname=pathname;
        if ~strcmp(obj.defaultpathname,pathname) %prob faster to not check this and just change
            obj.defaultpathname=pathname;
        end
    end
    FileNames=get(obj.gui_h.Filelistbox1, 'String');
    if ischar(filename)
        FileNames=[FileNames; fullfile(pathname,{filename})]; %don't invert because only one filestring
    else
        FileNames=[FileNames; fullfile(pathname,filename)']; % otherwise invert so always displays properly in gui
    end
    
    set(obj.gui_h.Filelistbox1, 'String', FileNames); % add filenames to the list of files
    set(obj.gui_h.Filelistbox1,'Value',1);% this is to incorporate work around because selection doesn't default switch after you've deleted the file
    obj.savedir=fullfile(obj.defaultpathname,'TrackLengthGroups'); % set the object save dir field to be a subfolder of the used path.
    set(obj.gui_h.saveGroupDiredit3, 'String', obj.savedir); % display the save dir in the SaveDir text box.
    
end
        
        function obj= ClearFiles_callback(obj, src, event)
            filelist=get(obj.gui_h.Filelistbox1,'String'); %
            if isempty(filelist)
                return
            end
            index_toDelete = get(obj.gui_h.Filelistbox1,'Value'); % use get to obtain selected file index
            filelist(index_toDelete)=[]; % then delete that index from the filenames list.
            set(obj.gui_h.Filelistbox1, 'String', filelist); % update the filelistbox
            set(obj.gui_h.Filelistbox1,'Value',1); % this is to incorporate work around because selection doesn't default switch after you've deleted the file
            
        end
        
        function obj = GroupFiles_callback(obj, src, event) % call addGroup class method to make it a cell of cell array
        % callback functions for grouping files, adding files individually,
        % and setting the save directory
            groupname = get(obj.gui_h.GroupNameedit1, 'String'); % get group name
            if ~groupname
                display('Need to give this group a Name')
            else
            FileNames=get(obj.gui_h.Filelistbox1, 'String'); % get list of files
            % check that name isn't already in the list
            if sum(strcmp(groupname,obj.GroupNames))>0
                groupname=''; % if it is then set to empty string and continue to prompt user until they have actually added something
                while strcmp(groupname,'')
                groupname=inputdlg('Oops, you already have a group of this name. Please pick a valid new name.','Group Name',1);
                end
            end
            
            obj.addGroup(groupname,{FileNames})  % use class function to add group and all files
            
            %now clean up
            set(obj.gui_h.GroupstoPlotlistbox3, 'String', obj.GroupNames);
            set(obj.gui_h.Filelistbox1, 'String', '');
            set(obj.gui_h.GroupNameedit1, 'String','GroupName')
             set(obj.gui_h.GroupNameedit1, 'Value',1)
            set(obj.gui_h.Filelistbox1,'Value',1);
             set(obj.gui_h.Plotbutton10,'Enable','off') 
            end
        end
        
        function obj = AddIndividualFiles_callback(obj, src, event) % call addGroup class method as loop to add as cell(1 file) of cell array 
            % get filenames from queue
            % still doesn't check to make sure file name isn't already in
            % the list- needs to be fixed
            FileNames=get(obj.gui_h.Filelistbox1, 'String');
            if isempty(FileNames) % double check that there actually are files there - list box will be a cell - but an empty one.
                return;
            end
                           
                [PATHSTR,NAME,EXT]= cellfun(@(x) fileparts(x), FileNames, 'UniformOutput',false);
                % now call class method
                for ii=1:numel(NAME)
                 obj.addGroup(NAME(ii),{FileNames(ii)}); %input to addGroup is a cell - so that always a cell and never have to worry about chars.
                end
              % now clean up gui 
                set(obj.gui_h.GroupstoPlotlistbox3, 'String', obj.GroupNames);
                set(obj.gui_h.Filelistbox1, 'String', '');
                set(obj.gui_h.GroupNameedit1, 'String','GroupName')
                set(obj.gui_h.Filelistbox1,'Value',1);
                 set(obj.gui_h.Plotbutton10,'Enable','off') 
            end
       

        function obj = addGroup(obj, varargin) %-- should work that groupname and filestoadd are more than one thing for single files being passed in.
           % method of class that adds group to obj.GroupNames and obj.Groupings
           % groupname: string of groupname to add
           % filestoadd: cell array containig list of files to add. should
           % be ** dimension
           obj.GroupNames = [obj.GroupNames; varargin{1}]; % Add file name to Group names: cell index should match with Groupings
           obj.Groupings = [obj.Groupings; varargin{2}]; % list of filenames, cell index of list should match with GroupNames
           if size(varargin,2)==3 % also has grouptracklength info need to make sure this is appended correctly - cell?
               obj.GroupTrackLengths = [obj.GroupTrackLengths; {varargin{3}}] ;
           end
        end

        function obj=Edit_SaveDir(obj, src, event)
            % callback function for user to manually edit save directory (note directory only created when actually saving files)
            %
            savedir=get(obj.gui_h.saveGroupDiredit3,'String');
            obj.savedir=savedir;
        end
        
        function obj=GetSaveDir_callback(obj, src, event) 
            % callback function for user button to select the save directory (note directory only created when actually saving files)
            %
            if ~exist(obj.savedir) %check if there is already a save dir, default to start with this
                savedir=uigetdir('','Pick a Directory to Save Groups')
            else
                savedir=uigetdir(fullfile(obj.savedir),'Pick a Directory to Save Groups');
            end
            obj.savedir=savedir; % save this directory as save dir field of object
            set(obj.gui_h.saveGroupDiredit3, obj.savedir); % display the save dir in the SaveDir text box.
        end

        function obj = Analyze_callback(obj, src, event) % clears values from filegroupqueue at end
            % clears out GroupTrackLengths var
            % obj.GroupTrackLengts
            % Callback for analysing groups and saving results!
            %check to make sure this works without a cell?
            % first create empty 1xn cell array for holding all the data tracklengths. n is the number of groups
            
            % tracklengths are multiplied by frame rate so they are
            % reported in seconds!!
            obj.GroupTrackLengths={}; % variable for all the tracklengths of the group
            % loop through each group name
            for gettr=1:length(obj.GroupNames)
                % for each group name, loop through each file
                GroupTrackLengths=[];
                grplen=numel(obj.Groupings{gettr}); %get length of groupings index to determine if single file.
                for jj=1:grplen
                     % otherwise index into the group properly
                        sptObj = SPT.loadFile(obj.Groupings{gettr}{jj}); % load in .spt file
                   
                    tracklenvalscell=arrayfun(@(x) max(x.Frame)-min(x.Frame), sptObj.Tracks,'UniformOutput',false); % use array fun to get through all tracks
                    tracklenvals=cell2mat(tracklenvalscell); % make cell array into matrix
                    GroupTrackLengths=cat(2,GroupTrackLengths,tracklenvals*sptObj.ParamsGeneral.TimeStep); % concat matrix into cell array of all data at appropriate group index
                end
                obj.GroupTrackLengths{gettr}=GroupTrackLengths;
                
                
                % create save dir here if it doesn't exist
                if ~exist(obj.savedir)
                    mkdir(obj.savedir)
                end
                grouptracklengths=obj.GroupTrackLengths{gettr};
                groupfiles=obj.Groupings{gettr}; % also save the groupfiles- list of full file path and names for all files that made up tracklengths
                groupname=obj.GroupNames{gettr}; % also save the groupname to ensure that tracklengths are always linked to correct groupname
                % save the concatenated tracks together here and the groupname as a variable into a .mat file
                save(fullfile(obj.savedir,obj.GroupNames{gettr}),'grouptracklengths','groupfiles','groupname');
                % note that there is no modification to object here.
                % obj.GroupNames
            end
            % double check that obj.GroupTrackLengths is same size as obj.GroupNames
            
            % now automatically populate groups into Plot Group Comparison list for plotting. User can always delete later.
            set(obj.gui_h.GroupstoPlotlistbox3,'String', obj.GroupNames); % this shouldn't change anything so can just delete
            % now allow user to use the plot button command
            set(obj.gui_h.Plotbutton10,'Enable','on') ;
            
        end
        
%         function obj = Edit_Groupings(obj, src, event) % function for
%         making sure we can ... actually don't need- should delte
        function obj = LoadGroups_callback(obj, src, event)
    if numel(obj.Groupings)~=numel(obj.GroupTrackLengths)
        errordlg('There are un-analyzed files in the queue, please delete or analyze them before adding any new groups', 'Warning');
        return
    end
    if isempty(obj.defaultpathname)
        [filename, pathname, ~] = uigetfile('*.mat', 'MultiSelect', 'on');
        obj.defaultpathname=pathname;
    else
        [filename, pathname, ~] = uigetfile(fullfile(obj.defaultpathname,'*.mat'), 'MultiSelect', 'on');
        if ~strcmp(obj.defaultpathname,pathname) %prob faster to not check this and just change
            obj.defaultpathname=pathname;
        end
    end
    % now load in each group and add it to the obj.Group fields using
    % addGroup class function
    if ischar(filename)
        load(fullfile(pathname,filename))
        obj.addGroup(groupname,{groupfiles},grouptracklengths)
    else
        for ii=1:numel(filename)
            load(fullfile(pathname,filename{ii})); %- load the file
            % for each file load in, add to the obj. Group properties
            obj.addGroup(groupname,{groupfiles},grouptracklengths)
        end
    end
    
    set(obj.gui_h.GroupstoPlotlistbox3,'String', obj.GroupNames) % update the list of groups
    
    set(obj.gui_h.Plotbutton10,'Enable','on');  %set plot buttong
end

        % could make delete group button and call back
        function obj = ClearGroup_callback(obj,src,event)
            % function to clear a group from group queue and from object
            % list
            
            % first make sure there are actually groups listed - then clear
            % them
            if isempty(get(obj.gui_h.GroupstoPlotlistbox3,'String')) && isempty(obj.GroupNames)
                return
            end
            %get selection
           index_toDelete = get(obj.gui_h.GroupstoPlotlistbox3,'Value'); % use get to obtain selected file index
           obj.GroupNames(index_toDelete)=[];
           obj.Groupings(index_toDelete)=[];
           if ~isempty(obj.GroupTrackLengths)
               obj.GroupTrackLengths(index_toDelete)=[];
           end
           set(obj.gui_h.GroupstoPlotlistbox3, 'String', obj.GroupNames); % update the groupqueue list box
           set(obj.gui_h.Filelistbox1,'Value',1); % this is to incorporate work around because selection doesn't default switch after you've deleted the file
            set(obj.gui_h.GroupstoPlotlistbox3,'Value',1);
            if isempty(get(obj.gui_h.GroupstoPlotlistbox3, 'String'))
                set(obj.gui_h.Plotbutton10,'Enable','off') 
            end
            
            % empty group
            %update everything
        end
        
        function obj = PlotGroups_callback(obj, src, event) %function for plotting what's in the groupstoplot queue -- need to fix saving!!!
             
           %--------------------- code to take for plotting 
            % make histogram plot
            tracklentots=cell(1,length(obj.GroupNames));
            clrs=lines(length(obj.GroupNames));
            edges=0:.1:50; %----------------------max histogram readout here is fixed to 1000 frames: should change?-----------------------
            fhist=figure('Color',[1 1 1]);
            axhist=axes('Parent',fhist,'FontSize',16,'LineWidth',2);
            for gettr=1:length(obj.GroupNames) % make histogram
                [n,bin]=histc(obj.GroupTrackLengths{gettr}(obj.GroupTrackLengths{gettr}>0),edges);
                yval=n./sum(n);
                yval(yval==0)=NaN;
                plot(edges,yval,'Color',clrs(gettr,:),'LineWidth',2); hold on;
            end
            legnames=cellfun(@(x) strrep(x,'_','-'),obj.GroupNames,'UniformOutput',false);
            xlabel('Tracklength (sec)','FontSize',16,'FontWeight','bold');
            ylabel('Fraction of Total Tracks','FontSize',16,'FontWeight','bold');
            legend(legnames,'FontSize',10,'FontWeight','bold','Location', 'SouthEast');
            title([obj.plot_title ' Histogram'],'FontSize',16,'FontWeight','bold'); 
%             obj.fighist=fhist; don't save this as part of object delete
       
            % make CPA of values 
            fcpa=figure('Color',[1 1 1]);
            axcpa=axes('Parent',fcpa,'FontSize',16,'LineWidth',2);
            cpatable=zeros(length(obj.GroupTrackLengths),2);
            for ii=1:length(obj.GroupTrackLengths) %make CPD plot
                tracklens=obj.GroupTrackLengths{ii}(obj.GroupTrackLengths{ii}>0);
%                 tracklens=tracklens(tracklens>4); % change this to change which values are included
                sortracklens=sort(tracklens,2,'ascend');
                cpavals=(0:size(sortracklens,2)-1)./size(sortracklens,2);
                semilogx(sortracklens,cpavals,'Color',clrs(ii,:),'LineWidth',2,'Parent',axcpa); hold on;
                cpatable(ii,1)=sortracklens(find((cpavals-0.5)>0,1,'first'));
                cpatable(ii,2)=sortracklens(find((cpavals-0.8)>0,1,'first'));
            end
            xlabel('Tracklength (sec)','FontSize',16,'FontWeight','bold');
            ylabel('Cumulative Probability Distribution','FontSize',16,'FontWeight','bold');
            legend(legnames,'FontSize',10,'FontWeight','bold','Location', 'SouthEast');
            title([obj.plot_title ' CPD'],'FontSize',16,'FontWeight','bold'); 
%              obj.cpdhist=fcpa; don't save as part of object delete
            % make a table of results
             colnames={'50% Value (sec)','80% Value (sec)'};
            tablefig=figure('Color',[1 1 1]);             
            test=uitable('Parent',tablefig,'Data',cpatable,'ColumnName',colnames,'RowName',legnames');
             tbpos=get(test,'Position');
             tbext=get(test,'Extent');
             set(test,'Position',[tbpos(1:2) tbext(3:4)]);
             savedir=obj.savedir;    %fullfile(obj.defaultpathname,'Track Length Plots');
             if ~exist(savedir,'dir')
                 mkdir(savedir)
             end
             saveas(fhist,fullfile(savedir,[obj.plot_title ' hist']),'fig');
             saveas(fhist,fullfile(savedir,[obj.plot_title ' hist']),'png');
             
             saveas(fcpa,fullfile(savedir,[obj.plot_title ' CPD']),'fig');
             saveas(fcpa,fullfile(savedir,[obj.plot_title ' CPD']),'png');
             
             saveas(tablefig,fullfile(savedir,[obj.plot_title ' table']),'fig');
             saveas(tablefig,fullfile(savedir,[obj.plot_title ' table']),'png');
             
        end   
        
        function obj = Edit_PlotTitle(obj, src, event)
            %read in the value from the edit box.
            ptitle = get(obj.gui_h.PlotTitleedit2, 'String');
            obj.plot_title=ptitle;
        end
        
        function obj = ShowFiles_callback(obj, src, event)
            index_toShowFiles = get(obj.gui_h.GroupstoPlotlistbox3,'Value'); % use get to obtain selected file index
            set(obj.gui_h.Filelistbox1, 'String', obj.Groupings{index_toShowFiles}); % update the filelistbox to be the list of files
        end
        
    end
    
end

%bugs to fix:
% uncomment save function
% error if click clear file when there aren't any files left




%--- delete
%                 obj.GroupNames = [obj.GroupNames; NAME]; % Add file name to Group names: cell index should match with Groupings
% %                 names=cellfun(@(x) {x}, obj.FileNames','UniformOutput',false); % add in full file name for the Groupings so that the can be loaded later. keep as cells so indexing works
%                 
%                obj.Groupings = [obj.Groupings; ]; % list of filenames, cell index of list should match with GroupNames
                %--- delete
%                         FileNames = {};

%  obj.GroupNames = [obj.GroupNames; groupname]; % Group names as defined by user: cell index should match with Groupings
%            if ~isempty(obj.Groupings)
%            obj.Groupings = [obj.Groupings; {obj.FileNames}]; % list of filenames, cell index of list should match with GroupNames
%            else
%               obj.Groupings={obj.FileNames}; 
%            end