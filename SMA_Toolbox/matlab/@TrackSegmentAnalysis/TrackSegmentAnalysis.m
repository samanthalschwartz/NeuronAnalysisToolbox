%Mark J. Olah (mjo@cs.unm.edu)
% 01-2015
%
%
% A Class for track segment analysis


classdef TrackSegmentAnalysis < Pickle & GUIBuilder

    properties (Constant = true)
        %Inherited from Pickle
        saveFileExt='.tsa'; % The extension used for saved files
        SaveableDataFormats={'*.tsa', 'TrackSegementAnalysis (.tsa)'};
        LoadableDataFormats={'*.tsa;*.rpt;*.spdata;','All Loadable Sources (.tsa, .rpt, .spdata)';...
                             '*.tsa','TrackSegmentAnalysis file (.tsa)';'*.spdata','SPData file'; '*.rpt', 'RPT object files (.rpt)'};

        CSVFormatsLinux={'*.csv', 'Comma Seperated Value File (*.csv)'};
        CSVFormatsWin={'*.csv;*.xls','All exporable formats'; '*.csv', 'Comma Seperated Value File (*.csv)'; '*.xls','MicrosoftExcel File (.xls)'};
        MinDmleTrackLength=10; % The shortest track we will attempt to estimate a D_mle value for
        MinTrackLength=4; % The minimum track length to load in. Maybe remove this?
        %These define the ordering of table outputs
        TrackStatsFieldOrder={'NSegments', 'Classification','TotalTime','NumberLocalizations',...
                              'TotalDistance', 'NetDistance', 'MaxDistance','ConfinementRatio','Dmle',...
                              'MeanSpeed', 'NetSpeed', 'LinearityRatio',...
                              'InstSpeed','MeanVelocity','InstVelocity'};
        SegStatsFieldOrder={ 'Classification','TotalTime','NumberLocalizations',...
                              'TotalDistance', 'NetDistance', 'MaxDistance','ConfinementRatio','Dmle',...
                              'MeanSpeed', 'NetSpeed', 'LinearityRatio',...
                              'InstSpeed','MeanVelocity','InstVelocity'};
        LocStatsFieldOrder={'Displacements', 'TotalCumDistances', 'NetCumDistances', 'Speed','Direction','Angle'};

        % Table Properties for obj.trackStatsTable() track statistics table
        TrackStatsTableColumns={'TrackID','NSegments', 'Classification','TotalTime','NumberLocalizations',...
                                'TotalDistance', 'NetDistance', 'MaxDistance','ConfinementRatio','Dmle',...
                                'MeanSpeed', 'NetSpeed', 'LinearityRatio'};
        TrackStatsTableUnits={'', '', '', 's', '', 'um', 'um','um','','um^2/s','um/s','um/s',''};
        % Table Properties for obj.segmentStatsTable() segment statistics table
        SegStatsTableColumns={'TrackID','SegmentID', 'Classification','TotalTime','NumberLocalizations',...
                              'TotalDistance', 'NetDistance', 'MaxDistance','ConfinementRatio','Dmle',...
                              'MeanSpeed', 'NetSpeed', 'LinearityRatio'};
        SegStatsTableUnits={'', '', '', 's', '', 'um', 'um','um','','um^2/s','um/s','um/s',''};
        % Table Properties for obj.trackLocalizationsTable() track localization statistics table
        LocStatsTableColumns={'TrackID','SegmentID','LocalizationID','FrameIdx',...
                              't','x','y','I','bg','sigma',...
                              'Displacement','TotalCumDist','NetCumDist','Speed','Direction','Angle'};
        LocStatsTableUnits={'', '', '', '','s', 'um', 'um','photons','photons/px','',...
                            'um','um','um','um/s','deg','deg'};

    end % public constant properties

    properties 
        ROI; % The ROI in pixels [xmin xmax ymin ymax]
        ROIname;
        ROIPhysical; % The ROI in real coordinate units [xmin xmax ymin ymax]
        pixelSize;
        frameT;

        Paths=struct(... %All filenames are relative paths from the workingDir
             'saveFile',[],... %filename(w/extension) of the .tsa file where this will save to
             'data',[],... %relative path to the .spdata file we are associated with
             'RPT',[],...
             'SPT',[],...
             'Export','ExportedTracks'...
            );

        % A Track Segment is stored as a cell array of [nT x nTColumn] track matrices where each matrix
        %  represents a single track.  Cols are given by TColumns and related constants. 
        %  Rows are sequential localizations.  This format is inherited
        %  from RPT.
        % 't'         col 1: time in seconds
        % 'x'         col 2: x-position in micron
        % 'y'         col 3: y-position in micron
        % 'I'         col 4: wavelength in phonton count
        % 'bg'        col 5: background intenisty in photons/pixel in the ROI box
        % 'sigma'     col 6: Measured fit sigma (um)
        % 'SE_x'     col 7:  variance for x-position in (um)
        % 'SE_y'     col 8:  variance for y-position in (um)
        % 'SE_I'     col 9: variance for intensity in counts
        % 'SE_bg'    col 10: variance  for background parameter in counts
        % 'SE_sigma' col 11: variance for the sigma parameter (um)
        % 'frame'     col 12: frame idx, the index of the frame this
        %                     localization  corresoponds to
        Tracks; % cell(ntrack,1) -> cell(nseg,1) -> RPT track format matrix;
        trackStats;% struct(ntrack,1)
        segStats;% cell(ntrack,1) -> struct(nseg,1)
    end %set-restricted public properties
    
    properties (Dependent=true)
        nTracks; %Number of tracks
    end

    properties (Transient=true)
        data;        
    end

    properties (Hidden=true)
        version=3; %For future file format version changes
    end    
    
    methods
        %% Track Statistics Methods %%

        function obj=TrackSegmentAnalysis( varargin )
            % This constructor calls the obj.load method with the arguments
            % given.  Empty input arguments give an empty object, otherwise the load
            % method determines the valid input arguments
            if nargin>0
                obj.load(varargin{:});
            end
        end
        
        function load(obj, varargin)
            %
            % [Case1]: Load from saved .tsa file
            %    obj.load(tsaFilePath)
            %  tsaFilePath - Full path to a .tsa file to load
            %
            % [Case2]: New TSA from SPData and RPT/SPT file
            %    obj.load(spd,roi,spt)
            %
            %  spd - An SPData filepath or object
            %  roi - An roi index into the SPData object's list of ROI,
            %        set to empty for whole area.
            %  spt/rpt - [optional] An SPT/RPT filepath or object with tracks already computed.
            %        If omitted automatically look for .rpt then .spt file when use the first save .spt file in the ROI path.
            %
            % [Case3]: New TSA from RPT
            %  obj.load(rpt)
            %  
            %  rpt - An RPT object or full path to file
            % 
            in=varargin{1};
            if ischar(in)
                [~,~,ext]=fileparts(in);
                switch ext
                    case '.tsa'
                        obj.loadTSA(in);                        
                    case '.spdata'
                        obj.loadSPData(varargin{:});
                    case '.rpt'
                        obj.loadRPT(in);
                    otherwise
                        error('TrackSegmentAnalysis:load','Unknown file format extenstion "%s"',ext);
                end
            elseif isa(in,'SPData')
                obj.loadSPData(varargin{:});
            elseif isa(in,'RPT')
                obj.loadRPT(in);
            else
                error('TrackSegmentAnalysis:load','Unknown argument type "%s"',class(in));
            end
        end
        
%         function s = getSegment(obj, trackid, segmentids)
%             % segmentids - If empty or omitted, get all segments for this
%             % track id as a cell array.  If single index, get just that
%             % asegment as a matrix of localizations.   If multiple ids,
%             % then return selected ids from this trackid as a cell array.
%             T = obj.Tracks{trackid};
%             if nargin==2 || isempty(segmentids)
%                 s = T; % Get all segments as a cell array
%             elseif length(segmentids)==1
%                 s = T{segmentids}; %Get this one segment as a matrix of localizations
%             else
%                 s = T(segmentids); %Get selected ids as a cell array
%             end
%         end
        
        function [Ts, col_names, col_units, col_desc] = getTracks(obj, trackIdxs)
            % Get tracks with all segments concatenated in RPT format.
            % See RPT for more info on tracks formats.
            % [in]
            %   trackIdxs - [optional] An array of track indexes. [default=all tracks]
            % [out]
            %   Ts - Tracks as cellarray of matricies in RPT format
            %   col_names - [optional] the value of RPT.TrackColumnNames for convenience
            %   col_units - [optional] the value of RPT.TrackColumnUnits for convenience
            %   col_desc - [optional] the value of RPT.TrackColumnDescriptions for convenience
            if nargin<2
                trackIdxs = 1:obj.nTracks;
            end
            Ts = cellmap(@(t) obj.spliceSegments(t), obj.Tracks(trackIdxs));
            if nargout>1 % Provide these for convenience only
                col_names = RPT.TrackColumnNames;
                col_units = RPT.TrackColumnUnits;
                col_desc = RPT.TrackColumnDescriptions;
            end
        end
        
        function tableT = getTracksTable(obj, varargin)
            % This is a reformatting of the tracks format into a single table where all tracks have
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
            tableT = RPT.convertTracksTable(obj.getTracks(varargin{:}), varargin{:});
        end

        function [structT, field_units, field_desc] = getTracksStruct(obj, varargin)
            % This is a reformatting of the tracks format into a structure array where each element of the
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
            structT = RPT.convertTracksStruct(obj.getTracks(varargin{:}));
            if nargout>1 % Provide these for convenience only
                field_units = RPT.TrackColumnUnits;
                field_desc = RPT.TrackColumnDescriptions;
            end
        end

        function n = nSegments(obj, trackid)
            n = length(obj.Tracks{trackid});
        end

        function n = trackLength(obj, trackid)
            n = sum(cellfun(@(T) size(T,1) ,obj.Tracks{trackid}));
        end
        
        function n = segLength(obj, trackid, segid)
            n = size(obj.Tracks{trackid}{segid},1);
        end
        
        function deleteTrack(obj, trackids)
            if any(trackids<1 | trackids>obj.nTracks)
                error('TrackSegmentAnalysis:divideSegment','Bad trackids: %s ',arr2str(trackid));
            end
            obj.Tracks(trackids)=[];
            obj.trackStats(trackids)=[];
            obj.segStats(trackids)=[];
        end

        function trackIds = divideTrack(obj, trackid, segid, locid)
            %First track will include localization rows 
            % 1:splitidx  next track will include splitidx+1:end
            if trackid<1 || trackid>obj.nTracks
                error('TrackSegmentAnalysis:divideTrack','Bad trackid: %i ',trackid);
            elseif segid<1 || segid>obj.nSegments(trackid)
                error('TrackSegmentAnalysis:divideTrack','Bad segmentid:%i for trackid: %i ',segid,trackid);
            elseif locid<=1 || locid>=obj.segLength(trackid,segid);
                error('TrackSegmentAnalysis:divideTrack','locid: %i must be internal to segment length: %i',...
                       locid,obj.segLength(trackid,segid));
            end
            segs = obj.Tracks{trackid};
            stats = obj.segStats{trackid};
            %Construct new segment and statistics
            s1 = segs{segid}(1:locid,:);
            s2 = segs{segid}(locid:end,:);
            s1_stats = obj.computeSegStats(s1);
            s1_stats.Classification = stats(segid).Classification;
            s2_stats = obj.computeSegStats(s2);
            s2_stats.Classification = stats(segid).Classification;
            T1 = [segs(1:segid-1), s1];
            T2 = [s2, segs(segid+1:end)];
            T1_stats = obj.computeTrackStats(T1);
            T2_stats = obj.computeTrackStats(T2);
            old_seg_stats = obj.segStats{trackid};
            T1_seg_stats = [old_seg_stats(1:segid-1), s1_stats];
            T2_seg_stats = [s2_stats, old_seg_stats(segid+1:end)];
            % Modify the Datastructures
            obj.Tracks = [obj.Tracks(1:trackid-1), {T1, T2}, obj.Tracks(trackid+1:end)];
            obj.segStats = [obj.segStats(1:trackid-1), T1_seg_stats, T2_seg_stats, obj.segStats(trackid+1:end)];
            obj.trackStats = [obj.trackStats(1:trackid-1), T1_stats, T2_stats, obj.trackStats(trackid+1:end)];
            trackIds=[trackid,trackid+1];
        end
%         
%         function joinTrack(obj, trackids)
%             % trackids must be valid sequential (2 or more) indexes of segments to merge.
%             if length(trackids)<2 || any( trackids<1 | trackids>obj.nTracks )
%                 error('TrackSegmentAnalysis:joinTracks','trackids invalid');
%             end
%             for i = 1:length(trackids)-1
%                 if obj.Tracks{trackids(i)}{end}(end:end) >= obj.Tracks{trackids(i+1)}{1}(1:end)
%                     error('TrackSegmentAnalysis:joinTracks','Tracks overlap in time. Cannot join.');
%                 end
%             end
%             
%             %Construct new segment and statistics
%             new_track = [obj.Tracks{trackids}];            
%             for i = length(trackids):-1:1 %Merge the first and last segment of each track
%                 idx=trackids(i);
%                 new_track = [new_track(1:idx-1), cell2mat(new_track(idx:idx+1)), new_track(idx+2:end)];
% %                 classifications = [
%             end
% %             new_seg_stats = cellmap(@(s) obj.computeSegmentStats(s), new_track);
% %             for i = 1:length(trackids)
% %                 n=obj.nSegments(trackids(i));                
% %                 new_seg_stats(
%             new_track_stats = obj.computeTrackStats(new_track);
%             new_stats.Classification = stats{segids(1)}.Classification; %Copy classification
%             % Modify the Datastructures
%             obj.Tracks{trackid} = [segs(1:segids(1)-1), new_seg, segs(segids(end)+1:end)];
%             obj.segStats{trackid} = [stats(1:segids(1)-1), new_stats, stats(segids(end)+1:end)];
%             obj.trackStats{trackid}.NSegments = obj.nSegments(trackid);
%         end
         
        function divideSegment(obj, trackid, segid, locid)
            %locid - will be last point of previous segment and first
            %point on new segment
            if trackid<1 || trackid>obj.nTracks
                error('TrackSegmentAnalysis:divideSegment','Bad trackid: %i ',trackid);
            elseif segid<1 || segid>obj.nSegments(trackid)
                error('TrackSegmentAnalysis:divideSegment','Bad segmentid:%i for trackid: %i ',segid,trackid);
            elseif locid<=1 || locid>=obj.segLength(trackid,segid);
                error('TrackSegmentAnalysis:divideSegment','locid: %i must be internal to segment length: %i',...
                       locid,obj.segLength(trackid,segid));
            end
            segs = obj.Tracks{trackid};
            stats = obj.segStats{trackid};
            %Construct new segment and statistics
            s1 = segs{segid}(1:locid,:);
            s2 = segs{segid}(locid:end,:);
            s1_stats = obj.computeSegStats(s1);
            s1_stats.Classification = stats(segid).Classification;
            s2_stats = obj.computeSegStats(s2);
            s2_stats.Classification = stats(segid).Classification;
            % Modify the Datastructures
            obj.Tracks{trackid} = [segs(1:segid-1), s1, s2, segs(segid+1:end)];
            obj.segStats{trackid} = [stats(1:segid-1), s1_stats, s2_stats, stats(segid+1:end)];
            obj.trackStats(trackid).NSegments = obj.nSegments(trackid);
        end
        
        function joinSegment(obj, trackid, segids)
            % segids must be valid sequential (2 or more) indexes of segments to merge.
            segids=segids(:)';
            if trackid<1 || trackid>obj.nTracks
                error('TrackSegmentAnalysis:joinSegment','Bad trackid: %i ',trackid);
            elseif length(segids)<2 ...
               || length(segids(1):segids(end))~=length(segids)...
               || ~all( (segids == round(segids)) & (segids == sort(segids)) & (segids == segids(1):segids(end)))
                error('TrackSegmentAnalysis:joinSegment','segids:%s must be 2 or more contiguous segment ids',mat2str(segids));
            elseif segids(1)<1 || segids(end)>obj.nSegments(trackid);
                error('TrackSegmentAnalysis:joinSegment','segids:%s invalid for track:%i with %i segment(s)',mat2str(segids),trackid,obj.nSegments(trackid));
            end
            segs = obj.Tracks{trackid};
            stats = obj.segStats{trackid};
            %Construct new segment and statistics
            new_seg = obj.spliceSegments(segs(segids));
            new_stats = obj.computeSegStats(new_seg);
            new_stats.Classification = stats(segids(1)).Classification; %Copy classification
            % Modify the Datastructures
            obj.Tracks{trackid} = [segs(1:segids(1)-1), new_seg, segs(segids(end)+1:end)];
            obj.segStats{trackid} = [stats(1:segids(1)-1), new_stats, stats(segids(end)+1:end)];
            obj.trackStats(trackid).NSegments = obj.nSegments(trackid);
        end
        
        %% Track Statistics Retrival and Outputs

        function recomputeAllStats(obj)
            %This can be called to force a re-evaluation of all stats for
            %all tracks and segments

            %Preserve Classifications
%             if ~isempty(obj.trackStats)
%                 trackClassifications = {obj.trackStats(:).Classification};
%                 segmentClassifications = cellmap(@(t) {t(:).Classification}, obj.segStats);
%             else
%                 trackClassifications = [];
%                 segmentClassifications = [];
%             end
            %Recompute stats
            obj.trackStats = cellfun(@(T) obj.computeTrackStats(T), obj.Tracks);
            obj.segStats = cellmap(@(T) arrayfun(@(i) obj.computeSegStats(T{i}), 1:length(T)), obj.Tracks);
            %Restore classifications
%             if ~isempty(trackClassifications)
%                 [obj.trackStats(:).Classification] = trackClassifications{:};                 
%                  obj.segStats = cellmap(@segStatsMap, obj.segStats, segmentClassifications);
%             end
%             function T = segStatsMap(T,C)
%                 [T(:).Classification] = C{:};
%             end
        end
        
        function stats=computeTrackStats(obj,track)
            % track is a cellarray of segments
            if iscell(track)
                NSegments = length(track);
                track = TrackSegmentAnalysis.spliceSegments(track);
            else
                NSegments = 1;
            end
            stats = obj.computeTrackSummaryStats(track, obj.frameT);
            stats.Classification = 'None';
            stats.NSegments = NSegments;
            stats = orderfields(stats,obj.TrackStatsFieldOrder);
        end
        
        function stats = computeSegStats(obj, segment)
            %segment is a matrix in standerd RPT track format
            stats = obj.computeTrackSummaryStats(segment, obj.frameT);
            stats.Classification = 'None';
            stats = orderfields(stats,obj.SegStatsFieldOrder);
        end
        
        function tab = trackStatsTable(obj,trackIds)
            % [IN]
            %  trackIds [optional] - Array of indexes into the Tracks data structure to select.
            %                        if omitted selects all tracks
            % [OUT]
            %  tab - A table of track stats
            if nargin==1
                stats_arr = obj.trackStats;
                trackIds = 1:obj.nTracks;
            else
                stats_arr = obj.trackStats(trackIds);
            end
            if isempty(stats_arr)
                tab=table();
                return
            end
            stats_arr = structfilter(stats_arr,obj.TrackStatsTableColumns); %Select only the relevent columns
            tab = [array2table(trackIds(:)), struct2table(stats_arr(:))];
            tab.Properties.VariableNames = obj.TrackStatsTableColumns;
            tab.Properties.VariableUnits = obj.TrackStatsTableUnits;
            tab.Properties.Description = 'Track Statistics';
            tab.Properties.DimensionNames = {'Track', 'Statistics'};
        end
 
        function tab = segmentStatsTable(obj,varargin)
            % [IN]
            %  trackIds [optional] - Array of indexes into the Tracks data structure to select.
            %                        If omitted selects all segments of all tracks are selected.
            %  segmentIds [optional] - Array of same length as trackIds that lists the corresponding
            %                        segments to select.  If more than one segment from an track is to be
            %                        selected then there must be a seperate trackIds and segIds entry for
            %                        each segmetent.
            %                        If omitted then all segments of selected tracks are selected.
            % [OUT]
            %  tab - A table of segment stats
            [trackIds, segIds ] = obj.getTrackSegmentInputIds(varargin{:});
            stats_arr = arrayfun(@(t,s) obj.segStats{t}(s), trackIds, segIds);
            stats_arr = structfilter(stats_arr,obj.SegStatsTableColumns); %Select only the relevent columns
            tab = [array2table([trackIds(:),segIds(:)]), struct2table(stats_arr(:))];
            tab.Properties.VariableNames = obj.SegStatsTableColumns;
            tab.Properties.VariableUnits = obj.SegStatsTableUnits;
            tab.Properties.Description = 'Segment Statistics';
            tab.Properties.DimensionNames = {'Segment', 'Statistics'};
        end
        
        
        
        function tab=trackLocalizationsTable(obj,varargin)
            % [IN]
            %  trackIds [optional] - Array of indexes into the Tracks data structure to select.
            %                        If omitted selects all segments of all tracks are selected.
            %  segmentIds [optional] - Array of same length as trackIds that lists the corresponding
            %                        segments to select.  If more than one segment from an track is to be
            %                        selected then there must be a seperate trackIds and segIds entry for
            %                        each segmetent.
            %                        If omitted then all segments of selected tracks are selected.
            % [OUT]
            %  tab - A table of localizations stats
            [trackIds, segIds ] = obj.getTrackSegmentInputIds(varargin{:});
            function seg_arr=prepareLocStats(track_id, seg_id)
                seg=obj.Tracks{track_id}{seg_id};
                n = size(seg,1);
                track_idx = ones(n,1)*track_id;
                seg_idx = ones(n,1)*seg_id;
                loc_idx = (1:n)';
                cum_stats = orderfields(obj.computeTrackCumulativeStats(seg), obj.LocStatsFieldOrder);
                seg_arr=[track_idx, seg_idx, loc_idx, seg(:,12), seg(:,1:6), struct2array(cum_stats)];
            end
            stats_arr = cellmatfun(@(i) prepareLocStats(trackIds(i),segIds(i))', 1:numel(trackIds) )';
            tab = array2table(stats_arr);
            tab.Properties.VariableNames = obj.LocStatsTableColumns;
            tab.Properties.VariableUnits = obj.LocStatsTableUnits;
            tab.Properties.Description = 'Track Segment Localization Statistics';
            tab.Properties.DimensionNames = {'Localization', 'Statistics'};
        end

        function exportAllStats(obj, fileName)
            % One output file with all of the track, segment and localization level statistics in
            % a mulittab output file
            if isunix()
                error('TrackSegmentAnalysis:exportAllStats','Not availible on linux.');
            end
            [~,~,ext] = fileparts(fileName);
            if ~strcmp(ext,'.xls')
                error('TrackSegmentAnalysis:exportAllStats','Can only write All Stats to .xls file');
            end
            warning('off');
            writetable(obj.trackStatsTable(), fileName,'Sheet','Tracks');
            writetable(obj.segmentStatsTable(), fileName,'Sheet','Segments');
            writetable(obj.trackLocalizationsTable(), fileName,'Sheet','Localizations');
            warning('on');
        end

        function exportTrackStats(obj, fileName, varargin)
            tab = obj.trackStatsTable(varargin{:});
            writetable(tab, fileName);
        end
        
        function exportSegmentStats(obj, fileName, varargin)
            tab = obj.segmentStatsTable(varargin{:});
            writetable(tab, fileName);
        end
        
        function exportLocalizationStats(obj, fileName, varargin)
            tab = obj.trackLocalizationsTable(varargin{:});
            writetable(tab, fileName);
        end
        
        function [sqDistCdf, sqDistXs]=getSquaredDisplacementCDFOverall(obj, varargin)
            % [IN]
            %  trackIds [optional] - Array of indexes into the Tracks data structure to select.
            %                        If omitted selects all segments of all tracks are selected.
            %  segIds [optional] - Array of same length as trackIds that lists the corresponding
            %                        segments to select.  If more than one segment from an track is to be
            %                        selected then there must be a seperate trackIds and segIds entry for
            %                        each segmetent.
            %                        If omitted then all segments of selected tracks are selected.             
            [trackIds, segIds ] = obj.getTrackSegmentInputIds(varargin{:});
            [sqDistCdf, sqDistXs] = squaredDisplacementCDF(cellmap(@(t,s) obj.Tracks{t}{s}, trackIds, segIds));
        end
        
        function [sqDistCdf, sqDistXs]=getSquaredDisplacementCDFIndividual(obj, varargin)
            % [IN]
            %  trackIds [optional] - Array of indexes into the Tracks data structure to select.
            %                        If omitted selects all segments of all tracks are selected.
            %  segIds [optional] - Array of same length as trackIds that lists the corresponding
            %                        segments to select.  If more than one segment from an track is to be
            %                        selected then there must be a seperate trackIds and segIds entry for
            %                        each segmetent.
            %                        If omitted then all segments of selected tracks are selected. 
            [trackIds, segIds ] = obj.getTrackSegmentInputIds(varargin{:});
            [sqDistCdf, sqDistXs] = cellmap(@(t,s) squaredDisplacementCDF(obj.Tracks{t}{s}), trackIds, segIds);
        end
        
        function [msds, lags]=getMSDTemporal(obj, varargin)
            % [IN]
            %  trackIds [optional] - Array of indexes into the Tracks data structure to select.
            %                        If omitted selects all segments of all tracks are selected.
            %  segIds [optional] - Array of same length as trackIds that lists the corresponding
            %                        segments to select.  If more than one segment from an track is to be
            %                        selected then there must be a seperate trackIds and segIds entry for
            %                        each segmetent.
            %                        If omitted then all segments of selected tracks are selected. 
            
            [trackIds, segIds ] = obj.getTrackSegmentInputIds(varargin{:});
            track_arr = cellmap(@(t,s) obj.Tracks{t}{s}, trackIds, segIds);
            [msds, lags] = cellmap(@(T) msdtemporal(T), track_arr);
        end
        
        %% Track Visualization
        function h=view3DTrackSequence(obj,varargin)
            h = figure();
            obj.plot3DTrackSequence(varargin{:})            
            h.Position = h.Position + [-200 -200 400 200];
        end
        
        function h=view3DTrackTemporal(obj,varargin)
            h = figure();
            obj.plot3DTrackTemporal(varargin{:})
            h.Position = h.Position + [-200 -200 400 200];
        end
        
        function h = view3DTrackSpeed(obj,varargin)
            h = figure();
            obj.plot3DTrackSpeed(varargin{:})
            h.Position = h.Position + [-200 -200 400 200];
        end
        
        function h=view3DSegmentSequence(obj,varargin)
            h = figure();
            obj.plot3DSegmentSequence(varargin{:})            
            h.Position = h.Position + [-200 -200 400 200];
        end
        
        function h=view3DSegmentTemporal(obj,varargin)
            h = figure();
            obj.plot3DSegmentTemporal(varargin{:})
            h.Position = h.Position + [-200 -200 400 200];
        end
        
        function h = view3DSegmentSpeed(obj,varargin)
            h = figure();
            obj.plot3DSegmentSpeed(varargin{:})
            h.Position = h.Position + [-200 -200 400 200];
        end

        function h = view3DTrackMovie(obj, trackIds, varargin)
            if nargin==1 || isempty(trackIds)
                trackIds = 1:obj.nTracks;
            end
            tm = obj.makeTrackMovie(trackIds);
            h=tm.viewSequence(varargin{:});
        end
        
        function track_movie = makeTrackMovie(obj, trackIds)
            if nargin==1 || isempty(trackIds)
                trackIds = 1:obj.nTracks;
            end
            Ts = [obj.Tracks{trackIds}];
            track_movie = TrackMovie(Ts, obj.data.getFrames(obj.ROI), obj.ROIPhysical(1:4));
        end
        
        function plot3DTrackSequence(obj,trackIds, opts)
            %  trackIds [optional] - Array of indexes into the Tracks data structure to select.
            %                        If omitted selects all segments of all tracks are selected.
            %  opts [optional] - If provided this must be 2rd argument, so you must provide trackIds also.  
            %     It is a strcut which allow for setting of several options.
            %  
            if nargin==1 || isempty(trackIds)
                trackIds = 1:obj.nTracks;
            end
            Ts = [obj.Tracks{trackIds}];            
            Cs = cellmap(@(i) i*ones(size(Ts{i},1),1), 1:length(Ts));
            stats = table2struct(obj.segmentStatsTable(trackIds));           
            opts.cLabel = 'TrackID';
            opts.trackColorMap = @prism;
            opts.trackColorRange = length(trackIds);
            obj.plotTrackSegmentStats(Ts,Cs,stats,opts);
        end

        function plot3DTrackTemporal(obj,trackIds,opts)
            if nargin==1 || isempty(trackIds)
                trackIds = 1:obj.nTracks;
            end
            Ts = [obj.Tracks{trackIds}];
            Cs = cellmap(@(T) T(:,1), Ts);
            stats = table2struct(obj.segmentStatsTable(trackIds));           
            opts.cLabel= 'Time (s)';
            obj.plotTrackSegmentStats(Ts,Cs,stats,opts);
        end
        
        
        function plot3DTrackSpeed(obj,trackIds,opts)
            if nargin==1 || isempty(trackIds)
                trackIds = 1:obj.nTracks;
            end
            if nargin<=2 || ~isfield(opts,'winsize')
                opts.winsize= 6;
            end
            Ts = [obj.Tracks{trackIds}];
            Cs = cellmap(@(T) obj.estimateSpeedWindow(T, opts.winsize), Ts);
            stats = table2struct(obj.segmentStatsTable(trackIds));           
            opts.cLabel = 'Speed ($\mu\mathrm{m}/\mathrm{s}$)';
            obj.plotTrackSegmentStats(Ts,Cs,stats,opts);
        end
        
        function plot3DSegmentSequence(obj,varargin)
            % [IN]
            %  trackIds [optional] - Array of indexes into the Tracks data structure to select.
            %                        If omitted selects all segments of all tracks are selected.
            %  segIds [optional] - Array of same length as trackIds that lists the corresponding
            %                        segments to select.  If more than one segment from an track is to be
            %                        selected then there must be a seperate trackIds and segIds entry for
            %                        each segmetent.
            %                        If omitted then all segments of selected tracks are selected. 
            %  opts [optional] - If provided this must be 3rd argument, so you must provide trackIds and segIds
            %  also.  It is a strcut which allow for setting of several options.
            %  
            if numel(varargin)==3
                opts = varargin{3};
            else
                opts = struct();
            end
            [trackIds, segIds] = obj.getTrackSegmentInputIds(varargin{:});            
            Ts = cellmap(@(t,s) obj.Tracks{t}{s},trackIds,segIds);            
            Cs = cellmap(@(i) i*ones(size(Ts{i},1),1), 1:length(Ts));
            stats = table2struct(obj.segmentStatsTable(trackIds, segIds));           
            opts.cLabel = [];
            if ~isfield(opts,'trackColorMap')
                opts.trackColorMap = @prism;
            end
            opts.trackColorRange = length(Ts);
            obj.plotTrackSegmentStats(Ts,Cs,stats,opts);
        end
        
        function plot3DSegmentSpeed(obj,varargin)
            % [IN]
            %  trackIds [optional] - Array of indexes into the Tracks data structure to select.
            %                        If omitted selects all segments of all tracks are selected.
            %  segIds [optional] - Array of same length as trackIds that lists the corresponding
            %                        segments to select.  If more than one segment from an track is to be
            %                        selected then there must be a seperate trackIds and segIds entry for
            %                        each segmetent.
            %                        If omitted then all segments of selected tracks are selected.
            %  opts [optional] - If provided this must be 3rd argument, so you must provide trackIds and segIds
            %  also.  It is a strcut which allow for setting of several options.
            if numel(varargin)==3
                opts = varargin{3};
            else
                opts = struct();
            end
            if ~isfield(opts,'winsize')
                opts.winsize=6;
            end
            [trackIds, segIds] = obj.getTrackSegmentInputIds(varargin{:});            
            Ts = cellmap(@(t,s) obj.Tracks{t}{s},trackIds,segIds);            
            Cs = cellmap(@(T) obj.estimateSpeedWindow(T, opts.winsize), Ts);
            stats = table2struct(obj.segmentStatsTable(trackIds, segIds));           
            opts.cLabel = 'Speed ($\mu\mathrm{m}/\mathrm{s}$)';
            obj.plotTrackSegmentStats(Ts,Cs,stats,opts);
        end
        
        function plot3DSegmentTemporal(obj,varargin)
            % [IN]
            %  trackIds [optional] - Array of indexes into the Tracks data structure to select.
            %                        If omitted selects all segments of all tracks are selected.
            %  segIds [optional] - Array of same length as trackIds that lists the corresponding
            %                        segments to select.  If more than one segment from an track is to be
            %                        selected then there must be a seperate trackIds and segIds entry for
            %                        each segmetent.
            %                        If omitted then all segments of selected tracks are selected. 
            %  opts [optional] - If provided this must be 3rd argument, so you must provide trackIds and segIds
            %  also.  It is a strcut which allow for setting of several options.
            if numel(varargin)==3
                opts = varargin{3};
            else
                opts = struct();
            end
            [trackIds, segIds] = obj.getTrackSegmentInputIds(varargin{:});            
            Ts = cellmap(@(t,s) obj.Tracks{t}{s},trackIds,segIds);            
            Cs = cellmap(@(T) T(:,1), Ts);
            stats = table2struct(obj.segmentStatsTable(trackIds, segIds));           
            opts.cLabel = 'Time (s)';
            obj.plotTrackSegmentStats(Ts,Cs,stats,opts);
        end
        
        function h=viewSquaredDisplacementCDFOverall(obj,varargin)
            % [IN]
            %  trackIds [optional] - Array of indexes into the Tracks data structure to select.
            %                        If omitted selects all segments of all tracks are selected.
            %  segIds [optional] - Array of same length as trackIds that lists the corresponding
            %                        segments to select.  If more than one segment from an track is to be
            %                        selected then there must be a seperate trackIds and segIds entry for
            %                        each segmetent.
            %                        If omitted then all segments of selected tracks are selected. 
            [F, xs] = obj.getSquaredDisplacementCDFOverall(varargin{:});
            h = figure();
            if isempty(varargin)
                name = 'Combined Squared Displacement CDF (AllTracks)';
            elseif numel(varargin{1})==1
                name = sprintf('Squared Displacement CDF (Track: %i)',varargin{1}(1));
            else
                name = 'Combined Squared Displacement CDF (Selected Tracks)';
            end 
            semilogx(xs,F,'DisplayName',name);
            xlim([xs(1), xs(end)]);
            title('Overall Squared Displacement CDF','interpreter','latex');
            legend('location','best');
            xlabel('Squared Displacement ($\mu \mathrm{m}^2$)','interpreter','latex');
            ylabel('Cumulative Probability','interpreter','latex');
        end
        function h=viewSquaredDisplacementCDFIndividual(obj,varargin)
            % [IN]
            %  trackIds [optional] - Array of indexes into the Tracks data structure to select.
            %                        If omitted selects all segments of all tracks are selected.
            %  segIds [optional] - Array of same length as trackIds that lists the corresponding
            %                        segments to select.  If more than one segment from an track is to be
            %                        selected then there must be a seperate trackIds and segIds entry for
            %                        each segmetent.
            %                        If omitted then all segments of selected tracks are selected. 
            [F, xs] = obj.getSquaredDisplacementCDFIndividual(varargin{:});
            h = figure();
            ax = axes();
            hold('on');
            [trackids, segids] = obj.getTrackSegmentInputIds(varargin{:}); 
            for n = 1:numel(F)
                name = sprintf('Track %i Segment %i',trackids(n), segids(n));
                plot(xs{n},F{n},'DisplayName',name);                
            end
            hold('off');
            ax.XScale='log';
            xlim([min(cellfun(@min,xs)), max(cellfun(@max,xs))]);
            title('Individual Squared Displacement CDF','interpreter','latex');
            legend('location','best');
            xlabel('Squared Displacement ($\mu \mathrm{m}^2$)','interpreter','latex');
            ylabel('Cumulative Probability','interpreter','latex');
        end
        function h=viewMSDTemporal(obj,varargin)
            % [IN]
            %  trackIds [optional] - Array of indexes into the Tracks data structure to select.
            %                        If omitted selects all segments of all tracks are selected.
            %  segIds [optional] - Array of same length as trackIds that lists the corresponding
            %                        segments to select.  If more than one segment from an track is to be
            %                        selected then there must be a seperate trackIds and segIds entry for
            %                        each segmetent.
            %                        If omitted then all segments of selected tracks are selected. 
            [msds, lags] = obj.getMSDTemporal(varargin{:});
            h = figure();
            ax = axes();
            hold('on');
            [trackids, segids] = obj.getTrackSegmentInputIds(varargin{:}); 
            for n = 1:numel(msds)
                name = sprintf('Track %i Segment %i',trackids(n), segids(n));
                plot(lags{n},msds{n},'DisplayName',name);                
            end
            minlag = min(cellfun(@min, lags));
            maxlag = max(cellfun(@max, lags));
            startMSD = mean(cellfun(@(m) mean(m(1:2)) ,msds));
            lin_xs = logspace(log10(minlag), log10(maxlag), 10);
            plot(lin_xs, startMSD/minlag * lin_xs,'k:','DisplayName','(linear reference line)');
%             xlim([minlag,maxlag]);
            hold('off');
            ax.XScale='log';
            ax.YScale='log';
            title('MSD vs Time Lag','interpreter','latex');
            legend('location','best');
            ylabel('MSD($\tau$) ($\mu \mathrm{m}^2$)','interpreter','latex');
            xlabel('Time lag $\tau$ (s)','interpreter','latex');
        end
    end %public methods
    
    %% Dependant properties
    methods
        function n=get.nTracks(obj)
            n=length(obj.Tracks);
        end
    end %dependant properties accessors

    

    methods (Static=true)
        function obj=loadobj(propstruct)
            %Called by the matlab loading routine to initialize an object.  Anything that
            %needs to be re-initialized on load should go here
            obj = TrackSegmentAnalysis();
            switch propstruct.version
                case 1
                    propstruct=TrackSegmentAnalysis.updateVersion1(propstruct);
            end
            obj.reloadobj(propstruct);
        end
        
        function [speed, direction] = estimateSpeedWindow(T, wsize)
            if wsize<2
                error('TrackSegementAnalysis:estimateSpeedWindow','window size must be >=2');
            end
            N = size(T,1); %number of localiztions
            speed=zeros(N,1);
            direction=zeros(N,1);
            if N<=wsize
                disp = T(end,2:3)-T(1,2:3);
                speed(1:N) = sqrt(sum(disp.^2,2)) / (T(end,1)-T(1,1));
                direction(1:N) = (180/pi)*atan2(disp(:,2),disp(:,1));
                return
            end
            ts = T(:,1);
            ps = T(:,[2,3]); % 2D points
            S = floor(wsize/2)+1; % The start index for the window
            E = ceil(wsize/2)-1; % The end index for the window
            disp = ps(1:end-wsize+1,:) - ps(wsize:end,:);
            dt = ts(wsize:end,:) - ts(1:end-wsize+1,:); 

            
            speed(S:end-E) = sqrt(sum(disp.^2,2))./dt;            
            direction(S:end-E) = (180/pi)*atan2(disp(:,2),disp(:,1)); %degrees

            %Extend speed and direction to the same lenght as T, by copying first and last elements.
            speed(1:S-1) = speed(S);
            direction(1:S-1) = direction(S);
            
            speed(end-E+1:end) = speed(end-E);
            direction(end-E+1:end) = direction(end-E);            
        end

        function stats = computeTrackSummaryStats(T, frameT)
            % [IN]
            % T - A Track in RPT format (All distances in microns, times in seconds)
            % frameT - The frame exposure time in seconds.  (Used to compute motion-blur diffusion
            %             esimtation effects)
            % [OUT]
            % stats - structure with lots of statistics regarding the track as a whole
            N = size(T,1); %number of localiztions
            ts = T(:,1);
            ps = T(:,[2,3]); % 2D points
            steps = diff(ps); %displacments
            dts = diff(ts);
            dists = sqrt(sum(steps.^2,2)); %pointwise distances;
            vs = steps./repmat(dts,1,2); %stepwise velocities.
            speeds = dists./dts; %stepwise speeds
          
            stats.TotalTime = ts(end)-ts(1);
            stats.NumberLocalizations = N;
            if stats.NumberLocalizations >= TrackSegmentAnalysis.MinDmleTrackLength
                destr = DEstimator(ps,ts,T(:,[6,7]),frameT);
                stats.Dmle = destr.MLE();
            else
                stats.Dmle = 0;
            end
            if stats.NumberLocalizations==1
                stats.TotalDistance = 0;
                stats.NetDistance = 0; 
                stats.MaxDistance = 0;
                stats.ConfinementRatio = 0;
                stats.MeanSpeed = 0;
                stats.NetSpeed = 0;
                stats.MeanVelocity = 0;            
                stats.LinearityRatio = 0;
                stats.InstSpeed = 0;
                stats.InstVelocity = 0;
            else
                stats.TotalDistance = sum(dists);
                stats.NetDistance = sqrt(sum((ps(end,:)-ps(1,:)).^2)); %end to end
                stats.MaxDistance = max(arrayfun(@(i) max(sqrt(sum((ps-repmat(ps(i,:),N,1)).^2, 2)))  , 1:N));
                stats.ConfinementRatio = stats.NetDistance/stats.TotalDistance;
                stats.MeanSpeed = mean(speeds);
                stats.NetSpeed = stats.NetDistance/stats.TotalTime;
                stats.MeanVelocity = mean(vs);            
                stats.LinearityRatio = stats.NetSpeed/stats.MeanSpeed; %1=linear
                stats.InstSpeed = speeds;
                stats.InstVelocity = vs;
            end
        end

        function stats = computeTrackCumulativeStats(T)
            % Compute statistics that are relevent at the individual localization level
            N = size(T,1); %number of localiztions
            ps = T(:,[2,3]); % 2D points
            steps = diff(ps); %displacments
            dists = sqrt(sum(steps.^2,2)); %pointwise distances;            
          
            if N==1
                stats.Displacements = 0;
            else
                stats.Displacements = [0; dists];
            end
            stats.TotalCumDistances = cumsum(stats.Displacements);
            stats.NetCumDistances = sqrt(sum((ps(:,:)-repmat(ps(1,:),N,1)).^2,2)); % to end

            windowSize = 6;
            stats.Angle=zeros(N,1);
            if N>=6
                [stats.Speed, stats.Direction] = TrackSegmentAnalysis.estimateSpeedWindow(T, windowSize);
                stats.Angle(2:end)=mod(diff(stats.Direction),360)-180;
            else
                stats.Speed=zeros(N,1);
                stats.Direction=zeros(N,1);
            end
        end        
    end % Public static methods

    methods (Access=protected)        
        function loadTSA(obj, tsafile)
            % Load from a fullpath to a .tsa file
            [wDir,saveFileBase,~] = fileparts(tsafile); 
            s = load(tsafile,'-mat');
            obj.copyobj(s.obj);
            obj.workingDir = collapsepath(wDir);
            obj.Paths.saveFile = [saveFileBase obj.saveFileExt];
            obj.data = SPData(fullfile(obj.workingDir,obj.Paths.data));
            obj.initialized = true;
            obj.dirty = false;
        end

        function loadRPT(obj, rpt)
            % Load from  RPT object 
            if isa(rpt,'RPT')
                obj.data = rpt.data;
            elseif ischar(rpt)
                rpt = RPT(rpt);
                obj.data = rpt.data;
            else
                error('TrackSegmentAnalysis:loadRPT','Bad RPT argument.  Must be a .rpt filename or RPT object.');
            end
            obj.ROI = rpt.ROI;
            obj.ROIname = rpt.ROIname;
            obj.workingDir = obj.data.getFilePath('TrackSegmentAnalysis');
            obj.Paths.RPT = relativepath(obj.workingDir, rpt.saveFilePath);
            obj.initializeNewTracks(rpt.getTracks())
        end

        function loadSPData(obj, spd, roi_in, trackobj)
            % Load from  SPData and SPT objects given as files or objects.
            % trackobj -  can be an SPT or RPT filepath or object
            if isa(spd,'SPData')
                obj.data = spd;
            elseif ischar(spd)
                obj.data = SPData(spd);
            else
                error('TrackSegmentAnalysis:loadSPData','Bad SPData argument.  Must be a .spdata filename or SPData object.');
            end
            obj.ROI = obj.data.getROI(roi_in);
            obj.workingDir = obj.data.getFilePath('TrackSegmentAnalysis');
            obj.frameT = obj.data.frameT;
            obj.pixelSize = [1 1] .* obj.data.pixelSize;
            %Determine ROI name
            if isscalar(roi_in) && roi_in>0 && roi_in<=length(obj.data.ROI)  %Index into predefined ROI             
                obj.ROIname = obj.data.ROIname{roi_in};
            elseif nargin==2 || isempty(roi_in) %Full frame ROI
                obj.ROIname = 'FullFrame';
            else
                error('TrackSegmentAnalysis:loadSPData','Unable to determine ROI name');
            end
            
            %First try to find an rpt or spt from this dataset
            rpt=[];
            spt=[];
            if nargin<4
                trackobj = obj.findRPTFile();
                if isempty(trackobj)
                    trackobj = obj.findSPTFile();
                    if isempty(trackobj)
                        error('TrackSegmentAnalysis:loadSPData','Unable to find .rpt or .spt file for "%s_%s"',obj.data.saveFileBaseName,obj.ROIname);
                    end
                    obj.Paths.SPT = relativepath(obj.workingDir, trackobj);
                    spt = SPT(trackobj);                 
                else
                    rpt = RPT(trackobj);                    
                end
            elseif isa(trackobj,'RPT')
                rpt = trackobj;
            elseif isa(trackobj,'SPT')
                spt = trackobj;
            elseif ischar(trackobj)
                [~,~,ext] = fileparts(trackobj);
                switch ext
                    case '.spt'
                        obj.Paths.SPT = relativepath(obj.workingDir, trackobj);
                        spt = SPT(trackobj);
                    case '.rpt'
                        rpt = RPT(trackobj); 
                end
            end
            
            % Now do the loading 
            if ~isempty(rpt)
                obj.Paths.RPT = relativepath(obj.workingDir, rpt.saveFilePath);
                obj.initializeNewTracks(rpt.getTracks())
            elseif ~isempty(spt)
                tracks = sptTracksToCellArray(spt,obj.frameT,obj.pixelSize);
                obj.initializeNewTracks(tracks);
            else
                error('TrackSegmentAnalysis:loadSPData','Unable to open or get valid .rpt or .spt');
            end
        end
        
        function initializeNewTracks(obj, tracks)
            % 
            %Common code for both RPT and SPT initializations
            % [IN] tracks - RPT tracks format as cell array of matiricies
            %Remove short tracks
            
            [~, file_name, file_ext] = obj.data.getROIFileNameParts(obj.ROIname, 'TrackSegmentAnalysis');
            save_pattern = sprintf('%s*%s',file_name, file_ext);
            obj.Paths.saveFile = relativepath(obj.workingDir,Pickle.findUnusedFileName(obj.workingDir,save_pattern));
            obj.Paths.data = relativepath(obj.workingDir, obj.data.saveFilePath);

            obj.pixelSize = obj.data.pixelSize .* [1 1];
            obj.frameT = obj.data.frameT;       
            
            obj.ROIPhysical = [(obj.ROI(1:2)-1)*obj.pixelSize(1), obj.ROI(3:4)*obj.pixelSize(2),...
                               (obj.ROI(5)-1)*obj.frameT, obj.ROI(6)*obj.frameT];
            nTracks = find(cellfun(@length,tracks)>=obj.MinTrackLength, 1, 'last');
            obj.Tracks = num2cell(tracks(1:nTracks)); %Convert to single segment tracks   
            %compute stats
            obj.recomputeAllStats();
            obj.initialized = true;
            obj.dirty = true;
        end

        function spt_path = findSPTFile(obj)
            %Look automatically for .spt 
            file_pattern = sprintf('%s_%s*.spt',obj.data.saveFileBaseName, obj.ROIname);
            spt_files = Pickle.listExistingFileNames(obj.data.getFilePath('SPT'),file_pattern);
            if isempty(spt_files)
                spt_path = [];
            elseif length(spt_files)==1
                spt_path = spt_files{1};
            else
                formats = {file_pattern, ['SPT Files for: ' obj.ROIname]};
                spt_path = Pickle.selectExistingFileName(obj.data.getFilePath('SPT'),file_pattern, formats,...
                                                'Load SPT File to retrieve tracks from');
            end
        end
        
        function rpt_path=findRPTFile(obj)
            %Look automatically for .spt
            [rpt_path, file_name, file_ext] = obj.data.getROIFileNameParts(obj.ROIname, 'RPT');
            save_pattern = sprintf('%s*%s',file_name, file_ext);
            rpt_files = Pickle.listExistingFileNames(rpt_path, save_pattern);
            if isempty(rpt_files)
                rpt_path = [];
                return;
            elseif length(rpt_files)==1
                rpt_path = rpt_files{1};
            else
               formats = {file_pattern, ['RPT Files for: ' obj.ROIname]};
               rpt_path = Pickle.selectExistingFileName(obj.data.getFilePath('RPT'),file_pattern, formats,...
                                                'Load RPT File to retrieve tracks from');
            end
        end
        
        
        
        %% Plotting helpers
        function plotTrackSegmentStats(obj, Ts,Cs,Stats,opts)
            % Makes a viewer that uses dataCursorMode object to allow stats viewing 
            im = double(obj.data.getSumImage(obj.ROI));

            trackHs = RPT.plotTracks3D(Ts, Cs, obj.ROIPhysical,  im, opts);
            %Set the track user data for the data cursor to use.  Consists of a tuple of the stats and
            %the track itself.
            dcm_obj = datacursormode(gcf());
            dcm_obj.UpdateFcn = @obj.cursorStats;
            dcm_obj.DisplayStyle = 'window';
            dcm_obj.Enable = 'on';
            for i = 1:length(trackHs)
                if isa(trackHs(i),'matlab.graphics.primitive.Surface');
                    trackHs(i).UserData = {Stats(i),Ts{i}};
                end
                if isfield(opts,'ButtonDownFcn')
                    trackHs(i).ButtonDownFcn = opts.ButtonDownFcn;
                end
            end
        end

        
        %% Abstract methods inherited from Pickle
        function val = getProtectedProperty(obj, name)
            %This is necessary for Pickle functionality to be able to access subclass protected variables
            val = obj.(name);
        end

        function modifyProtectedProperty(obj, name, newval)
            %This is necessary for Pickel functionality to be able to change subclass variables
            obj.(name) = newval;
        end

    end % protected methods

    methods (Static=true, Access=protected)        
        function newS = spliceSegments(segs)
            % Given an cell array of segments concatenate them, removing
            % the repeated localizations at the segment boundaries.            
            newS = cell2mat([cellmatfun(@(s) s(1:end-1,:), segs(1:end-1)'); segs(end)]);
        end
    
        function msg = cursorStats(~,event)
                % This is a GUI callback for the data cursor mode UpdateFcn callback
                pos = event.Position;
                
                v = event.Target.UserData;
                if isempty(v)
                    msg = 'Select a Track or Segment';
                    return
                end
                stats=v{1};
                ts=v{2};
                idx = binarysearch(ts(:,1),pos(3));
                msg={sprintf('TrackID: %i',stats.TrackID),...
                     sprintf('SegID: %i',stats.SegmentID),...
                     sprintf('LocID: %i',idx),...
                     sprintf('X: %.5g um',pos(1)),...
                     sprintf('Y: %.5g um',pos(2)),...
                     sprintf('T: %.5g s',pos(3)),...
                     sprintf('Classification: %s',stats.Classification),...
                     sprintf('TotalTime: %.3g s',stats.TotalTime),...
                     sprintf('#Localizations: %i',stats.NumberLocalizations),...
                     sprintf('TotalDistance: %.3g um',stats.TotalDistance),...
                     sprintf('NetDistance: %.3g um',stats.NetDistance),...
                     sprintf('MaxDistance: %.3g um',stats.MaxDistance),...
                     sprintf('ConfinementRatio: %.3g',stats.ConfinementRatio),...
                     sprintf('MeanSpeed: %.3g um/s',stats.MeanSpeed),...
                     sprintf('NetSpeed: %.3g um/s',stats.NetSpeed),...
                     sprintf('LinearityRatio: %.3g',stats.LinearityRatio),...
                     sprintf('Dmle: %.5g um^2/s',stats.Dmle)};                
        end
    end %Protected static methods
    
    methods (Access=protected)
        function [trackIds,segIds] = getTrackSegmentInputIds(obj, trackIds, segIds)
            %This is a helper function to read in a list of trackIDs 
            % and segmentIDs and check that they are valid and return them
            % for actual processing.
            if nargin==1 || isempty(trackIds)
                trackIds = 1:obj.nTracks;
            end
            if nargin<3 || isempty(segIds)
                %Make an entry for each segment of each track
                segIds = cellmatfun(@(i) 1:obj.nSegments(i), trackIds(:)');
                trackIds = cellmatfun(@(i) repmat(i,1,obj.nSegments(i)), trackIds(:)');
            end
            if length(trackIds)~=length(segIds)
                error('TrackSegmentAnalysis:segmentStatsTable','Mismatched lengths for track / segment indexes');
            end
        end
    end %Protected methods
end %classdef
    
    
