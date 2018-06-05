classdef RegisteredPairAnalysis < PairAnalysis
    % This extends the PairAnalysis class to specifically work with registered 2-color data.  At the moment
    % we only support 2-colors, but theroretically could support more.
    %
    % A ReisgesterdPairAnalysis object can be created with 2 tracked RPT files together with a channel registration
    % object.  A central concept is that both RPT's come from the same dataset which is described and accessed with
    % a SPData object.
    %
    % In particular this class provides special visualization tools that are specific to 2-color registered data.
    properties (Constant=true, Hidden=true)
        %Abstract properties inherited from Pickle
        saveFileExt = '.regpairs'
        SaveableDataFormats = {'*.regpairs', 'RegisteredPairAnalysis (.regpairs)'};
        LoadableDataFormats = {'*.regpairs;*.spdata','All Loadable Formats (*.regpairs, *.spdata)';...
                                '*.regpairs', 'RegisteredPairAnalysis (.regpairs)';...
                                '*.spdata', 'SPData (.spdata)'};       
        NChannels=2; %Only 2 channels handled for now until we have a mircoscope that can do 4-channels...
        DefaultChannelWavelengths = [585 655];
    end

    properties (Dependent=true)      
        nChannelTracks; % number of tracks is a 2x1 cell array giving tracks for each color
    end %abstract dependent public properties

    properties
        %These propoerties help to organize the mapping of tracks to colors (wavelengths) according to their
        %respective RPT's
        channelWavelengths = [585 655]; % 1xNChannels vector giving the wavelengths of each of the two registered channels.
        trackChannels; %A vector mapping trackId to channels;
        trackChannelIdxs; %A 1xNChannels cell array where the two cells are vectors giving the indexes of all tracks 
                             %from that wavelength (coorresponding to the specific RPT)
        sumImage; % A size:[NChannels,1] cell array with the background sum image for both channels
        frameSize; %A size:[NChannels,1] cell array with the frame size in pixels for each channel
        pixelSize; % The [pixelSize_x, pixelSize_y]
        chMap; % The optimal channel mapping structure with the following fields:
               %  mapFunction - A function mapping Nx2 arrays of point locations to Nx2 arrays of point locations
               %                where the mapping translates ch2->ch1 points.  This is what is used to perform all
               %                mappings
               %  algorithsm - The name of the algorithm used
               %  params - A structure of parameters used for the given algorithm
               %  trainData - A cell array of index vectors describing which refPts were used for training
               %  testData - A cell array of index vectors describing which refPts were used for training
               %  RMSE - The root-mean-squared-error for the mapping as measured with the test data
               %  maxError - The maximum displacement error observed over the test data set.     
        channelPhysicalBounds; % cell array [NChannels,1], where for each channel we record the smallest axis-aligned rectangular bounding that encompasses the registered channels physical ROI coordinates.
        Paths=struct(... %All filenames are relative paths from the workingDir
         'saveFile',[],... %filename(w/extension) of the .tsa file where this will save to
         'data',[],... %relative path to the .spdata file we are associated with
         'rpt_channel1',[],... %relative path to the .rpt file for channel 1
         'rpt_channel2',[],... %relative path to the .rpt file for channel 2
         'regAnalysis',[]... %relative path to the channel .reganalysis file
         );

        minPairLength = 10; %Minimum lenght for a pair to be considered valid
        minPairDistance = 0.1; %Minimum approach distance necessary for a pair to be considered valid
    end
    
    %These properties are for the plotting and colormaps
    properties (Hidden=true)
        FigureSize = [1200,900];
        ChannelTrackColors = [0.5,1.0,0.0; 0.5,0.0,1]; % Green=channel1 and Magenta=channel2
        ImageColorMaps; % a 1xNChannels cell array of track colormaps
        TrackColorMethod = 'Wavelength';
        TrackColorMap; % set by obj.setTrackColors()
        ColorBarTitle; % set by obj.setTrackColors()
        ImageColorRange = 256; %The size of the image color maps.
        TrackColorRange = 256; %The size of the track color map (for appropriate TrackColorMethods)
        TrackValueRange = []; % A 2x1 vector giving min and max values for the tracks for internal display of the colorbar ticks
        ColorBarOn=true; % Show the colorbar?
        ColorMap; % Final colormap which is the concatenation of the image colormaps and thje track color map
        AxesFontSize = 12;
        TrackLineWidth = 1;
        TrackVertexMarker = 'none';
        TitleOn = true;
        TitleFontSize = 12;
        VideoWriterProfile = 'MPEG-4';
        VideoWriterFrameRate = 25;
        video;
        MovieSize=[800, 800];
        RenderNormIms;
        RenderImCoords;
        RenderTracks;
        RenderTrackCs;
        figureH;
    end
    

    
    properties (Transient=true)
        rpt_channel1;
        rpt_channel2;
        regAnalysis;
    end
    
    properties (Dependent=true, Hidden=true)
        nImageColors; % Lenght of image portion of color map
    end
    methods
        function obj=RegisteredPairAnalysis(varargin)
            % RegisteredPairAnalysis
            %    Arguments are the same as for a the call to obj.load()
            %
            %   A Registered pair analysis is a class to represent the combined data from a pairing of
            % 2xChannel split image registered tracking data.  (AKA 2-color tracking).  This class
            % provides a gui to aid in interaction.  [ Call: obj.gui() ]
            %
            %  An object of this type is created by passing the RPT tracking object used to track ther
            % ROIs for both channels.


            obj@PairAnalysis(); % call superclass
            if nargin>0
                obj.load(varargin{:});
            end
        end
        
        function load(obj, varargin)
            % [Case1]: Call with saved .regpairs file
            %  [in]: regpairs_file -> filename or object of RegisteredPairAnalysis
            %
            % [Case2]: Call with rpt files and channel registration
            %  [in] rpt1 - The tracked RPT for the first color.  This is used as the base of the coordinate system.
            %             [This argument may be passed as an RPT object or a complete path to a .rpt file]
            %  [in] rpt2 - The tracked RPT for the second color.  This is the registered channel that will be mapped
            %               to the coordinates of the primary channel.
            %             [This argument may be passed as an RPT object or a complete path to a .rpt file]
            %  [in] regAnalysis - The RegistrationAnalysis registration object (or complete path to saved .reganalysis file).
            %  [in] wavelengths - 1x2 vector giving the wavelength (in nm) of the respective color channels.
            % [Casee]: Call with spdata file and channel registration
            %   [in] spdata - An SPData filename or object which has the two channel ROIs
            %   [in] regAnalysis - The RegistrationAnalysis registration object (or complete path to saved .reganalysis file).
            %   [in] wavelengths - 1x2 vector giving the wavelength (in nm) of the respective color channels
            %   [in] channelROI - [optional] (default: {1,2}) - A 1x2 cell array giving the indexes or names of 
            %                                                   the ROI's for each channel in the SPData file.
            
            
            if nargin==2
                obj.loadRegPairs(varargin{1}); %Case1 Load from saved .regpairs file
            elseif ischar(varargin{1})
                [~,~,ext] = fileparts(varargin{1});
                switch ext
                    case SPData.saveFileExt
                        obj.loadSPData(varargin{:});
                    case RPT.saveFileExt
                        obj.loadRegisteredRPT(varargin{:});
                    otherwise
                        error('RegisteredPairAnalysis:load','Unknown load format');
                end
            elseif isa(varargin{1},'SPData')
                obj.loadSPData(varargin{:});
            elseif isa(varargin{1},'RPT')
                obj.loadRegisteredRPT(varargin{:});
            end
            obj.initialized = true;
        end
        
        function spdata = getData(obj)
           % Return the SPData for this set of registered pairs.  We retrive it from the RPT object
%            if ~isa(obj.rpt_channel1,'RPT')
%                error('RegisteredPairAnalysis:RPTs not set.  Cannot retrieve SPData.');
%            end
            if isa(obj.rpt_channel1,'RPT')
                spdata = obj.rpt_channel1.data;
            else
                spdata=[];
            end
        end
               
        function identifyPairs(obj, min_length, min_distance)
            % min_length: integer>=2.  Minimum number of simultaneous localizations for a pair to be
            %                          considered
            % min_distance: (um).  Minimum distance of approach for pair to be considered.
            obj.clearPairs(); % Clear all old pairs
            nT = obj.nChannelTracks;
            startFrame = cellfun(@(T) T(1,end), obj.tracks(obj.trackChannelIdxs{2}));
            obj.pairIds = zeros(nT{1}*nT{2},2);
            count = 0;
            for n = 1:nT{1}
                idx_ch1 = obj.trackChannelIdxs{1}(n);
                t1 = obj.tracks{idx_ch1};
                if size(t1,1)<min_length
                    continue; % Too short
                end
                for m = find(t1(end,end) - startFrame+1 >= min_length);
                    idx_ch2 = obj.trackChannelIdxs{2}(m);                   
                    t2 = obj.tracks{idx_ch2};
                    [~, loc_idx1, loc_idx2] = intersect(t1(:,end),t2(:,end));
                    if numel(loc_idx1) < min_length
                        continue; % Not enough simultaneous localizations
                    end
                    min_approach_dist = sqrt(min(sum((t1(loc_idx1,2:3)-t2(loc_idx2,2:3)).^2,2)));
                    if  min_distance <  min_approach_dist;
                        continue; % Do not come close enough
                    end
                    count = count+1;
                    obj.pairIds(count,1) = idx_ch1;
                    obj.pairIds(count,2) = idx_ch2;
                end                
            end
            obj.pair_min_length = min_length;
            obj.pair_min_distance = min_distance;
            obj.pairIds   = obj.pairIds(1:count,:);
            obj.pairs     = obj.makePair(1:count);
            obj.pairStats = cellfun(@obj.makePairStats, obj.pairs); % Pair stats will be a struct array
            obj.pairDists = cellmap(@obj.makePairDists, obj.pairs); % Pair dists will be a cell array of distance maticies
            obj.pairsIdentified = true;
        end



        function pair_mats = makePair(obj,pairidxs)
            Npairs = numel(pairidxs);
            pair_mats = cell(Npairs,1);
            for n=1:Npairs 
                T1 = obj.tracks{obj.pairIds(pairidxs(n),1)};
                T2 = obj.tracks{obj.pairIds(pairidxs(n),2)};
                [~,i1,i2] = intersect(T1(:,end), T2(:,end));
                N = numel(i1);
                wv1 = ones(N,1)*obj.channelWavelengths(1);
                wv2 = ones(N,1)*obj.channelWavelengths(2);
                zs = zeros(N,1);
                pair_mats{n} = [T1(i1,1),T1(i1,2:3),wv1,T1(i1,4),T2(i2,2:3),wv2,T2(i2,4),...
                                T1(i1,7:8),zs,T1(i1,9),T2(i2,7:8),zs,T2(i2,9), T1(i1,end)];
            end
        end
        
        function pair_mats = makePairTrajectories(obj,pairidxs)
            % Return the smaller pair format for the interaction change point analysis
            % [t x1 y1 x2 y2 SEx1 SEy1 SEx2 SEy2]
            Npairs = numel(pairidxs);
            pair_mats = cell(Npairs,1);
            for n=1:Npairs 
                T1 = obj.tracks{obj.pairIds(pairidxs(n),1)};
                T2 = obj.tracks{obj.pairIds(pairidxs(n),2)};
                [~,i1,i2] = intersect(T1(:,end), T2(:,end));
                pair_mats{n} = [T1(i1,1), T1(i1,2:3),T2(i2,2:3),T1(i1,7:8),T2(i2,7:8)];
            end
        end
        
        function stats = makePairStats(obj, pair)
            stats.length = size(pair,1);
            stats.duration_frames = pair(end,end)-pair(1,end)+1;
            stats.duration_secs = pair(end,1)-pair(1,1);
            stats.min_distance = sqrt(min(sum( (pair(:,2:3) - pair(:,6:7)).^2,2)));
            stats.bbox = obj.pairPhysicalBounds(pair);
            dest = DEstimator(pair(:,2:3), pair(:,1), pair(:,10:11), obj.frameT);
            stats.Dmle1 = dest.MLE();
            dest.initializeTrack(pair(:,6:7), pair(:,1), pair(:,14:15), obj.frameT);
            stats.Dmle2 = dest.MLE();
        end

        function icp = makeInteractionChangePoint(obj, pairid)
            D_A = 0.020;
            D_B = 0.020;
            D_AB = 0.010;
            exposureT = obj.rpt_channel1.data.frameT;
            rho = 0.050;
            icp = InteractionChangePoint(D_A, D_B, D_AB, exposureT, rho);
            pair_mat = obj.makePair(pairid);
            icp.initializePair(pair_mat);
        end

        function initializeColorMaps(obj)
            % Call this to reset the colormaps
            delta = linspace(0,1,obj.ImageColorRange)';
            zs = zeros(obj.ImageColorRange,1);
            logistic_rate = 10;
            logistic_center = 0;
            
            obj.ImageColorMaps{1} = [zs, delta, zs]; %green            
%             delta = 0.2 + 0.8*sqrt(delta);
            delta = 1./(1+exp(-logistic_rate*(delta-logistic_center))); %Logistic function to make S response curve
            delta = delta./max(delta);
            obj.ImageColorMaps{2} = [delta, zs, delta]; %magenta
        end
        
        function vw = initializeVideoSequence(obj, filename)
            
            obj.video = VideoWriter(filename, obj.VideoWriterProfile);
            obj.video.FrameRate = obj.VideoWriterFrameRate;
            obj.video.open();
%             obj.figureH=figure('Units','Pixels','Position',[0 0, obj.MovieSize]);
            
            vw = obj.video;
        end
        
        function renderActiveTracks(obj, Tidxs, method, pBounds, tBounds)
%             boarder= obj.FigureBoarderSize;
%             if obj.ColorbarOn %adjust right figure boarder if color bar is used
%                 boarder(3) = boarder(3) + 70;
%             end
            
            if nargin<3
                method = 'Wavelength';
            end
            if nargin<2 || isempty(Tidxs);
                Tidxs = 1:obj.nTracks;
                pBounds = obj.physicalBounds;
                tBounds = obj.timeBounds;
                fBounds = obj.frameBounds;
            else
                if nargin<4 || isempty(pBounds)
                    pBounds = RPT.trackPhysicalBounds(obj.tracks(Tidxs));
                end
                if nargin<5 || isempty(tBounds)
                    tBounds = RPT.trackTimeBounds(obj.tracks(Tidxs));
                    fBounds = RPT.trackFrameBounds(obj.tracks(Tidxs));
                else
                    fBounds = round(tBounds./obj.frameT);
                end
            end
            obj.figureH=figure('Units','pixels','Position',[5,5,obj.MovieSize]);
            obj.figureH.Color=[0 0 0];
            whitebg(obj.figureH);
            axH=axes();
            obj.RenderTracks = obj.tracks(Tidxs);
            obj.RenderTrackCs = obj.setTrackColors(Tidxs,method);
            obj.initializeAxes(axH, pBounds, tBounds);
            axH.Units='pixels';
            obj.figureH.Color = [0 0 0];
            GUIBuilder.positionAxes(axH, [0, 0, obj.MovieSize], [0 10 130 80]); 
            [nim1, nim2] = obj.normalizeFrameColors(obj.rpt_channel1.getFrames(), obj.rpt_channel2.getFrames());
            obj.RenderNormIms={nim1,nim2};
            obj.RenderImCoords = obj.makeImageSurfaceCoords(tBounds(1));
            lag =0;
%             imHs = obj.plotBackgroundFrames(fBounds(1));
%             trackHs = obj.plotTracksFrame(fBounds(1));
%             locHs = obj.plotLocsFrame(fBounds(1));
            if obj.TitleOn
                title(sprintf('Frame:%i Time %.3f (s)',fBounds(1),(fBounds(1)-1)*obj.frameT),'interpreter','latex','fontsize',obj.TitleFontSize);
            end
            frameIdxs = fBounds(1):fBounds(2);
            deltaV = [30-0, 10-90]./50;
            v = [0,90];
            for i=frameIdxs
                imH=obj.plotBackgroundFrames(i);
                if lag>0
                    trackHs=obj.plotTracksFrame(i,i-lag);
                else
                    trackHs=obj.plotTracksFrame(i);
                end
                locHs=obj.plotLocsFrame(i);
                if obj.TitleOn
                    title(sprintf('Frame:%i Time %.3f (s)',i,(i-1)*obj.frameT),'interpreter','latex','fontsize',obj.TitleFontSize);
                end
                drawnow();
                F=getframe(obj.figureH);
                if i>=100 && i<150
                    v=v+deltaV;
                    view(v);
                elseif i>=200 && i<250
                    v=v-deltaV;
                    view(v);
                end
                
                obj.video.writeVideo(F);
%                 pause(0.1);
                delete(imH);
                delete(trackHs);
                delete(locHs);
            end
        end
        
        function finalizeVideoSequence(obj)
            obj.video.close();
            close(obj.figureH);
        end
    
        function plotPairs(obj, axH, pairIdxs, colorMethod)
            % physicalBounds and tBounds are set to include all tracks for all pairIds
            if nargin<2 || isempty(axH) || ~ishandle(axH)
                F=figure('Position',[10,10,obj.FigureSize]);
                whitebg(F);
                F.Color = [0,0,0];
                axH = axes();
            end
            if nargin<3 || isempty(pairIdxs)
                pairIdxs = 1:obj.nParis;
            end
            if nargin<4
                colorMethod = 'Wavelength';
            end
            Tidx = obj.getPairTrackIdx(pairIdxs);
            obj.plotTracks(axH, Tidx, colorMethod);
        end
        
        function viewPairMovie(obj, pairIdxs, colorMethod)
            if nargin<2 || isempty(pairIdxs)
                pairIdxs = 1:obj.nParis;
            end
            if nargin<3
                colorMethod = 'Wavelength';
            end
            Tidx = obj.getPairTrackIdx(pairIdxs);
            obj.viewTracksMovie(Tidx, colorMethod);
        end
   
        function viewTracks(obj, varargin)
            obj.plotTracks([], varargin{:});
        end

        function plotAllPairsWavelength(obj, axH)
            Tidxs = unique(obj.pairIds(:));
            obj.plotTracks(axH,Tidxs,'Wavelength');
        end
        
        function plotTracks(obj, axH, Tidxs, colorMethod, pBounds, tBounds)
            if obj.nTracks<1
                return;
            end
            if nargin<2 || isempty(axH) || ~ishandle(axH)
                figure('Position',[10,10,obj.FigureSize],'Color',[0 0 0]);
                axH = axes();
            end
            if nargin<4
                colorMethod = 'Wavelength';
            end
            if nargin<3 || isempty(Tidxs);
                Tidxs = 1:obj.nTracks;
                pBounds = obj.physicalBounds;
                tBounds = obj.timeBounds;
            else
                if nargin<5 || isempty(pBounds)
                    pBounds = RPT.trackPhysicalBounds(obj.tracks(Tidxs));
                end
                if nargin<6 || isempty(tBounds)
                    tBounds = RPT.trackTimeBounds(obj.tracks(Tidxs));
                end
            end
            Cs = obj.setTrackColors(Tidxs, colorMethod);            
            [nim1, nim2] = obj.normalizeFrameColors(obj.sumImage{1}, obj.sumImage{2});
            imCoords = obj.makeImageSurfaceCoords(tBounds(1));
            obj.initializeAxes(axH, pBounds, tBounds);
%             surface(imCoords{1,:},zeros(size(nim1)),'FaceColor','texturemap',...
%                 'EdgeColor','none','CDataMapping','direct','FaceAlpha',1);
            alphaData = nim2;
            alphaData(:) = alphaData(:)-min(alphaData(:));
            alphaData(:) = 64*(0.8*alphaData(:)/max(alphaData(:))+0.2);
            surface(imCoords{1,:},nim1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct',...
                    'FaceAlpha',1);
            surface(imCoords{2,:},nim2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct',...
                    'FaceAlpha','texturemap','AlphaData',alphaData,'AlphaDataMapping','direct');       
            for i = 1:numel(Tidxs)
                RPT.drawTrackSurface(obj.tracks{Tidxs(i)}, Cs{i}, 'LineWidth', obj.TrackLineWidth);                
            end
            view([-10,30]);
        end
    end %public methods

    methods
        function nT = get.nChannelTracks(obj)
            nT = cellmap(@numel, obj.trackChannelIdxs);
        end
        function nColors = get.nImageColors(obj)
            % Length of image portion of color map
            nColors = 1+2*obj.ImageColorRange;
        end
    end % Dependent properties

    methods (Access=protected)        
        function loadRegPairs(obj, regpairs_file)
            [wDir,saveFileBase,~] = fileparts(regpairs_file); 
            s = load(regpairs_file,'-mat');
            obj.copyobj(s.obj);
            obj.workingDir = collapsepath(wDir);
            obj.Paths.saveFile = [saveFileBase obj.saveFileExt];
            obj.rpt_channel1 = RPT(obj.getFilePath('rpt_channel1'));
            obj.rpt_channel2 = RPT(obj.getFilePath('rpt_channel2'));
            obj.regAnalysis = RegistrationAnalysis(obj.getFilePath('regAnalysis'));
            obj.initialized = true;
            obj.dirty = false;
        end
        
        function loadSPData(obj, spdata, regAnalysis, wavelengths, channelROI)
            % channelROI [optional] - a cell array of ROI names or ROI
            if ischar(spdata)
                spdata = SPData(spdata);
            elseif ~isa(spdata, 'SPData')
                error('RegisteredPairAnalysis:load', 'Invalid SPData object');
            end
            if nargin<5
                channelROI = {1,2};
            elseif ~iscell(channelROI)
                channelROI = num2cell(channelROI);
            end
            rpt_path1 = spdata.getROIFiles(channelROI{1}, 'RPT', true);
            rpt_path2 = spdata.getROIFiles(channelROI{2}, 'RPT', true);
            if isempty(rpt_path1) || isempty(rpt_path2)
                error('RegisteredPairAnalysis:load', 'SPData object is not fully tracked.  Cannot find .rpt files');
            end
            obj.loadRegisteredRPT(rpt_path1, rpt_path2, regAnalysis, wavelengths);                
        end
                
        function loadRegisteredRPT(obj, rpt1, rpt2, regAnalysis, wavelengths)
            % Create a new RegisteredPairAnalysis object with 2 rpt files and correpsonding channel registration
            %  [in] rpt1 - The tracked RPT for the first color.  This is used as the base of the coordinate system.
            %             [This argument may be passed as an RPT object or a complete path to a .rpt file]
            %  [in] rpt2 - The tracked RPT for the second color.  This is the registered channel that will be mapped
            %               to the coordinates of the primary channel.
            %             [This argument may be passed as an RPT object or a complete path to a .rpt file]
            %  [in] channelReg - The channel registration object (or complete path to saved file).
            %  [in] wavelengths - 1x2 vector giving the wavelength (in nm) of the respective color channels.           
                        
            %load rpt1
            if ischar(rpt1)
                obj.rpt_channel1 = RPT(rpt1);
            elseif isa(rpt1,'RPT')
                obj.rpt_channel1 = rpt1;
            else
                error('RegisteredPairAnalysis:load','Unkown RPT format for rpt1');
            end
            %load rpt2
            if ischar(rpt2)
                obj.rpt_channel2 = RPT(rpt2);
            elseif isa(rpt1,'RPT')
                obj.rpt_channel2 = rpt2;
            else
                error('RegisteredPairAnalysis:load','Unkown RPT format for rpt2');
            end
            %load SPData object
            obj.Paths.data = obj.rpt_channel1.data.workingDir;
            %Load channel registration
            if ischar(regAnalysis)
                obj.regAnalysis = RegistrationAnalysis(regAnalysis);
            elseif isa(regAnalysis,'RegistrationAnalysis')
                obj.regAnalysis = regAnalysis;
            else
                error('RegisteredPairAnalysis:load','Channel registration requires full path to registration .mat file in same directory tree');
            end
            %Check RPT files are OK
            if ~strcmp( obj.rpt_channel1.getFilePath('data'),  obj.rpt_channel2.getFilePath('data'))
                error('RegisteredPairAnalysis:load','RPT1 and RPT12 are not from the same SPData');
            else
                obj.rpt_channel1.data = obj.rpt_channel2.data; %make sure we only have 1 data open.
            end            
            if ~strcmp(obj.rpt_channel1.phase,'Tracked')
                error('RegisteredPairAnalysis:load','RPT1 is not tracked');
            end
            if ~strcmp(obj.rpt_channel2.phase,'Tracked')
                error('RegisteredPairAnalysis:load','RPT2 is not tracked');
            end
            %Setup and record paths
            obj.initialized = false;
            data = obj.rpt_channel1.data;
            obj.workingDir = fullfile(data.workingDir,'RegisteredPairAnalysis');
            obj.Paths.data = relativepath(obj.workingDir, data.saveFilePath);
            obj.Paths.rpt_channel1 = relativepath(obj.workingDir, obj.rpt_channel1.saveFilePath);
            obj.Paths.rpt_channel2 = relativepath(obj.workingDir, obj.rpt_channel2.saveFilePath);
            obj.Paths.regAnalysis = relativepath(obj.workingDir, obj.regAnalysis.saveFilePath);
            [~,saveFileBase,~] = fileparts(obj.Paths.data);
            roi1_name = obj.rpt_channel1.ROIname;
            roi2_name = obj.rpt_channel2.ROIname;
            obj.Paths.saveFile = sprintf('%s_%s_%s%s',saveFileBase, roi1_name, roi2_name, obj.saveFileExt);

            %Check wavelengths
            assert(numel(wavelengths)==2 && all(wavelengths>0) && wavelengths(1)~=wavelengths(2));
            obj.channelWavelengths = wavelengths;
            %Get tracks and register them
            nT1 = obj.rpt_channel1.nTracks;
            nT2 = obj.rpt_channel2.nTracks;
            if nT1==0
               error('RegisteredPairAnalysis:load','RPT1 has no tracks');
            end
            if nT2==0
               error('RegisteredPairAnalysis:load','RPT2 has no tracks');
            end
            T1 = obj.rpt_channel1.getTracks();
            T2 = obj.rpt_channel2.getTracks();
            obj.pixelSize(1:2) = data.pixelSize(:)'; % enusre it is a row vector
            obj.sumImage{1} = obj.rpt_channel1.sumImage;            
            obj.sumImage{2} = obj.rpt_channel2.sumImage;            
            [~,obj.chMap] = obj.regAnalysis.getOptimalMapMicrons();
            reg_T2 = cellmap(@obj.registerTrack, T2); % registered channel 2 tracks
            obj.tracks = [T1, reg_T2]; %all tracks together in a single list
            obj.trackChannels = [ones(nT1,1); 2*ones(nT2,1)];
            obj.trackChannelIdxs = {1:nT1, nT1+1:nT1+nT2};

            %Check ROI Bounds
            if obj.rpt_channel1.ROI(2)>obj.rpt_channel2.ROI(1)
                error('RegisteredPairAnalysis:load','Pixel-based ROIs intersect for RPT1 and RPT2 on x-dimension');
            end
            %Record physical bboxes for ROI's
            obj.frameSize{1} = obj.rpt_channel1.frameSize;
            obj.frameSize{2} = obj.rpt_channel2.frameSize;
            obj.channelPhysicalBounds{1} = obj.rpt_channel1.ROIPhysical;
            obj.channelPhysicalBounds{2} = obj.rpt_channel2.ROIPhysical;
            ch2ROI = obj.chMap.mapFunctionMicrons(obj.channelPhysicalBounds{2}([1,3;1,4;2,3;2,4]));
            reg_channelPhysicalBounds = [min(ch2ROI(:,1)),max(ch2ROI(:,1)),min(ch2ROI(:,2)),max(ch2ROI(:,2))];
            obj.physicalBounds = [min(obj.channelPhysicalBounds{1}(1),reg_channelPhysicalBounds(1)),... %xmin
                                  max(obj.channelPhysicalBounds{1}(2),reg_channelPhysicalBounds(2)),... %xmax
                                  min(obj.channelPhysicalBounds{1}(3),reg_channelPhysicalBounds(3)),... %ymin
                                  max(obj.channelPhysicalBounds{1}(4),reg_channelPhysicalBounds(4))];   %ymax
            obj.frameT = data.frameT;
            obj.frameBounds = RPT.trackFrameBounds(obj.tracks);
            obj.timeBounds = RPT.trackTimeBounds(obj.tracks);
            obj.identifyPairs(obj.minPairLength, obj.minPairDistance);
            obj.initializeColorMaps();
            obj.dirty = true;
        end
        
        function reg_track = registerTrack(obj, track_ch2)
            %register a single track
            reg_track = track_ch2;
            reg_track(:,[2,3]) = obj.chMap.mapFunctionMicrons(reg_track(:,[2,3]));
        end
        

        function [nframes1,nframes2] = normalizeFrameColors(obj, frames1, frames2)
            stretch = @(F) 1./(1+exp(-4*(cosmicNorm(F)-0.)));
            stretch = @cosmicNorm;
            nframes1 = cosmicNorm(frames1)*(obj.ImageColorRange-1)+1;
            nframes2 = stretch(frames2)*(obj.ImageColorRange-1)+1+obj.ImageColorRange;
%             obj.imageSurfaceCoords{1}={
        end
        
        function imCoords = makeImageSurfaceCoords(obj, tmin)
            %Compute the coordinate grid for the pixels of frames from channel 1 and channel 2 so as to be
            % inputs for the surface() function.
            xs = linspace(obj.channelPhysicalBounds{1}(1), obj.channelPhysicalBounds{1}(2), obj.frameSize{1}(2)+1);
            ys = linspace(obj.channelPhysicalBounds{1}(3), obj.channelPhysicalBounds{1}(4), obj.frameSize{1}(1)+1);
            [AX, AY] = meshgrid(xs,ys);
            AZ = repmat(tmin, obj.frameSize{1}+1);
            xs = linspace(obj.channelPhysicalBounds{2}(1), obj.channelPhysicalBounds{2}(2), obj.frameSize{2}(2)+1);
            ys = linspace(obj.channelPhysicalBounds{2}(3), obj.channelPhysicalBounds{2}(4), obj.frameSize{2}(1)+1);
            [BX, BY] = meshgrid(xs,ys);
            pts = [BX(:), BY(:)];
            rpts = obj.chMap.mapFunctionMicrons(pts);
            BX(:) = rpts(:,1);
            BY(:) = rpts(:,2);
            BZ = repmat(tmin, obj.frameSize{2}+1);
            imCoords = {AX, AY, AZ; BX, BY, BZ};
        end



        function Cs = setTrackColors(obj, Tidxs, method)
            if nargin==2
                if ~isempty(obj.TrackColorMethod)
                    method=obj.TrackColorMethod;
                else
                    method='Wavelength';
                end
            end
            colorstart = 2*obj.ImageColorRange;
            nT = numel(Tidxs);
            switch method
                case 'Wavelength'
                    Cs = cellmap(@(i) ones(size(obj.tracks{i},1),1)*obj.trackChannels(i)+colorstart, Tidxs);
                    obj.TrackValueRange=[1, 2];
                    obj.TrackColorMap = obj.ChannelTrackColors;
                    obj.ColorBarTitle='';
                case 'Sequence'
                    Cs = cellmap(@(i) ones(size(obj.tracks{Tidxs(i)},1),1)*i+colorstart, 1:nT);
                    obj.TrackValueRange=[1, nT];
                    obj.TrackColorMap = [0,1,0; 1,0,1; prism(nT-2)];
                    obj.ColorBarTitle='Sequence';
                case 'Speed'
                    Cs = cellmap(@(T) TrackSegmentAnalysis.estimateSpeedWindow(T,6), obj.tracks(Tidxs));
                    mx = max(cellfun(@max,Cs));
                    obj.TrackValueRange=[0, mx];
                    Cs = cellmap(@(C) C./(mx/obj.TrackColorRange)+colorstart+1, Cs);
                    obj.TrackColorMap = jet(obj.TrackColorRange);
                    obj.ColorBarTitle='Speed (um/s)';
                case 'Temporal'
                    Cs = cellmap(@(T) T(:,1), obj.tracks(Tidxs));
                    mx = max(cellfun(@max,Cs));
                    obj.TrackValueRange = [0, mx];
                    Cs = cellmap(@(C) C./(mx/obj.TrackColorRange)+colorstart+1, Cs);
                    obj.TrackColorMap = jet(obj.TrackColorRange);
                    obj.ColorBarTitle='Time (s)';
                otherwise
                    error('RegisteredPairAnalysis:setTrackColors',  'Unknown track color method: %s', method);
            end
            obj.ColorMap = [obj.ImageColorMaps{1}; obj.ImageColorMaps{2}; obj.TrackColorMap];
            obj.TrackColorMethod = method;
        end
        
        function initializeAxes(obj, axH, bbox, tbounds)
            axH.Units='Pixels';
            axH.SortMethod='ChildOrder';
            axH.YDir='reverse';
            axH.TickDir='out';
            axH.XMinorTick='on';
            axH.YMinorTick='on';
            axH.ZMinorTick='on';
            grid(axH,'on');
            grid(axH,'minor');
            axH.Box='on';
            axH.BoxStyle='full';
            axH.Projection='Orthographic';

            axis([bbox, tbounds]);
            asp(1) = (bbox(2)-bbox(1))/(bbox(4)-bbox(3));
            asp(2) = 1;
            asp(3) = mean(asp(1:2)); % golden ratio
            pbaspect(asp);
            view(0,90);

            zlabel('t (s)','interpreter','latex','fontsize',obj.AxesFontSize);
            xlabel('x (um)','interpreter','latex','fontsize',obj.AxesFontSize);
            ylabel('y (um)','interpreter','latex','fontsize',obj.AxesFontSize);
            
            colormap(axH, obj.ColorMap);
            if obj.ColorBarOn
                colorbarH = colorbar();
                colorbarH.YLim = [0 size(obj.TrackColorMap,1)]+obj.nImageColors;
                colorbarH.Label.String = obj.ColorBarTitle;
                colorbarH.Label.Interpreter = 'latex';
                colorbarH.Label.FontSize = obj.AxesFontSize;
                switch obj.TrackColorMethod
                    case 'Wavelength'
                        colorbarH.Ticks = obj.TrackValueRange-0.5+obj.nImageColors;
                        colorbarH.TickLabels = {sprintf('Ch1:%inm',obj.channelWavelengths(1)),...
                                                sprintf('Ch2:%inm',obj.channelWavelengths(2))};
                    case 'Sequence'
                        nT =  obj.TrackValueRange(2);
                        ticks = round(linspace(1, nT, min(nT,20)));
                        colorbarH.Ticks = ticks + obj.nImageColors-0.5;
                        colorbarH.TickLabels = cellmap(@num2str,ticks);
                    otherwise
                        ticks = linspace(obj.TrackValueRange(1), obj.TrackValueRange(2), 12);
                        colorbarH.Ticks = (ticks-obj.TrackValueRange(1)).* (obj.TrackColorRange./obj.TrackValueRange(2)) + obj.nImageColors;
                        colorbarH.TickLabels = cellmap(@(n) num2str(n,'%.3f'),ticks);
                end
            end
            hold(axH,'on');

        end
        
        function Hs = plotBackgroundFrames(obj,frameIdx)
            Hs = cell(2,1);
            alphaData = obj.RenderNormIms{2}(:,:,frameIdx);
            alphaData(:) = alphaData(:)-min(alphaData(:));
            alphaData(:) = (alphaData(:)/max(alphaData(:))).^.7;
            
            
            Hs{1} = surface(obj.RenderImCoords{1,:},obj.RenderNormIms{1}(:,:,frameIdx),'FaceColor','texturemap',...
                    'EdgeColor','none','CDataMapping','direct','FaceAlpha',1);%'texturemap','AlphaData',obj.RenderNormIms{1}(:,:,frameIdx),'AlphaDataMapping','scaled');
            Hs{2} = surface(obj.RenderImCoords{2,:},obj.RenderNormIms{2}(:,:,frameIdx),'FaceColor','texturemap',...
                    'EdgeColor','none','CDataMapping','direct','FaceAlpha','texturemap','AlphaData',alphaData,'AlphaDataMapping','scaled');
            Hs=[Hs{:}];
        end
        
        function trackHs = plotTracksFrame(obj, maxframe, minframe)
            % IN:
            %  maxframe - [optional] the maximum frame to show tracks for [default: last frame]
            %  minframe - [optional] the minimum frame to show tracks for [default: first frame]
            if nargin<3
                minframe = obj.frameBounds(1);
            end
            if nargin<2
                maxframe = obj.frameBounds(2);
            end
            nT = numel(obj.RenderTracks);
            trackHs = cell(1,nT);
            for i = 1:nT
                T = obj.RenderTracks{i};
                Cs = obj.RenderTrackCs{i};
                minIdx = find(T(:,end)>=minframe, 1, 'first');
                maxIdx = find(T(:,end)<=maxframe, 1, 'last');
                if isempty(minIdx) || isempty(maxIdx); continue; end %Track not visible
                if minIdx==maxIdx; continue; end %Only a single localization.  Too short to draw
                trackHs{i} = RPT.drawTrackSurface(T(minIdx:maxIdx,:), Cs(minIdx:maxIdx),...
                                     'LineWidth',obj.TrackLineWidth,'Marker',obj.TrackVertexMarker);
            end
            trackHs = [trackHs{:}];
        end
        
        function locHs = plotLocsFrame(obj, frameIdx)
            %Make a half-sphere of given size with alpha shading based on distance
            %This is the object we will use to represent localizations
            S = BaseRPT.makeLocalizationSurface();
            nT = numel(obj.RenderTracks);
            locHs = cell(1,nT);
            for i = 1:nT
                T = obj.RenderTracks{i};
                Cs = obj.RenderTrackCs{i};
                Lidx = find(T(:,end)==frameIdx, 1, 'first');
                if isempty(Lidx); continue; end %Track does not have a localization this frame
                locHs{i} = RPT.drawLocalizationSurface(T(Lidx,[2,3]), 3*T(Lidx,[7,8]), T(1,1), Cs(Lidx), S);
            end
            locHs = [locHs{:}];
        end

        %% Abstract methods inherited from Pickle
        function val = getProtectedProperty(obj, name)
            %This is necessary for Pickle functionality to be able to access subclass protected variables
            val = obj.(name);
        end

        function modifyProtectedProperty(obj, name, newval)
            %This is necessary for Pickel functionality to be able to change subclass variables
            obj.(name)=newval;
        end
    end % protected methods

    methods (Static = true)
        function pairs_files=batchProcess(spdatapath, spdatafile_patterns, regAnalysis, wavelengths, overwriteFlag) 
            % Make and save a .regpairs for each of the .spdata's found using a specific RegistrationA
            % 
            % This will also auto track the files based on the default parameters.
            %
            % [IN]
            %  spdataPath - A path to a directory where .spdata files are to be registered as pairs
            %  spdatafile_patterns - A cell-array of file patterns that can use wildcard '*'
            %                      to search for (a pattern can also just be a filename with no wildacard)
            %                      Ex: {'2015-01-01*.spdata', '2015-01-03-condition2.spdata'}
            %  registrationAnalysis - The RegistrationAnalysis object or path to use for channel
            %                  registration
            %  defaultParams - An example RegisteredPairsAnalysis object or filename or preservedParams 
            %                  struct.
            %  overwriteFlag - [optional] Integer 0=Do not overwrite; [Default] 
            %                                     1=Warn with selection dialog before overwrite;
            %                                     2=Force overwrite (caution!);
            % [OUT]
            %  regpairs_files - Cell array of full-paths to all .regpairs files corresponding to the given 
            %             data file patterns. 
            %             (This list inculdes the default and non-overwitten files, as we assum you will
            %             want to process all these files similarly in the next phase of batch
            %             processing).
            
            if nargin<4 || numel(wavelengths)~=2
                wavelengths = RegisteredPairAnalysis.DefaultChannelWavelengths;
            end
            if nargin<5
                overwriteFlag = 0;
            end
            spdatafile_patterns = makecell(spdatafile_patterns);            
            spd_filenames = cellmap(@(p) Pickle.listExistingFileNames(spdatapath,p), spdatafile_patterns);
            spd_filenames = [spd_filenames{:}];
            if isa(regAnalysis,'RegistrationAnalysis')
                reg = regAnalysis;
            else
                reg = RegistrationAnalysis(regAnalysis);
            end
            
            Ndatafiles = numel(spd_filenames);
            Nprocessed = 0;
            Nerror = 0;
            Nexisting = 0;
            pairs_files = cell(1,Ndatafiles);
            H=waitbar(0,'Batch Processing ...');
            for n=1:Ndatafiles
                waitbar(n/(Ndatafiles+1),H,sprintf('Batch Processing Datafile [%i/%i] ...', n,Ndatafiles));                    
                fprintf('\n*** [%i/%i] Batch Processing Dataset:%s\n',n,Ndatafiles,spd_filenames{n}); 
                spd = SPData(spd_filenames{n});
                rpt1 = spd.getROIFiles(1, 'RPT');
                rpt2 = spd.getROIFiles(2, 'RPT');
                if isempty(rpt1) || isempty(rpt2)
                    fprintf('>>>(oops)<<< RPT files do not exist for "%s".  Please Batch track both ROIs first.\n',...
                            spd.saveFileBaseName);
                    Nerror = Nerror +1;
                end
                file_path = fullfile(spd.workingDir, 'RegisteredPairAnalysis',...
                                        [spd.saveFileBaseName, RegisteredPairAnalysis.saveFileExt]);
                if ~isempty(file_path) && ~overwriteFlag
                    pairs_files{n} = file_path; % Record file
                    Nexisting = Nexisting+1;
                    continue %Don't overwrite the existing files
                end
                try
                    pairs = RegisteredPairAnalysis(rpt1{1},rpt2{1}, reg, wavelengths);
                    pairs.save();
                    pairs_files{n} = pairs.saveFilePath;
                catch err
                    fprintf('>>>(oops)<<<  Batch processing error.\n');
                    disp(getReport(err))
                    pairs_files{n} = [];
                    Nerror = Nerror+1;
                    continue;
                end
                Nprocessed = Nprocessed +1;
            end            
            close(H);
            Ntotal = Nprocessed+Nexisting+Nerror;
            fprintf('\n*** Batch Process Complete. [Num Total: %i | Num Processed: %i | Num Existing: %i | Num Error: %i]\n',...
                        Ntotal,Nprocessed,Nexisting,Nerror);   
        end
    end %Public static methods
end

