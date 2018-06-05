% RPT.m
% Mark J. Olah (mjo@cs.unm.edu)
% 11/13/14
%
% Robust Particle Tracking.

classdef RPT < BaseRPT
    properties (Constant=true, Hidden=true)
        %Abstract properties inherited from Pickle
        saveFileExt = '.rpt'
        SaveableDataFormats = {'*.rpt', 'RPT (.rpt)'};
        LoadableDataFormats = {'*.spdata;*.rpt','All Loadable Sources (.spdata,.rpt)';...
                             '*.spdata','SPData file (.rpt)'; '*.rpt', 'RPT file (.rpt)'};

        %Abstract properties inherited from BaseRPT

        %These methods allow the generic methods of BaseRPT to do alot of common work
        DataClass = 'SPData';
        DataFileExt = SPData.saveFileExt;

        % Emitters are an format mainly for internal use.  
        % The information is the similar to the Localizations, but we leave more information in pixel
        % coordinates and retain mapping to boxCoords index in boxIdx column.  This allows association
        % between the emitter information and the box information
        NEmitterColumns = 13;
        EmitterColumnNames = {'x', 'y', 'I', 'bg', 'sigma', 'SE_x', 'SE_y', 'SE_I', 'SE_bg', 'SE_sigma','LLH','boxIdx','frameIdx'};
        EmitterColumnUnits = {'px', 'px', 'photons', 'photons/px', 'px', 'px', 'px', 'photons', 'photons/px', 'px','','index','index'};
        EmitterColumnDescriptions = {'x-Position estimate', 'y-Position estimate', 'Intensity estimate',...
                                   'Mean background intensity per pixel estimate','Apparent gaussian sigma estimate',...
                                   'Standard error of x-position estimate','Standard error of y-position estimate',...
                                   'Standard error of intensity estimate', 'Standard error of background intensity estimate',...
                                   'Standard error of apparent gaussian sigma estimate','Log Likelihood of fit',...
                                   'Index of box in BoxCoords structure','Frame index'};

        % Table/Matrix columns format for Localization results
        % This is the prefered output format for getting the localization information without the track
        % assoications
        NLocalizationColumns = 12;
        LocalizationColumnNames = {'t', 'x', 'y', 'I', 'bg', 'sigma', 'SE_x', 'SE_y', 'SE_I', 'SE_bg', 'SE_sigma','frame'};
        LocalizationColumnUnits = {'s', 'um', 'um', 'photons', 'photons/px', 'um', 'um', 'um', 'photons', 'photons/px', 'um','index'};
        LocalizationColumnDescriptions = {'Time', 'x-Position estimate', 'y-Position estimate', 'Intensity estimate',...
                                   'Mean background intensity per pixel estimate','Apparent gaussian sigma estimate',...
                                   'Standard error of x-position estimate','Standard error of y-position estimate',...
                                   'Standard error of intensity estimate', 'Standard error of background intensity estimate',...
                                   'Standard error of apparent gaussian sigma estimate',...
                                   'Frame index'};

        % Table/Matrix columns format for Track results
        % This defines the "RPT" track format although we can covert between
        % (1) cell-array of matricies. One per track
        % (2) single table with a track index
        % (3) structure array
        NTrackColumns = 12;
        TrackColumnNames = {'t', 'x', 'y', 'I', 'bg', 'sigma', 'SE_x', 'SE_y', 'SE_I', 'SE_bg', 'SE_sigma','frame'};
        TrackColumnUnits = {'s', 'um', 'um', 'photons', 'photons/px', 'um', 'um', 'um', 'photons', 'photons/px', 'um','index'};        
        TrackColumnDescriptions = {'Time', 'x-Position estimate', 'y-Position estimate', 'Intensity estimate',...
                                   'Mean background intensity per pixel estimate','Apparent gaussian sigma estimate',...
                                   'Standard error of x-position estimate','Standard error of y-position estimate',...
                                   'Standard error of intensity estimate', 'Standard error of background intensity estimate',...
                                   'Standard error of apparent gaussian sigma estimate',...
                                   'Frame index'};
    end

    properties
        %Parameters controlling [phase=3]: The identification of fetures representing possible point
        % emitters, buy looking for local maxima in a filtered image.
        ParamsFindMaxima=struct(... %phaseIdx=3
            'method','DoG',... %Options ['LoG'=Laplacian of Gaussian, 'DoG'=Difference of Gaussian'
            'filterSigmas',[1, 1.4, 1.9],... %In multiples of PSF
            'maximaNeighborhoodSize',7, ... %Options odd number >=3
            'scaleNeighborhoodSize',7 ... %Options odd numbner >=3
            );

        %Parameters controlling [phase=4]: The filtering of maxima to remove the background peaks and
        % form a set of candidate imitter boxes or sub-images that are to be fit in next stage.
        ParamsFilterMaxima=struct(...%phaseIdx=4
            'maximaThreshold',-1 ... %If empty or negative we will estimate this
            );

        %Parameters controlling [phase=5]: The localization (fitting) of the
        % identified candidate emitters from the boxes identified in phase 4
        ParamsLocalizeEmitters=struct(...%phaseIdx=5
            'model','Gauss2DsMAP',... %The class name of the emitter model to use
            'estimator','Newton'...   %The estimation technique to use
            );

        %Parameters controlling [phase=6]: The filtration of the fitted
        %emitters for quality.
        ParamsFilterEmitters=struct(...%phaseIdx=6
            'minIntensity',30,...
            'minSigma',0.70,... %Sigma is allowed to go as small as 0.5 pixels.
            'maxSigma',3.0,...
            'maxPositionSE',0.4,... %Maximum SE in either of the position arguments as computed by sqrt(crlb)
            'certVsUniformModel', -1,... % <1  0.95 = 95% certainty the emiitter model is correct;
            'certVsNoiseModel', -1, ...
            'overlapDistance',2.5 ... %Overlap 
            );
        
        %Parameters controlling [phase=7]: The tracking of the fitted
        %emitters.
        ParamsTrack=struct(...%phaseIdx=7
            'D',0.5,...          % pixels^2/frame
            'Kon',0.1,...      % 1/frame
            'Koff',0.1,...     % 1/frame
            'MaxSpeed',2.5,... % px/frame
            'MaxGapCloseFrames',20, ... %frames
            'MinGapCloseTrackLength',1,... 
            'MinFinalTrackLength',5 ...
            );  
    end

    properties (Dependent=true)       
        frameSize; %size [Y X] of an idivudal frame from the ROI
        nFrames;
        ROIOrigin; 
        ROIPhysical;
        ROIPhysicalOrigin;
    end
    
    properties (Hidden=true)
        version=1; %For future file format version changes
    end

    methods
        function obj=RPT( varargin )
            % Input (options):
            %  (1) <empty>
            %  (2) rpt_filepath
            %  (3) data_filepath, roi_in, roiname (optional)
            %  (4) data_obj, roi_in, roiname (optional)
            %  
            % roi_in: (optional) allows the choice of ROI to be saved selected.   The roi variable if given can be:
            %   (I) integer index into already created ROI (from data.ROI)
            %   (II) 2D - 1x4 integer array with [xmin, xmax, ymin, ymax]
            %        HS - 1x6 integer array with [xmin, xmax, ymin, ymax, Lmin, Lmax]
            %        This option uses the BaseData object's globalTBounds property to set the time bounds
            %                         ROIname will default to 'ManualROI'
            %   (III) 2D - 1x6 integer array with [xmin, xmax, ymin, ymax, tmin, tmax]
            %         HS - 1x8 integer array with [xmin, xmax, ymin, ymax, Lmin, Lmax, tmin, tmax]
            %                         ROIname will default to 'ManualROI'
            %   (IV) name string of already created ROI (from data.ROInames)
            %   (V) [] empty vector - auto generates a full frame ROI.
            %                    ROIname will default to 'FullFrameROI'
            % For case II, III, V an ROI will be generated with the name as specified for the particular case
            %
            % roiname: string (optional). This will force RPT.ROIname to roiname.
            %
            % Option (1) makes an empty RPT with no filename Result: [phase==1]
            % Option (2) Opens a saved RPT object from an saved file. Result: phase==saved object's phase
            % Option (3) Makes a new RPT object using the Data object stored in the file.  The new RPT 
            %            will be named by the spdata's ROIname.  The roi is optional and defaults to the 
            %            whole frame.  The roi_name is also optional and will be given a sensible value
            %            if ommitted.
            % Option (4) The same as option (4), except the SPData is given as an
            %            object instead of a filename.  
            %
            % These argurments are the same as the load method which can be used
            % to reload a different RPT into an already created object.
            % 
            if nargin>0
                obj.load(varargin{:});
            end
        end

        function stats=getStats(obj)
            stats=struct();
            if obj.phaseIdx >= 2
               stats.ROI = obj.ROI;
               stats.nFrames = obj.nFrames;
               stats.PixelSizeMicron = obj.data.pixelSize;
               stats.TimeStepSeconds = obj.data.frameT;
               stats.FrameSizePixels = flip(obj.frameSize);
               stats.FrameSizeMicron = flip(obj.frameSize)*obj.data.pixelSize;
            end
            if obj.phaseIdx >= 3
                stats.NumRawMaxima = length(obj.ResultsFindMaxima.rawMaxima);
            end
            if obj.phaseIdx >= 4
                stats.nMaxima = obj.nMaxima;
                stats.MeanMaximaPerFrame = length(obj.ResultsFilterMaxima.maxima)/obj.nFrames;
            end
            if obj.phaseIdx >= 6
                stats.nLocalizations = obj.nLocalizations;
            end
            if obj.phaseIdx >= 7
                stats.nTracks = obj.nTracks;
                stats.NSingletonTracks = sum(1==obj.ResultsTrack.trackLengths);
                stats.TrackMaxLength = max(obj.ResultsTrack.trackLengths);
            end

        end
        
        %% Parameter check functions
        % These will throw an error if something is wrong.  Also will check a new set of params
        % Still thinking if there is a better method, but this will get alot of cruft out of the actual
        % computation functions which should be tight.

        function P=checkParamsFindMaxima(obj,P)
            %
            if nargin==1
                P=obj.ParamsFindMaxima;
            end
        end

        function P=checkParamsFilterMaxima(obj,P)
            %
            if nargin==1
                P=obj.ParamsFilterMaxima;
            end
        end

        function P=checkParamsLocalizeEmitters(obj,P)
            %
            if nargin==1
                P=obj.ParamsLocalizeEmitters;
            end
        end

        function P=checkParamsFilterEmitters(obj,P)
            %
            if nargin==1
                P=obj.ParamsFilterEmitters;
            end
        end

        function P=checkParamsTrack(obj,P)
            %
            if nargin==1
               P=obj.ParamsTrack;
            end
        end

        function findMaxima(obj)
            obj.checkPhase(2);
            P = obj.checkParamsFindMaxima();
            obj.updateWaitbar(0,'Phase: Find Maxima');
            tic;

            obj.initializeBoxxer(true); %force reset            
            fframes = obj.getFilteredFrames();
            if obj.nScales == 1
                [R.rawMaxima, R.rawMaximaVals] = obj.boxxer.enumerateImageMaxima(fframes, P.maximaNeighborhoodSize);
                R.filteredSumImage = sumImage2D(fframes(:,:,1)); %use only first scale to make filtered image
            else
                switch P.method
                    case 'LoG'
                        [R.rawMaxima, R.rawMaximaVals] = obj.boxxer.scaleSpaceLoGMaxima(obj.getFrames(), P.maximaNeighborhoodSize, P.scaleNeighborhoodSize);
                    case 'DoG'
                        [R.rawMaxima, R.rawMaximaVals] = obj.boxxer.scaleSpaceDoGMaxima(obj.getFrames(), P.maximaNeighborhoodSize, P.scaleNeighborhoodSize);
                    otherwise
                        error('RPT:findMaxima','Unknown filter method "%s"',P.method);
                end
                R.rawMaxima([3,4],:) = R.rawMaxima([4,3],:); % flip [x y s t] to [x y t s] coords
                R.filteredSumImage = sumImage2D(fframes(:,:,:,1)); %use only first scale to make filtered image
            end
            R.rawMaximaImage = obj.computeMaximaImage(R.rawMaxima, R.rawMaximaVals);
            obj.ResultsFindMaxima = R;

            obj.updateWaitbar(1);
            obj.times.findMaxima = toc;
            fprintf('Find Maxima Time: %.3fs\n',obj.times.findMaxima);
            obj.setPhase(3);
        end
                
        function filterMaxima(obj)
            obj.checkPhase(3); 
            P = obj.checkParamsFilterMaxima();
            obj.updateWaitbar(0,'Phase: Filter Maxima');
            tic;

            Maxima = obj.ResultsFindMaxima;
            if P.maximaThreshold<=0; %auto threshold
                smax = sort(Maxima.rawMaximaVals,1,'descend');
                [~,R.maximaThreshold] = triThres(smax);
            else
                R.maximaThreshold = P.maximaThreshold;
            end
            obj.updateWaitbar(0.2);

            R.filter = Maxima.rawMaximaVals>=R.maximaThreshold;
            R.maxima = Maxima.rawMaxima(:,R.filter);
            R.maximaVals = Maxima.rawMaximaVals(R.filter);
            R.maximaImage = obj.computeMaximaImage(R.maxima, R.maximaVals);
            obj.updateWaitbar(0.4);

            % boxCoords is a BoxCoords object which organizes all the information about
            % a group of boxes including scale and frame information.
            obj.initializeBoxxer()
            R.boxCoords = obj.boxxer.generateBoxCoords(R.maxima, obj.nFrames);
            obj.updateWaitbar(0.7);
            
            [R.emitterImages, R.emitterFrames] = R.boxCoords.makeROI(obj.getFrames());
            obj.ResultsFilterMaxima = R;

            obj.updateWaitbar(1);
            obj.times.filterMaxima = toc;
            fprintf('Filter Maxima Time: %.3fs\n',obj.times.filterMaxima);
            obj.setPhase(4);
        end
        
        function localizeEmitters(obj)
            obj.checkPhase(4); 
            P = obj.checkParamsLocalizeEmitters();
            obj.updateWaitbar(0,'Phase: Localize Emitters');
            tic;

            FMaxima = obj.ResultsFilterMaxima;
            im_list = FMaxima.emitterImages;
            boxCoords = obj.getBoxCoords();
            obj.initializeEmitterModel();

            %For each box size category we need to gather all the ROI and fit them
            %with one of the n boxCoords emitter models.
            N = boxCoords.NsizeCategories;
            etheta= cell(1,N);
            crlb = cell(1,N);
            llh = cell(1,N);
            thetaInit = cell(1,N);
            frameIdx = cell(1,N);
            boxIdx = cell(1,N);
            for n=1:N
                idxs = boxCoords.sizeIndexes{n};
                boxIdx{n} = idxs;
                frameIdx{n} = boxCoords.boxFrameIdx(idxs);
                scales = boxCoords.scaleSigmas(1,boxCoords.boxScaleIdx(idxs));
                positions = 0.5+double(boxCoords.boxMaxima(:,idxs) - boxCoords.boxOrigin(:,idxs));
                thetaInit{n} = double([positions; zeros(2,numel(idxs)); scales]);
                if obj.emitterModel{n}.nParams==4
                    thetaInit{n} = thetaInit{n}(1:4,:); %For the 4-parameter models don't include the sigma-scales.
                end
%                 [etheta{n}, crlb{n}, llh{n}] = obj.emitterModel{n}.estimate(im_list{n}, P.estimator, thetaInit{n});
                if isempty(im_list{n})
                    etheta{n}=[];
                    crlb{n}=[];
                    llh{n}=[];
                else
                    [etheta{n}, crlb{n}, llh{n}] = obj.emitterModel{n}.estimate(im_list{n}, P.estimator);
                end
                obj.updateWaitbar(0.1+0.8*n/N);
            end
            obj.ResultsLocalizeEmitters.rawTheta = etheta; %rawTheta is in cell-based format for easy re-use with the emitter Model
            obj.ResultsLocalizeEmitters.thetaInit = thetaInit; %save the theta init for later debugging.

            %Form localizations
            etheta = [etheta{:}];
            crlb = [crlb{:}];
            llh = vertcat(llh{:})';
            boxIdx = double([boxIdx{:}]);
            frameIdx = double([frameIdx{:}]);
            E = [etheta; sqrt(crlb); llh; boxIdx; frameIdx]'; %internal emitter format
            shift = double(boxCoords.boxOrigin([2,1],boxIdx)')-1; %shift switches x/y since boxxer deals in row/col and locs are in x/y
            E(:,1:2) = E(:,1:2) + shift + repmat(obj.ROIOrigin,size(E,1),1); %correct for box coords

            %Sorted emitters by frame Idx to maintain relationship with localizations
            [~,sidx] = sort(E(:,end));
            E = E(sidx,:);

            obj.ResultsLocalizeEmitters.rawEmitters = E;            

            obj.updateWaitbar(1);
            obj.times.localizeEmitters = toc;
            fprintf('Localize Emitters Time: %.3fs\n',obj.times.localizeEmitters);
            obj.setPhase(5);
        end
        
        function filterEmitters(obj)
            obj.checkPhase(5);
            P = obj.checkParamsFilterEmitters();
            R = obj.ResultsLocalizeEmitters;
            ims = obj.ResultsFilterMaxima.emitterImages;
            obj.updateWaitbar(0,'Phase: Filter Emitters');
            tic;
            
            E = obj.getRawEmitters;
            boxes = obj.getBoxCoords(); % the box coords object
            N = size(E,1); % number of raw localziations            
            filter = zeros(N,1); %true if we will keep the localizations
            filterDescriptions = {'MinIntensity', 'MinSigma', 'MaxSigma', 'MaxPositionSE', 'UniformModelComparison', 'NoiseModelComparison','OverlapDistance'};
            filterStepCode = 1; % increment this after each test.  This indicated which filter method was repsonsible for filtering point
            if ~isempty(P.minIntensity)
                filter(E(:,3)<P.minIntensity) = filterStepCode;
            end
            filterStepCode = filterStepCode+1;
            if ~isempty(P.minSigma) && P.minSigma > 0
                filter(~filter & E(:,5)<P.minSigma) = filterStepCode;
            end
            filterStepCode = filterStepCode+1;
            if ~isempty(P.maxSigma) && P.maxSigma > 0
                filter(~filter & E(:,5)>P.maxSigma) = filterStepCode;
            end
            filterStepCode = filterStepCode+1;
            if ~isempty(P.maxPositionSE) && P.maxPositionSE > 0
                filter(~filter & (E(:,6)>P.maxPositionSE | E(:,7)>P.maxPositionSE)) = filterStepCode;
            end
            filterStepCode = filterStepCode+1;
            if ~isempty(P.certVsUniformModel) && P.certVsUniformModel>=0
                pass = cell(obj.nScales,1);
                for n=1:obj.nScales
                    pass{n} = obj.emitterModel{n}.modelComparisonUniform(P.certVsUniformModel, ims{n}, R.rawTheta{n});
                end
                pass = cell2mat(pass);
                filter(~pass)=filterStepCode;
            end
            filterStepCode = filterStepCode+1;
            if ~isempty(P.certVsNoiseModel) && P.certVsNoiseModel>=0
                pass = cell(obj.nScales,1);
                for n=1:obj.nScales
                    pass{n} = obj.emitterModel{n}.modelComparisonNoise(P.certVsNoiseModel, ims{n}, R.rawTheta{n});
                end
                pass = cell2mat(pass);
                filter(~pass)=filterStepCode;
            end
            filterStepCode = filterStepCode+1;
            if ~isempty(P.overlapDistance) && P.overlapDistance > 0              
                for frame = 1:obj.nFrames
                    fidx = boxes.frameIndexes{frame};
                    fidx = fidx(~filter(fidx)); % Select only non-filtered points.
                    frameE = E(fidx,:);
                    nFrameE = numel(fidx);
                    if nFrameE > 1
%                         tooClose = false(nFrameE,1);
                        for n = 1:nFrameE-1
                            dists = sqrt(sum((repmat(frameE(n,[1,2]),nFrameE-n,1)-frameE(n+1:end,[1,2])).^2,2));                            
                            nearby = [n; n+find(dists <= P.overlapDistance)];
                            if numel(nearby)>1
                                nearby = fidx(nearby);
                                posSE = sum(E(nearby,6:7),2);
                                [~,sidx] = sort(posSE);
                                discard = nearby(sidx(2:end));
                                filter(discard) = filterStepCode;
                            end
                        end
                    end
                end
            end

            R.filter = filter;
            R.filterDescriptions = filterDescriptions;
            R.emitters = E(~filter,:);
            
            %Use the filter to make the rawTheta for each sizeCategory into a selected theta within the
            %size category
            R.theta = cellmap(@(k) obj.ResultsLocalizeEmitters.rawTheta{k}(:,~filter(boxes.scaleIndexes{k})), 1:obj.nScales);

            obj.ResultsFilterEmitters = R;

            obj.updateWaitbar(1);
            obj.times.filterEmitters = toc;
            fprintf('Filter Emitters Time: %.3fs\n',obj.times.filterEmitters);
            obj.setPhase(6);
        end
        
        function trackEmitters(obj)
            obj.checkPhase(6);
            P = obj.checkParamsTrack();
            obj.updateWaitbar(0,'Phase: Track Emitters');
            tic;

            
            E = obj.getEmitters();
            position = E(:,[1,2]);
            SE_position = E(:,[6,7]);
            frameIdx = E(:,13)+obj.ROI(5)-1;
            obj.tracker = BaseTrack(obj.ROI, frameIdx, position, SE_position);
            obj.tracker.D = P.D;
            obj.tracker.kon = P.Kon;
            obj.tracker.koff = P.Koff;
            obj.tracker.maxSpeed = P.MaxSpeed;
            obj.tracker.maxGapCloseFrames = P.MaxGapCloseFrames;
            obj.tracker.minGapCloseTrackLength = P.MinGapCloseTrackLength;
            obj.tracker.minFinalTrackLength = P.MinFinalTrackLength;

            obj.tracker.doLAP();
            obj.updateWaitbar(0.9,'Track Emitters');

            R.tracks = obj.tracker.makeTracksArray(obj.getLocalizations());
            %Sort by the track length
            R.trackLengths = cellfun(@(t) size(t,1), R.tracks);
            [R.trackLengths, sidx] = sort(R.trackLengths,2,'descend');
            R.tracks = R.tracks(sidx);

            obj.ResultsTrack = R;

            obj.updateWaitbar(1);
            obj.times.trackEmitters = toc;
            fprintf('Track Emitters Time: %.3fs\n',obj.times.trackEmitters);
            obj.setPhase(7);
        end

        %% Data Retrieval Methods
        function [L, col_names, col_units, col_desc] = getLocalizations(obj)
            % Localizations are sorted by frame
            % [out]
            %   L - Localizations in matrix format with column giving parameters rows giving localizations.  
            %    col_names - [optional] the value of obj.LocalizationColumnNames for convenience
            %    col_units - [optional] the value of obj.LocalizationColumnUnits for convenience
            %    col_desc  - [optional] the value of obj.LocalizationColumnDescriptions for convenience
            E = obj.getEmitters();
            % For now only handle data with uniform pixelSize
            assert(isscalar(obj.data.pixelSize) || obj.data.pixelSize(1) == obj.data.pixelSize(2));
            sz = obj.data.pixelSize(1);
            frameT = obj.data.frameT;
            frameIdx = E(:,13)+obj.ROI(5)-1;
            L = [frameIdx.*frameT,...%t(s)
                 E(:,1:2).*sz,... %x(um) y(um)
                 E(:,3:4),...%I bg
                 E(:,5).*sz,...%sigma(um)
                 E(:,6:7).*sz,... % SE_x(um) SE_y(um)
                 E(:,8:9),... % SE_I SE_bg
                 E(:,10).*sz,... % SE_sigma(um)
                 frameIdx]; %frameIdx
            %sort by frame index
            [~,sidx] = sort(L(:,end));
            L = L(sidx,:);
            if nargout>1 % Provide these for convenience only
                col_names = obj.LocalizationColumnNames;
                col_units = obj.LocalizationColumnUnits;
                col_desc = obj.LocalizationColumnDescriptions;
            end
        end
        
        function [Ts, col_names, col_units, col_desc] = getTracksPixels(obj,trackIdxs)
            % The RPT track format is defined by the RPT constant properties:
            %  TrackColumnNames, TrackColumnUnits, and TrackColumnDescriptions
            %
            % *NOTE*: This method uses pixels units instead of microns to report back positions and positionSE.  This
            % is a convenience method provided for those who wish to work in pixel units.  Note that this            
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
            %     col_units - [optional] The units of each column.  This is in pixels where appropriate which
            %                            differs from obj.TrackColumnUnits
            %     col_desc  - [optional] the value of obj.TrackColumnDescriptions for convenience
            obj.checkPhase(7);
            Ts = obj.ResultsTrack.tracks;
            if nargin==2
                Ts = Ts(trackIdxs);%select tracks if requested
            end
            for n=1:numel(Ts)
                % Convert x,y,sigma, SEx, SEy, SEsigma to pixels
                Ts{n}(:,[2,3,6,7,8,11])=Ts{n}(:,[2,3,6,7,8,11])./obj.data.pixelSize;
            end
            if nargout>1 % Provide these for convenience only
                col_names = obj.TrackColumnNames;
                col_units = obj.TrackColumnUnits;
                col_units([2,3,6,7,8,11]) = 'px'; % Change units to pixels
                col_desc = obj.TrackColumnDescriptions;
            end
        end
        
        %% Testing Methods

        function testFit(obj, idx)
            obj.checkPhase(5); 
            L = obj.ResultsLocalizeEmitters.rawLocalizations;
            boxCoords = obj.emitterBoxCoords;
            loc = L(idx,:);
            boxIdx = loc(end-1);
%             loc(1:2) = loc(1:2) - double(boxCoords.boxOrigin(:,loc(end-1)'))-1;
            roi_im = boxCoords.makeROIsingle(obj.getFrames(),boxIdx);
            sizeIdx = boxCoords.boxSizeCategoryIdx(boxIdx);
            model = obj.emitterModel{sizeIdx};
            position = 0.5+double(boxCoords.boxMaxima(:,boxIdx)-boxCoords.boxOrigin(:,boxIdx));
            scale = boxCoords.scaleSigmas(1,boxCoords.boxScaleIdx(boxIdx));
            theta_init =  double([position', 0, 0, scale]);
            obj.viewMaximizedDipFig(roi_im);
            [etheta, crlb, llh, stats, seq, seq_llh] = model.estimateMAPDebug(roi_im, 'Newton',theta_init)  
            fit_im = model.modelImage(etheta);
            obj.viewMaximizedDipFig(fit_im);
        end

       
        %% Plotting and Visualization
        % All methods make a new figure if called from without axes
        function plotSumImage(obj,axH)
            obj.checkPhase(2);
            if nargin==1
                figure();
                axH = axes();
            else
                axes(axH);
            end
            BaseRPT.plot_imagesc(axH,obj.sumImage,obj.ROIOrigin)
        end

        function plotFilteredSumImage(obj,axH)
            obj.checkPhase(3);
            if nargin==1
                figure(); 
                axH=axes();
            else
                axes(axH);
            end
            BaseRPT.plot_imagesc(axH,obj.ResultsFindMaxima.filteredSumImage,obj.ROIOrigin)
        end


        function plotRawMaximaImage(obj,axH)
            obj.checkPhase(3);
            if nargin==1
                figure();
                axH=axes();
            else
                axes(axH);
            end
            BaseRPT.plot_surface(axH,obj.ResultsFindMaxima.rawMaximaImage,obj.ROIOrigin);
        end

        function plotRawMaximaPerFrame(obj,axH)
            obj.checkPhase(3);
            if nargin==1
                figure();
                axH=axes();
            end
            frame_series = {obj.ResultsFindMaxima.rawMaxima(3,:)};
            series_names = {'Raw Maxima'};
            obj.plotEventsPerFrame(axH,frame_series, series_names);
            title('Raw Maxima Frame Distribution');
        end

        function plotRawMaximaPerScale(obj,axH)
            obj.checkPhase(3);
            if nargin==1
                figure();
                axH=axes();
            end
            scale_series = {obj.ResultsFindMaxima.rawMaxima(4,:)};
            series_names = {'Raw Maxima'};
            obj.plotEventsPerScale(axH,scale_series, series_names);
            title('Raw Maxima Scale Distribution');
        end

        function plotThreshold(obj,axH)
            obj.checkPhase(3);
            if nargin==1
                figure();
                axH=axes();
            end
            Maxima = obj.ResultsFindMaxima;
            smax = sort(Maxima.rawMaximaVals,1,'descend');
            triThres(smax, axH);
        end

        function plotEmitterSumImage(obj,axH)
            obj.checkPhase(4);
            if nargin==1
                figure(); 
                axH=axes();
            else
                axes(axH);
            end
            R = obj.ResultsFilterMaxima;
            im = sumImage2D(Boxxer.plotEmitterMovie(R.boxCoords,R.emitterImages));
            BaseRPT.plot_imagesc(axH,im,obj.ROIOrigin)
        end
        
        function plotMaximaPerFrame(obj,axH)
            obj.checkPhase(4);
            if nargin==1
                figure();
                axH=axes();
            end
            frame_series = {obj.ResultsFindMaxima.rawMaxima(3,:), ...
                            obj.ResultsFilterMaxima.maxima(3,:)};
            series_names = {'Raw Maxima', 'Filtered Maxima'};
            obj.plotEventsPerFrame(axH,frame_series, series_names);
            title('Filtered Maxima Frame Distribution');
        end
 
        function plotMaximaPerScale(obj,axH)
            obj.checkPhase(4);
            if nargin==1
                figure();
                axH=axes();
            end
            scale_series = {obj.ResultsFilterMaxima.maxima(4,:)};
            series_names = {'Filtered Maxima'};
            obj.plotEventsPerScale(axH,scale_series, series_names);
            title('Filtered Maxima Scale Distribution');
        end

        function plotMaximaImage(obj,axH)
            obj.checkPhase(4);
            if nargin==1
                figure();
                axH=axes();
            else
                axes(axH);
            end
            BaseRPT.plot_surface(axH,obj.ResultsFilterMaxima.maximaImage,obj.ROIOrigin);
        end

        function plotBoxesImage(obj,axH)
            obj.checkPhase(4);
            if nargin==1
                figure();
                axH=axes();
            else
                axes(axH);
            end
            obj.updateWaitbar(0.1,'Generating Boxes 2D');
            boxCoords = obj.getBoxCoords();
            im = single(boxCoords.plotBoxes(obj.getFrames()));
            sumim = sqrt(1+mean(im,3));
            im = joinchannels('RGB',sumim(:,:,:,1), sumim(:,:,:,2),sumim(:,:,:,3));
            axes(axH);
            image([.5,size(im,2)-.5]+obj.ROIOrigin(1),[.5,size(im,1)-.5]+obj.ROIOrigin(2),mat2gray(single(im)));
            xlabel('X (px)');
            ylabel('Y (px)');
            title('Boxes as a 2D image');
            obj.updateWaitbar(1);
        end
        
        %% Phase 5 LocalizeEmitters Plotting

        function plotRawLocalizationPosDist(obj,axH)
            obj.checkPhase(5);
            if nargin==1
                figure();
                axH=axes();
            else
                axes(axH);
            end
            boxCoords = obj.getBoxCoords();
            R=obj.getResults(5);
            distx = cell(1,boxCoords.NsizeCategories);
            disty = cell(1,boxCoords.NsizeCategories);
            for n=1:boxCoords.NsizeCategories
                sizes = double(boxCoords.sizeCategories(:,n))/2;
                distx{n} = R.rawTheta{n}(1,:)-sizes(1);
                disty{n} = R.rawTheta{n}(2,:)-sizes(2);
            end
            dists = {cell2mat(distx), cell2mat(disty)};
            names={'X','Y'};
            BaseRPT.plot_distributions(axH,dists, names);
            xlabel('Displacement from box center (px)');
            title('Raw Localization Poistion Distribution');
        end

        function plotRawLocalizationIDist(obj,axH)
            obj.checkPhase(5);
            if nargin==1
                figure();
                axH=axes();
            else
                axes(axH);
            end
            R=obj.getResults(5);
            boxCoords = obj.getBoxCoords;
            distI = cell(1,boxCoords.NsizeCategories);
            distbg = cell(1,boxCoords.NsizeCategories);
            for n=1:boxCoords.NsizeCategories
                sizes = double(boxCoords.sizeCategories(:,n));
                distI{n}  = R.rawTheta{n}(3,:);
                distbg{n} = R.rawTheta{n}(4,:)*prod(sizes);
            end
            dists = {cell2mat(distI), cell2mat(distbg)};
            names={'Intensity','Background/Box'};
            BaseRPT.plot_distributions(axH,dists, names);
            set(gca(),'XScale','log','YScale','linear');
            xlim([1,max(dists{1})]);
            xlabel('Intenisty (photons)');
            title('Raw Localization Intensity Distribution');
        end

        function plotRawLocalizationSigmaDist(obj,axH)
            obj.checkPhase(5);
            if nargin==1
                figure();
                axH=axes();
            else
                axes(axH);
            end
            R=obj.getResults(5);
            if size(R.rawTheta{1},1)<5
                error('RPT:ModelError','Emitters were fit with 4-parameter model that does not predict sigma');
            end
            dists = {cellmatfun(@(T) T(5,:), R.rawTheta)};
            names={'Sigma'};
            BaseRPT.plot_distributions(axH,dists, names);
            set(gca(),'YScale','log');
            xlim([0,max(cellfun(@max,dists))]);
            xlabel('Sigma (apparent/psf)');
            title('Raw Localization Sigma Distribution');
        end

        function plotRawLocalizationLLHDist(obj,axH)
            obj.checkPhase(5);
            if nargin==1
                figure();
                axH=axes();
            else
                axes(axH);
            end
            boxCoords = obj.getBoxCoords();
            E = obj.getRawEmitters();
            dists = cellmap(@(s) E(boxCoords.boxScaleIdx==s, end-2), 1:obj.nScales);
%             dists = {E(:,end-2)};
            names= cellmap(@(s) sprintf('LLH - Scale:%.2f', obj.ParamsFindMaxima.filterSigmas(s)), 1:obj.nScales);
            BaseRPT.plot_distributions(axH,dists, names);
            xlabel('Log Likelihood');
            title('Raw Localization Log Likelihood Distribution');
        end

        %% Phase 6 LocalizeEmitters Plotting
        function plotLocalizationSumImage(obj,axH)
            obj.checkPhase(6);
            if nargin==1
                figure();
                axH=axes();
            else
                axes(axH);
            end
            E = obj.getEmitters();
            BaseRPT.plot_imagesc(axH,obj.sumImage,obj.ROIOrigin);
            colormap('hot');
            hold('on');
            plot(E(:,1),E(:,2),'LineStyle','none','Marker','o','MarkerFaceColor',[0,1,0],...
                  'MarkerEdgeColor','none','MarkerSize',1.5);
        end

        function plotSelectedBoxesImage(obj,axH)
            obj.checkPhase(6);
            if nargin==1
                figure();
                axH=axes();
            else
                axes(axH);
            end
            obj.updateWaitbar(0.1,'Generating Boxes 2D');
            boxCoords = obj.getBoxCoords();
            filter = obj.ResultsFilterEmitters.filter;
            im = single(boxCoords.plotBoxes(obj.getFrames(),filter));
            sumim = sqrt(1+mean(im,3));
            im = joinchannels('RGB',sumim(:,:,:,1), sumim(:,:,:,2),sumim(:,:,:,3));
            axes(axH);
            image([.5,size(im,2)-.5]+obj.ROIOrigin(1),[.5,size(im,1)-.5]+obj.ROIOrigin(2),mat2gray(single(im)));
            xlabel('X (px)');
            ylabel('Y (px)');
            title('Boxes as a 2D image');
            obj.updateWaitbar(1);
        end

        function plotLocalizationPosDist(obj,axH)
            obj.checkPhase(6);
            if nargin==1
                figure();
                axH=axes();
            else
                axes(axH);
            end
            theta = obj.ResultsFilterEmitters.theta;
            boxCoords = obj.filteredEmitterBoxCoords;
            distx = cell(boxCoords.NsizeCategories,1);
            disty = cell(boxCoords.NsizeCategories,1);
            for n=1:boxCoords.NsizeCategories
                sizes = double(boxCoords.sizeCategories(:,n))/2;
                distx{n}  = theta{n}(1,:)'-sizes(2);
                disty{n} = theta{n}(2,:)'-sizes(1);
            end
            dists = {cell2mat(distx), cell2mat(disty)};
            names={'X','Y'};
            BaseRPT.plot_distributions(axH,dists, names);
            xlabel('Displacement from box center (px)');
            title('Filtered Localization Poistion Distribution');
        end

        function plotLocalizationIDist(obj,axH)
            obj.checkPhase(6);
            if nargin==1
                figure();
                axH=axes();
            else
                axes(axH);
            end
            theta = obj.ResultsFilterEmitters.theta;
            boxCoords = obj.filteredEmitterBoxCoords;
            distI = cell(boxCoords.NsizeCategories,1);
            distbg = cell(boxCoords.NsizeCategories,1);
            for n=1:boxCoords.NsizeCategories
                sizes = double(boxCoords.sizeCategories(:,n));
                distI{n}  = theta{n}(3,:)';
                distbg{n} = theta{n}(4,:)'*prod(sizes);
            end
            dists = {cell2mat(distI), cell2mat(distbg)};
            names={'Intensity','Background/Box'};
            BaseRPT.plot_distributions(axH,dists, names);
            set(gca(),'XScale','log','YScale','linear');
            xlabel('Intenisty (photons)');
            title('Localization Intensity Distribution');
        end

        function plotLocalizationSigmaDist(obj,axH)
            obj.checkPhase(6);
            if nargin==1
                figure();
                axH=axes();
            else
                axes(axH);
            end
            R=obj.getResults(6);
            dists = {cellmatfun(@(T) T(5,:), R.theta)};
            names={'Sigma'};
            BaseRPT.plot_distributions(axH,dists, names);
            xlabel('Sigma (apparent/psf)');
            title('Localization Sigma Distribution');
        end

        function plotLocalizationLLHDist(obj,axH)
            obj.checkPhase(6);
            if nargin==1
                figure();
                axH=axes();
            else
                axes(axH);
            end
            E = obj.getEmitters();
            dists = {E(:,end-2)};
            names = {'LLH'};
            BaseRPT.plot_distributions(axH,dists, names);
            xlabel('Log Likelihood');
            title('Raw Localization Log Likelihood Distribution');
        end


        function plot3DTrackSequence(obj, opts)
            %  opts [optional] - If provided this must be 2rd argument, so you must provide trackIds also.  
            %     It is a strcut which allow for setting of several options.
            %
            Ts = obj.getTracks();            
            Cs = cellmap(@(i) i*ones(size(Ts{i},1),1), 1:length(Ts));
            opts.cLabel = 'TrackID';
            opts.trackColorMap = @prism;
            opts.trackColorRange = max(1,length(Ts));
            obj.plotTracks3D(Ts, Cs, obj.ROIPhysical, obj.sumImage, opts);
        end

        function plot3DTrackTemporal(obj,opts)
            Ts = obj.getTracks();            
            Cs = cellmap(@(T) T(:,1), Ts);
            opts.cLabel= 'Time (s)';
            obj.plotTracks3D(Ts, Cs, obj.ROIPhysical, obj.sumImage, opts);
        end
        
        
        function plot3DTrackSpeed(obj,opts)
            if nargin==1 || ~isfield(opts,'winsize')
                opts.winsize= 6;
            end
            Ts = obj.getTracks();            
            Cs = cellmap(@(T) TrackSegmentAnalysis.estimateSpeedWindow(T, opts.winsize), Ts);          
            opts.cLabel = 'Speed ($\mu\mathrm{m}/\mathrm{s}$)';
            obj.plotTracks3D(Ts, Cs, obj.ROIPhysical, obj.sumImage, opts);
        end

        %% Phase 3 Views

        function f = viewRawMaximaMovie(obj)
            obj.checkPhase(3);
            obj.updateWaitbar(0.1,'Generating Raw Maxima Movie');
            R = obj.ResultsFindMaxima;            
            if obj.nScales == 1
                im = Boxxer.plotMaxima(obj.getFrames(), R.rawMaxima);
            else  
                im = Boxxer.plotScaleMaximaRPT(obj.getFrames(), R.rawMaxima);
            end
            f = obj.viewMaximizedDipFig(im);
            f.Name = 'Raw Maxima over Frames';
            obj.updateWaitbar(1);
        end

        function f = viewRawMaximaFilteredMovie(obj)
            obj.checkPhase(3);
            obj.updateWaitbar(0.1,'Generating Raw Maxima Movie');
            R = obj.ResultsFindMaxima;
            if obj.nScales == 1
                im = Boxxer.plotMaxima(obj.getFilteredFrames(), R.rawMaxima);
            else  
                im = Boxxer.plotScaleMaximaRPT(obj.getFilteredFrames(), R.rawMaxima);
            end
            f = obj.viewMaximizedDipFig(im);
            f.Name = 'Raw Maxima over Filtered Frames';    
            obj.updateWaitbar(1);
        end
        
        %% Phase 4 Views
        
        function f = viewMaximaMovie(obj)
            obj.checkPhase(4);
            obj.updateWaitbar(0.1,'Generating Maxima Movie');
            R = obj.ResultsFindMaxima;
            F = obj.ResultsFilterMaxima;
            if obj.nScales == 1
                im = Boxxer.plotMaxima(obj.getFrames(), R.rawMaxima, F.filter);
            else  
                im = Boxxer.plotScaleMaximaRPT(obj.getFrames(), R.rawMaxima);
            end
            f = obj.viewMaximizedDipFig(im);
            f.Name = 'Selected Maxima over Frames'; 
            obj.updateWaitbar(1);
        end

        function f = viewMaximaFilteredMovie(obj)
            obj.checkPhase(4);
            obj.updateWaitbar(0.1,'Generating Maxima Movie');
            R = obj.ResultsFindMaxima;
            F = obj.ResultsFilterMaxima;
            if obj.nScales == 1
                im = Boxxer.plotMaxima(obj.getFilteredFrames(), R.rawMaxima, F.filter);
            else
                im = Boxxer.plotScaleMaximaRPT(obj.getFilteredFrames(), R.rawMaxima, F.filter);
            end
            f = obj.viewMaximizedDipFig(im);
            f.Name = 'Selected Maxima over Filtered Frames';   
            obj.updateWaitbar(1);
        end

        function f = viewBoxesMovie(obj)
            obj.checkPhase(4);
            obj.updateWaitbar(0.1,'Generating Boxes Movie');
            boxCoords = obj.getBoxCoords();
            im = boxCoords.plotBoxes(obj.getFrames());
            f = obj.viewMaximizedDipFig(im);
            f.Name = 'Boxed Maxima over Frames';
            obj.updateWaitbar(1);
        end

        function f = viewBoxesFilteredMovie(obj)
            obj.checkPhase(4);
            obj.updateWaitbar(0.1,'Generating Boxes Movie');
            R = obj.ResultsFilterMaxima;           
            im = Boxxer.plotBoxCoordsDIP(obj.getFilteredFrames(), obj.getBoxCoords(),  R.maxima, R.maximaVals);
            f = obj.viewMaximizedDipFig(im);
            f.Name = 'Boxed Maxima over Filtered Frames';
            obj.updateWaitbar(1);
        end

        function f = viewEmitterMovie(obj)
            obj.checkPhase(4);
            obj.updateWaitbar(0.1,'Generating Emitter Movie');
            R = obj.ResultsFilterMaxima;
            im = Boxxer.plotEmitterMovie(obj.getBoxCoords(), R.emitterImages);
            f = obj.viewMaximizedDipFig(dip_image(im));
            dipmapping(f, 'global');
            %colormap('hot');
            f.Name = 'Emitters ROIs as a Movie';
            obj.updateWaitbar(1);
        end
        
        function f = viewEmitterImages(obj)
            obj.checkPhase(4);
            obj.updateWaitbar(0.1,'Generating Emitter Images');
            R = obj.ResultsFilterMaxima;
            imcube = cosmicNorm(R.boxCoords.makeROIcube(obj.getFrames()));
            [~,sidx] = sort(R.maximaVals,1,'descend');
            f = obj.viewMaximizedDipFig(dip_image(imcube(:,:,sidx)));
            dipmapping(f, 'global');
            colormap('hot');
            f.Name = 'Boxed Emitters indexed by Maxima Value';
            obj.updateWaitbar(1);
        end
 
        function f=viewBoxesSlices(obj)
            obj.checkPhase(4);
            f = figure('Name','Identified Candidate Emitters');
            obj.plotBoxesSlices(axes());
        end

        function plotBoxesSlices(obj, axH)
            obj.checkPhase(4);
            obj.updateWaitbar(0.1,'Generating Boxes 3D');
            [im,cm]=Boxxer.plotBoxCoordsOverlay(cosmicNorm(obj.getFrames()), obj.ResultsFilterMaxima.boxCoords);
            roi=obj.ROI;
            xs=((roi(1):roi(2)+1)-1)*obj.data.pixelSize;
            ys=((roi(3):roi(4)+1)-1)*obj.data.pixelSize;
            ts=((roi(5):roi(6)+1)-1)*obj.data.frameT;
            thin=ceil((roi(6)-roi(5)+1)/500);
            ts=ts(1:thin:end);
            [X,Y,Z]=meshgrid(xs,ys,ts);
            C=zeros(length(ys),length(xs),length(ts));
            im_thin=im(:,:,1:thin:end); 
            C(1:end-1,1:end-1,1:size(im_thin,3))=im_thin;
            obj.updateWaitbar(0.2);
            axes(axH);
            set(axH,'SortMethod','childorder');
            hold('on');
            colormap(cm);
            for t=1:length(ts)
                surface('XData',X(:,:,t),...
                        'YData',Y(:,:,t),...
                        'ZData',Z(:,:,t),...
                        'CData',C(:,:,t),...
                        'AlphaData',C(:,:,t),...
                        'FaceColor','flat','FaceAlpha',0.5,'EdgeColor','none',...
                        'FaceAlpha','flat','AlphaDataMapping','scaled');
                if mod(t,10)==0
                    obj.updateWaitbar(0.2+0.7*t/length(ts));
                end
            end
            axis([xs(1) xs(end) ys(1) ys(end) ts(1) ts(end)]);
            x_stretch=(xs(end)-xs(1))/(length(xs)-1);
            y_stretch=(ys(end)-ys(1))/(length(ys)-1);
            t_stretch=max(x_stretch,y_stretch);
            daspect([x_stretch y_stretch t_stretch]);
            view(135,54);
            hold('off');
            obj.updateWaitbar(0.9);
            set(axH,'Box','on','BoxStyle','full');
            set(axH,'XGrid','on','XMinorGrid','on','XMinorTick','on',...
                   'YGrid','on','YMinorGrid','on','YMinorTick','on',...
                   'ZGrid','on','ZMinorGrid','on','ZMinorTick','on');            
            xlabel('x ($\mu$m)','interpreter','latex');
            ylabel('y ($\mu$m)','interpreter','latex');
            zlabel('Time (s)','interpreter','latex');
            obj.updateWaitbar(1);
        end
    
        %% Phase 5 Views
        function f = viewRawLocalizationMovie(obj)
            obj.checkPhase(5);
            obj.updateWaitbar(0.1,'Generating Raw Localizations Movie');
            ims = cellmap(@(model, theta) model.modelImage(theta), obj.emitterModel, obj.ResultsLocalizeEmitters.rawTheta);
            f = obj.viewMaximizedDipFig(Boxxer.plotEmitterMovie(obj.getBoxCoords(), ims));
            dipmapping(f, 'global');
            f.Name = 'Emitters Model Fits as a Movie';
            obj.updateWaitbar(1);
        end
        
        function f = viewRawEmitterSuperResGauss(obj)
            obj.checkPhase(5);
            obj.updateWaitbar(0.1,'Generating Super Res Image');
            imSizePx = 8192;
            im = obj.makeEmitterSuperResGauss(obj.getRawEmitters(), imSizePx);
            obj.updateWaitbar(0.9,'Generating Super Res Image');
            f = obj.viewMaximizedDipFig(dip_image(im));
            colormap(hot);
            dipmapping(f,[0, prctile(im(:),99.9)]); %scale the color mapping to keep it from being too dim
            f.Name = 'Emitters Super-res Gaussian';
            obj.updateWaitbar(1);
        end            

        function f = viewEmitterModelComparison(obj)
            obj.checkPhase(5);
            obj.updateWaitbar(0.1,'Generating Fit Comparison');
            emitter_ims = obj.ResultsFilterMaxima.emitterImages;
            imsizes = cellmatfun(@(model) model.imsize', obj.emitterModel);
            maximsize = max(imsizes,[],2);
            ims = cell(obj.nScales,1);
            Mims = cellmap(@(model, theta) model.modelImage(theta), obj.emitterModel, obj.ResultsLocalizeEmitters.rawTheta);
            for n = 1:obj.nScales
                theta = obj.ResultsLocalizeEmitters.rawTheta{n};
                [~,sidx] = sort(theta(3,:),2,'descend');
                ims{n} = zeros(maximsize(1), 2*maximsize(2),numel(sidx));
                ims{n}(1:imsizes(1,n), 1:imsizes(2,n), :) = emitter_ims{n}(:,:,sidx);
                ims{n}(1:imsizes(1,n), maximsize(2): maximsize(2)+imsizes(2,n)-1, :) = Mims{n};
            end
            ims = dip_image(cosmicNorm(cat(3,ims{:})));
            f = obj.viewMaximizedDipFig(ims);
            dipmapping(f, 'global');
            colormap('hot');
            f.Name = 'Data vs Emitter Mdel Image';
            obj.updateWaitbar(1);
        end

        function f = viewRawLocaliztions3D(obj)
            obj.checkPhase(5);
        end

        %% Phase 6 Views      
        function f = viewEmitterSuperResGauss(obj)
            obj.checkPhase(6);
            if isempty(obj.srimage)
                obj.makeSRImage();
            end
            f = obj.viewMaximizedDipFig(dip_image(obj.srimage));
            colormap(hot);
            f.Name = 'Emitters Super-res Gaussian';
            obj.updateWaitbar(1);
        end            

        function plotEmitterSuperResGauss(obj,axH)
            if isempty(obj.srimage)
                obj.makeSRImage();
            end
            BaseRPT.plot_imagesc(axH,obj.srimage,obj.ROIOrigin);
        end

        function makeSRImage(obj)
            obj.updateWaitbar(0.1,'Generating Super Res Image');
            imSizePx = 8192;
            im = obj.makeEmitterSuperResGauss(obj.getEmitters(),imSizePx);
            im(:) = min(im(:), prctile(im(:),99.9));
            obj.srimage = im;
            obj.updateWaitbar(1);
        end

        function f = viewLocalizationMovie(obj)
            obj.checkPhase(6);
            obj.updateWaitbar(0.1,'Generating Localizations Movie');
            boxCoords = obj.filteredEmitterBoxCoords;
            ims = cellmap(@(model, theta) model.modelImage(theta), obj.emitterModel, obj.ResultsFilterEmitters.theta);
            obj.updateWaitbar(0.3,'Generating Localizations Movie');
            movie = Boxxer.plotEmitterMovie(boxCoords, ims);
            obj.updateWaitbar(0.9,'Generating Localizations Movie');
            f = obj.viewMaximizedDipFig(movie);            
            dipmapping(f, 'global');
            f.Name = 'Selected Emitter Model Fits as a Movie';
            obj.updateWaitbar(1);
        end
        
        function f = viewSelectedBoxesMovie(obj)
            obj.checkPhase(4);
            obj.updateWaitbar(0.1,'Generating Boxes Movie');
            boxCoords = obj.getBoxCoords();
            filter = obj.ResultsFilterEmitters.filter;
            filterDescriptions = ['Accepted', obj.ResultsFilterEmitters.filterDescriptions];
            im = boxCoords.plotBoxes(obj.getFrames(),filter,filterDescriptions);
            f = obj.viewMaximizedDipFig(im);
            f.Name = 'Boxed Maxima over Frames';
            obj.updateWaitbar(1);
        end

        function f = view3DTrackSequence(obj, varargin)
            obj.checkPhase(7);
            f=figure('Name','Tracks 3D Sequence');
            whitebg(f);
            obj.plot3DTrackSequence(varargin{:})
        end

        function f = view3DTrackSpeed(obj, varargin)
            obj.checkPhase(7);
            f=figure('Name','Tracks 3D Speed');
            whitebg(f);
            obj.plot3DTrackSpeed(varargin{:})
        end

        function f = view3DTrackTemporal(obj, varargin)
            obj.checkPhase(7);
            f=figure('Name','Tracks 3D Temporal');
            whitebg(f);
            obj.plot3DTrackTemporal(varargin{:})
        end

        function f = viewTrackMovie(obj, varargin)
            % [in] trackIds - [optional] list of track ids to add to movie
            obj.checkPhase(7);
            tm= TrackMovie(obj, varargin{:});
            tm.setTrackColorMethod('Sequence');
            f = tm.viewSequence();
            f.Name='Tracks Movie';            
        end
    end % Public methods

        
    %% Dependent properties
    methods
        function val = get.frameSize(obj)
            % Frame size [Y X] in pixels
            if isempty(obj.ROI)
                val = [];
            else
                val = [obj.ROI(4)-obj.ROI(3)+1, obj.ROI(2)-obj.ROI(1)+1];
            end
        end
        
        function N = get.nFrames(obj)
            if isempty(obj.ROI)
                N = 0;
            else
                N = obj.ROI(6)-obj.ROI(5)+1;
            end
        end
        
        function origin = get.ROIOrigin(obj)
            % shift [X,Y] in pixels that should be added to local ROI pixel
            % coords to translate into global ROI coords
            if isempty(obj.ROI)
                origin = [];
            else
                origin = obj.ROI([1,3])-1;
            end
        end
        
        function pROI = get.ROIPhysical(obj)
            % The ROI in physical coordinates of um and s instead of pixels
            % and frames
            if isempty(obj.ROI)
                pROI = [];
            else
                pS = obj.data.pixelSize;
                dT = obj.data.frameT;
                pROI = [(obj.ROI(1)-1)*pS, obj.ROI(2)*pS,...
                        (obj.ROI(3)-1)*pS, obj.ROI(4)*pS,...
                        (obj.ROI(5)-1)*dT, obj.ROI(6)*dT];
            end
        end
        
        function porigin = get.ROIPhysicalOrigin(obj)
            % shift [X,Y] in um that should be added to local physical
            % coords to translate into global physical coords
            if isempty(obj.ROI)
                porigin = [];
            else
                porigin = (obj.ROI([1,3])-1)*obj.data.pixelSize;
            end
        end
    end % Public methods

    methods (Static = true)
        function tableT = convertTracksTable(Ts, trackIds)
            % This is static method to preform the reformatting of the RPT tracks format into a single table 
            % where all tracks have been concatenated and distinguished by an extra trackID column.  
            % This is slower an more inefficient than the default RPT Tracks format, but can be convenient as Matlab
            % tables have the ability to store the names, units, and descritptions of all columns.
            % [in]
            %   Ts - A cellarray of track matricies in normal RPT Tracks format
            %   trackIds - [optional] if Ts was a subset of all tracks selected, the trackIds selected should
            %               also be provided to ensure that the first column can be correctly indexed. [default is
            %               a linear track indexing which will only match this RPT file if all tracks were
            %               selected]
            % [out] 
            %  tableT - The tracks converted into a table giving each tracks information using an inital column for trackID
            if nargin<2
                trackIds = 1:numel(Ts);
            end           
            tableT = array2table(cellcat(cellmap(@(id,t) [repmat(id,size(t,1),1), t], num2cell(trackIds), Ts)'));
            tableT.Properties.DimensionNames = RPT.TrackTableDimensionNames;
            tableT.Properties.Description = RPT.TrackTableTitle;
            %Add a Track index column to the table output.
            tableT.Properties.VariableNames = ['trackID',RPT.TrackColumnNames];
            tableT.Properties.VariableDescriptions = ['Track Index', RPT.TrackColumnDescriptions];
            tableT.Properties.VariableUnits = [{''}, RPT.TrackColumnUnits];
        end
        
        function structT = convertTracksStruct(Ts)
            % This is static method to preform the reformatting of the RPT tracks format into a structure array,
            % where each track is stored as a structure with field corresponding to the columns of the normal RPT
            % track format.
            % [in]
            %   Ts - A cellarray of track matricies in normal RPT Tracks format
            % [out] 
            %  tableT - The tracks converted into a structre array with one structure per track
            cols = cellmap(@(i) cellmap(@(t) t(:,i),Ts), 1:RPT.NTrackColumns);
            args = [RPT.TrackColumnNames; cols];
            structT = struct(args{:});
        end
        function D_mle = diffusionConstMLE_individual(Ts)
            % Compute the maximum likelihood estimate of D for each trajectory in Ts individually
            % [in]
            %   Ts - A cellarray of track matricies in normal RPT Tracks format
            % [out] 
            %  D_mle - [um^2/s] A vector of the MLE estimate of D for each trajectory.
            N=numel(Ts);
            D_mle=zeros(N,1);
            est = DEstimator();
            for k=1:N
                est.initializeTrack(Ts{k}(:,2:3),Ts{k}(:,1),Ts{k}(:,6:7));
                D_mle(k) = est.MLE();
            end
        end
        function [D_mle, mle_LLH, confInt] = diffusionConstMLE_ensemble(Ts,conf_alpha)
            % Compute the maximum likelihood estimate of D for all trajectorys together
            % [in]
            %   Ts - A cellarray of track matricies in normal RPT Tracks format
            %   conf_alpha - [optional] [default=0.05] the 100*(1-alpha)% confidence itervals to return.
            %                                              this should be 0<conf_alpha<1
            % [out] 
            %  D_mle - [um^2/s] A vector of the MLE estimate of D for all trajectories together
            %  mle_LLH - log likelihood at the D_mle
            %  confInt - [um^2/s] Lower and Upper bounds on D, for 100(1-conf_alpha)% confidence interval
            if nargin<2
                conf_alpha = 0.05;
            end
            [D_mle, mle_LLH, confInt] = DEstimator.computeEnsembleMLE(cellmap(@(T) T(:,[1,2,3,6,7]),Ts),[],conf_alpha);
        end
        

        function trackHs = plotTracks3D(Ts, Cs, ROIphysical, im, opts)
            % This operates in physical units of micron(um) and seconds.  Draws tracks using the surface
            % command which is OpenGL accelerated.
            %
            % [IN]
            % Ts - Tracks in RPT Track format - cell array of matricies of localizations
            % Cs - Color values for each track - cell array of 1D vectors of values the same length as 
            %      each track
            % ROIphysical - Bounding box of drawing in 3D space.  [xmin xmax ymin ymax tmin tmax]
            %             This is in phyiscal units of microns and seconds.
            % im - An image to display in the background in B/W, dimensions should conform to
            %         ROIphysical.
            % opts.trackColorMap [optional] - the function handle for a Tracks color map defaults to @jet
            % opts.imageColorMap [optional] - the function handle for an image color map defaults to @gray
            % opts.cLabel [optional] - The label for the colorbar.  If omitted no colorbar shown.
            % [OUT]
            %  trackHs - handles to the surface object correponding to each track.  Use these to modify
            %            the tracks afterwards with userdata/callbacks etc
            if nargin == 4
                opts=struct();
            end
            if isfield(opts,'trackColorMap')
                track_cmap = opts.trackColorMap;
            else
                track_cmap = @jet; 
            end
            if isfield(opts,'trackColorRange')
                track_cRange = opts.trackColorRange;
            else
                track_cRange = 256; %Number of different color values to display in tracks
            end
            if isfield(opts,'imageColorMap')
                image_cmap = opts.imageColorMap;
            else
                image_cmap = @gray;
            end
            if isfield(opts,'imageColorRange')
                image_cRange = opts.imageColorRange;
            else
                image_cRange = 256;%Number of different color values to display in bg image
            end
            if isfield(opts,'markerSize')
                markerSize = opts.markerSize;
            else
                markerSize = 2.5;%Size of marker at end of tracks
            end
            endPointColor = [1,0,1]; %Magenta
            singletonPointColor = [0,1,1]; %Cyan
            imAlpha = 1; %Alpha shading for bg image
            colormap([track_cmap(track_cRange);image_cmap(image_cRange)]);%Colormap is combined for track and image
            nTracks = numel(Ts);
            %Determine color range
            minCs = min(cellfun(@min,Cs));
            maxCs = max(minCs+1,max(cellfun(@max,Cs))); %assure span is at least 1
            spanCs = maxCs-minCs;
            if isempty(spanCs)
                minCs = 0;
                maxCs = 1;
                spanCs = 1;
            end

            %Plot Image
            im = im - min(im(:)); %shift to 0
            im = im ./ max(im(:)); %Scale to [0,1]
            im = (im.*spanCs.*(image_cRange/track_cRange))+maxCs*1.0001; %Shift image up to the top of the colorspace so it gets mapped to gray
            minT = min(cellfun(@(T) min(T(:,1)),Ts)); %Min time
            if isempty(Ts)
                minT = ROIphysical(5);
            end
            BX = repmat(ROIphysical(1:2),2,1);
            BY = repmat(ROIphysical(3:4)',1,2);
            BZ = repmat(minT,2,2);
            surface('XData',BX,'YData',BY,'ZData',BZ,'CData',im,...
                    'FaceColor','texturemap','EdgeColor','none','FaceAlpha',imAlpha);
            hold('on');
            set(gca,'YDir','reverse','TickDir','out');
            set(gca,'ZMinorTick','on','YMinorTick','on','XMinorTick','on',...
                    'ZGrid','on','XGrid','on','YGrid','on','Box','on','BoxStyle','full',...
                    'Projection','Orthographic');

            %Set aspect ratio of axes
            asp=pbaspect();
            asp(1) = (ROIphysical(2)-ROIphysical(1))/(ROIphysical(4)-ROIphysical(3));
            asp(2) = 1;
            asp(3) = max(asp(1),asp(2));
            pbaspect(asp);
            view(0,90);
            axis('tight');

        
            %Plot the tracks
            if nTracks>0
                for i = nTracks:-1:1
                    ts = Ts{i}(:,1)';
                    xs = Ts{i}(:,2)';
                    ys = Ts{i}(:,3)';
                    C = Cs{i}(:)';
                    if length(ts)>1
                        trackHs(i) = surface([xs;xs],[ys;ys],[ts;ts],[C;C],'EdgeColor','interp','LineWidth',1.0);
                        plot3([xs(1),xs(end)], [ys(1),ys(end)],[ts(1),ts(end)],'o','MarkerSize',markerSize,...
                                    'MarkerFaceColor',endPointColor,'MarkerEdgeColor',[0,0,0]);
                    else
                        plot3(xs(1), ys(1),ts(1),'o','MarkerSize',markerSize,'MarkerFaceColor',singletonPointColor,'MarkerEdgeColor',[0,0,0]);
                    end
                end
            end

            hold('off');
%             view(-37.5,30);
            %Fix grids and ticks
%             fixROIAxesTicks('XTick',ROIphysical(1), ROIphysical(2));
%             fixROIAxesTicks('YTick',ROIphysical(3), ROIphysical(4));
%             fixROIAxesTicks('ZTick',ROIphysical(5), ROIphysical(6));
            %Labels
            zlabel('t (s)','interpreter','latex');
            xlabel('x (um)','interpreter','latex');
            ylabel('y (um)','interpreter','latex');
            
            if isfield(opts,'cLabel') && ~isempty(opts.cLabel) %Do colorbar
                cbh = colorbar();
                cbh.YLim = [minCs maxCs];
                ylabel(cbh, opts.cLabel,'interpreter','latex');
            else
                colorbar('off');
            end

            function fixROIAxesTicks(prop,minval,maxval)
                %An inner function to fix the axes ticks
                ticks=get(gca,prop);
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
        end

        function checkTracks(tracks)
            %check the tracks are OK
            if ~iscell(tracks)
                error('RPT:checkTracks','Tracks must be a cell-array track format');
            elseif ~all(cellfun(@(t) size(t,2),tracks) == RPT.NTrackColumns)
                error('RPT:checkTracks','All tracks must have %i colums',RPT.NTrackColumns);
            end
        end

        function rpt_files=batchProcess(spdatapath, spdatafile_patterns, roi_selection, defaultParams, overwriteFlag) 
            % Make and save a .rpt or .hsrpt for each of the ROIs of each of the data files matching 
            % a patterns found in spdatapath. 
            % 
            % This will also auto track the files based on the default parameters.
            %
            % [IN]
            %  spdataPath - A path to a directory where .spdata files are to be tracked
            %  spdatafile_patterns - A cell-array of file patterns that can use wildcard '*'
            %                      to search for (a pattern can also just be a filename with no wildacard)
            %                      Ex: {'2015-01-01*.spdat', '2015-01-03-condition2.spdata'}
            %  roi_selection - Selection of 1 or more ROIs.  Roi selection is either 1 or more roi
            %                  indexes.  Or one or more ROInames.  If blank then select all ROI.
            %  defaultParams - An example RPT object or filename or preservedParams struct.  Or a cell
            %                  array of defaultParams matching each ROI by name.  If empty we open a
            %                  dialog to ask for an example file.
            %  overwriteFlag - [optional] Integer 0=Do not overwrite; [Default] 
            %                                     1=Warn with selection dialog before overwrite;
            %                                     2=Force overwrite (caution!);
            % [OUT]
            %  rpt_files - Cell array of full-paths to all .rpt files corresponding to the given 
            %             data file patterns. 
            %             (This list inculdes the default and non-overwitten files, as we assum you will
            %             want to process all these files similarly in the next phase of batch
            %             processing).
            if nargin==4
                overwriteFlag = 0;
            end
           
            spdatafile_patterns = makecell(spdatafile_patterns);
            roi_selection = makecell(roi_selection);
            defaultParams = makecell(defaultParams);
            if isempty(defaultParams) && numel(roi_selection)==1
                defaultParams = Pickle.selectExistingFileName(spdatapath,RPT.saveFileExt,RPT.SaveableDataFormats,...
                                        'Select .RPT with desired Default Properties');
            end
            if isempty(defaultParams)
                error('RPT:batchProcess','Unable to find default file or object');
            end
            spd_filenames = cellmap(@(p) Pickle.listExistingFileNames(spdatapath,p), spdatafile_patterns);
            spd_filenames = [spd_filenames{:}];
            try
                overwrite = (overwriteFlag==2) || RPT.confirmOverwriteDialog(); % Determine if we need to overwrite
            catch EX
                switch EX.identifier
                    case 'Pickle:OverwriteCancel'
                        rpt_files = {};
                        return;
                    otherwise
                        rethrow(EX);
                end
            end
                
            Ndatafiles = numel(spd_filenames);
            Ntracked = 0;
            Nerror = 0;
            Nexisting = 0;
            rpt_files = cell(1,Ndatafiles);
            H=waitbar(0,'Batch Processing ...');
            for n=1:Ndatafiles
                waitbar(n/(Ndatafiles+1),H,sprintf('Batch Processing Datafile [%i/%i] ...', n,Ndatafiles));
                spd = SPData(spd_filenames{n});
                spd.disableWaitbar=true;
                if isempty(roi_selection)
                    roi_idx = 1:numel(spd.ROI);
                else
                    [~,roi_idx] = cellmap(@(roi) spd.getROI(roi), roi_selection);
                end
                Nroi = numel(roi_idx);
                rpt_files{n} = cell(1,Nroi);
                for k=1:Nroi
                    fprintf('\n*** [%i/%i] Batch Tracking Dataset:%s ROIname:%s\n',n,Ndatafiles,spd.saveFileBaseName, spd.ROIname{roi_idx{k}}); 
                    file_path = spd.getROIFiles(roi_idx{k}, 'RPT');
                    if ~isempty(file_path) && ~overwrite
                        rpt_files{n}{k} = file_path{1}; % Record file
                        Nexisting = Nexisting + 1;
                        continue %Don't overwrite the existing files
                    end
                    try
                        rpt = spd.trackRPT(roi_idx{k}, defaultParams{max(k,numel(defaultParams))}, true);
                        rpt.disableWaitbar = true;
                        rpt.autoTrack();
                        rpt.save();
                    catch err
                        fprintf('>>>(oops)<<< Batch Tracking Error.\n');
                        disp(getReport(err))
                        rpt_files{n}{k} = [];
                        Nerror = Nerror+1;
                        continue;
                    end
                    Ntracked = Ntracked +1;
                    rpt_files{n}{k} = rpt.saveFilePath;
                end
            end
            if ~isempty(rpt_files) && isscalar(roi_selection) % Collapse down cell array of rptfiles if only one ROI processed per file
                rpt_files = cellmap(@(f) f{1}, rpt_files);
            end
            close(H);
            Ntotal = Ntracked+Nexisting+Nerror;
            fprintf('\n*** Batch Track Complete. [Num Total: %i | Num Tracked: %i | Num Existing: %i | Num Error: %i]\n',...
                        Ntotal,Ntracked,Nexisting,Nerror);
        end
    end %public static methods


    methods (Access = protected)
        function initializeBoxxer(obj, force)
            obj.checkPhase(2); 
            P = obj.checkParamsFindMaxima();
            if nargin==1
                force = false;
            end
            if force || isempty(obj.boxxer)
                sigmas = repmat( P.filterSigmas.*obj.data.psf, 2,1);
                obj.boxxer = Boxxer2D(obj.frameSize, sigmas);
                %clear filtered frames
                obj.filtered_frames_=[];
                obj.filtered_frames_loaded=false;
            end
        end
        
        function initializeEmitterModel(obj)
            % make sure the emitterModel transient property is initialized
            obj.checkPhase(4);
            R = obj.ResultsFilterMaxima;
            P = obj.ParamsLocalizeEmitters;
            emitterConstructor = str2func(P.model);
            size_cats = R.boxCoords.sizeCategories;
            Nsize_cats = R.boxCoords.NsizeCategories;
            obj.emitterModel = cellmap(@(i) emitterConstructor(size_cats(:,i), obj.data.psf), 1:Nsize_cats);
        end
        
        function loadData(obj, data, roi_in, roi_name)
            % roi_name: string (optional). Force overwrite of obj.ROIname to roi_name
            if isempty(data.psf)
                error('RPT:loadSPData','SPData.psf is not set');
            end
            if isempty(data.pixelSize)
                error('RPT:loadSPData','SPData.pixelSize is not set');
            end
            [obj.frames_, obj.ROI] = data.getFrames(roi_in);
            if isempty(obj.frames_)
                error('RPT:loadSPData','unable to load frames');
            end

            obj.frames_loaded = true;
            obj.data = data;
            obj.workingDir = obj.data.getFilePath('RPT');
            obj.Paths.data = relativepath(obj.workingDir, obj.data.saveFilePath);
            if isscalar(roi_in) && roi_in>0 && roi_in<=length(obj.data.ROI)  %Index into predefined ROI             
                obj.ROIname = obj.data.ROIname{roi_in};
                [~, file_name, file_ext] = obj.data.getROIFileNameParts(roi_in, 'RPT');
                obj.Paths.saveFile = [file_name, file_ext];
            elseif isempty(roi_in) %Full frame ROI
                obj.ROIname = 'FullFrameROI';
                obj.Paths.saveFile = sprintf('%s_%s%s',obj.data.saveFileBaseName,obj.ROIname,obj.saveFileExt);               
            else %Manual ROI
                obj.ROIname = 'ManualROI';
                pattern = sprintf('%s_%s%%i%s',obj.data.saveFileBaseName,obj.ROIname,obj.saveFileExt);
                obj.Paths.saveFile = Pickle.findUnusedFileName(obj.workingDir,pattern);
            end
            obj.sumImage = sumImage2D(data.getFrames(obj.ROI));
            if nargin == 4
                obj.ROIname = roi_name;
            end
            obj.initialized = true;
            obj.setPhase(2); %Initialized
        end
        
        function filterFrames(obj)
            % Actually do the filtering of the frames and set the properties
            obj.checkPhase(2); 
            P=obj.checkParamsFindMaxima();

            obj.initializeBoxxer();
            if obj.nScales == 1
                switch P.method
                    case 'LoG'
                        obj.filtered_frames_ = obj.boxxer.filterLoG(obj.getFrames());
                    case 'DoG'
                        obj.filtered_frames_ = obj.boxxer.filterDoG(obj.getFrames());
                    otherwise
                        error('RPT:filterFrames','Unknown filter method "%s"',P.method);
                end
            else
                switch P.method
                    case 'LoG'
                        obj.filtered_frames_ = obj.boxxer.filterScaledLoG(obj.getFrames());
                    case 'DoG'
                        obj.filtered_frames_ = obj.boxxer.filterScaledDoG(obj.getFrames());
                    otherwise
                        error('RPT:filterFrames','Unknown filter method "%s"',P.method);
                end
                % The C++ library has [x y t s] we want it as [x y s t]
                obj.filtered_frames_ = permute(obj.filtered_frames_,[1,2,4,3]);
            end
            if ~isempty(obj.filtered_frames_)
                obj.filtered_frames_loaded=true;
            end
        end

        %Plotting helper functions
        function im = computeMaximaImage(obj, maxima, maximaVals)
            mean_im = zeros(obj.frameSize);
            max_im = zeros(obj.frameSize);
            maximaVals = cosmicNorm(maximaVals);
            for mi=1:length(maximaVals)
                yi = maxima(1, mi);
                xi = maxima(2, mi);
                %si = maxima(3, mi); %scale index
                mean_im(yi,xi) = mean_im(yi,xi)+maximaVals(mi);
                max_im(yi,xi) = max(max_im(yi,xi),maximaVals(mi));
            end
            max_im = cosmicNorm(max_im);
            mean_im = mean_im/obj.nFrames;
            mean_im = mean_im./max(mean_im(:));
            im = 0.5*(mean_im+max_im);
        end

        function im = makeEmitterSuperResGauss(obj, E, imSizePx)
            % [in] E - emitters in internal emitters format
            physical_roi = obj.ROI(1:4);
            physical_roi(1) = physical_roi(1)-1;
            physical_roi(3) = physical_roi(3)-1;
            srr = SRRender2D(physical_roi,'single');
            points = E(:,[3,1,2,6,7,13]); %I x y sigmaX sigmaY frameIdx
            im = srr.renderGauss(points, imSizePx);
        end


        %% Abstract methods inherited from Pickle
        function val = getProtectedProperty(obj, name)
            %This is necessary for Pickle functionality to be able to access subclass protected variables
            val = obj.(name);
        end

        function modifyProtectedProperty(obj, name, newval)
            %This is necessary for Pickle functionality to be able to change subclass protected variables
            obj.(name)=newval;
        end
    end % Protected methods


        
end

