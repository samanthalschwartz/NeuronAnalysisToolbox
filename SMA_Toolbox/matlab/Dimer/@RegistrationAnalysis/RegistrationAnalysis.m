% RegistrationAnalysis - A class for loading/calibrating/localizing and
%  mapping channel registration data sets.
%
%  Mark J. Olah (mjo@cs.unm.edu)
%  2015
%
%

classdef RegistrationAnalysis < Pickle & GUIBuilder
    % RegistrationAnalysis - A class for loading/calibrating/localizing and
    % mapping channel registration data sets.  This class can load data 
    % aquired by many different systems including just sequences of static
    % beads, or other arbitrary point emitters.
    %
    % Images and Coordinate systems: 
    %  * Frist Read this: https://www.mathworks.com/help/images/image-coordinate-systems.html
    %  * The upper left corner of the image is taken to be (0,0).
    %  * This is similar to the matlab "intrinsic" image coorindate system, with the exception
    %    that our origin is at (0,0) while matlab uses (0.5, 0.5) which seems less natural to us.
    %  * Units are always in absoulte coordinates from the upper left
    %    of the image with the y-axis going down and x-axis to the right.
    %  * Most internal data structures work in pixels as the unit, although
    %    we provide mappings in either pixels or microns, where the pixelSize property determines the conversion to microns.
    %  

    properties        
        % General User specified properties
        frameSize; %[#rows #col] i.e., [sizeY sizeX] - The frame size in pixels.  All frames from all must be of this size
        pixelSize; %[sizeX, sizeY] - The pixel size in microns.
        CCDBackground; % The EMCCD camera background constant [ADU]
        CCDGain; % The EMCCD camera gain constant in [e-/ADU]

        %Advanced parameters controlling the localization and filtering
        Params = struct(...
            'BoxSize', [9,9], ...     % Size of boxes in pixels
            'MinBoxes', 0, ... % Attempt to detect at least this many boxes for every frame
            'NeighborhoodSize',int32(7),...   % Size of filtering neighborhood >=3 and odd.
            'PSFSigma', [1,1], ...    % This should remain [1;1] because of the localization
            'EstimationMethod','Newton',...  % The method to use in localization.
            'BoundaryWidth',5, ... %Distance away from edge that is considered within the ROI in pixels.
            'PixelSaturationThreshold', 2^14-2^10,... %Pixel value from raw image above which saturation is met.
            'MinIntensity',400,...  % Minimum intensity of fit emitter
            'MinSigma',0.8,...       % Minimum sigma of fit emitter (px)
            'MaxSigma',3,...         % Maximum sigma of fit emitter (px)
            'MaxPositionSE',0.2,... % In max x,y SE in pixels            
            'MinLocalizationDist',4, ... % Minimum overlap distance for nearby localizations
            'MaxRefDisplacement',4 ... % Maximum allowed displacement for reference point pairing
            );

        data = struct(... % make a 0x0 structure array with fixed fields
                'sourceName',{},... % An identifiing name for the data source (normally filename, but sometimes a filedname also)
                'dateTime',{},... % A DateTime for aquisition (if availible)
                'Nframes',{},... %Integer count of frames
                'rawSequence',{},...  % The raw sequence of data
                'saturatedPixels',{},... % Array Nx3 cols:[X,Y,T] the locationsr of saturated pixels detected
                'NSaturatedPixels',{},... % The total number of saturated pixels detected
                'valid',{},... %True if any valid frames
                'validFrames',{},...  % boolean:[Nframes,1] true=keep frame false=discard frame
                'sequence',{},...     % The gain calibrated sequence of data
                'rawBoxCoords',{},...    % A BoxCoords object recording our box locations
                'rawLocalizations',{},... % A Nx11 array of the rawLocalizations in absolute frame (pixel) coordinates cols:[x y I bg sigma SEx SEy SEI SEbg SEsigma FrameIdx]
                'localizationFilter',{},... % 0=keep non-zero indicates reason for removal
                'NLocalizations',{},... % A count of valid localizations from both channels
                'localizations',{},... % A Nx11 array of the filtered localizations in absolute frame (pixel) coordinates cols:[x y I bg sigma SEx SEy SEI SEbg SEsigma FrameIdx]
                'boxCoords',{},... % A boxcoords object for the filtered localizations
                'NChEmitters',{},... % Nframes x2 matrix of number of final emitters founf in left and right channels for each frame
                'NRefPts',{},... % A Nframes x 1 vector counting valid reference points
                'refPts',{},...    % A 1x2 cell array of Nx2 refrence points positions cols:[x, y];
                'refPtsSE',{},...    % A 1x2 cell array of Nx2 refrence point standard errors cols:[SEx, SEy];
                'refPtsIdx',{},...  % An Nx2 array of refrence points indexes (rowIdxs) into the localizations array cols:[idxL, idxR];
                'RMSDisp',{},... % (px) Root mean squared displacement of the reference points with simple shift mapping
                'MaxDisp',{}... % (px) Max displacement of the reference points with simple shift mapping
                );

        maps = struct(...; % structure of named channel mappings along with relevent parameters links to training and test data
                'name',{},...      %  name - A name for this mapping
                'algorithm',{},... %  algorithsm - The name of the algorithm used
                'valid',{}, ...    % true if a valid map
                'params',{},...    %  params - A structure of parameters used for the given algorithm
                'trainRefPts',{},...  % The refererence point pairs extracted from the train data: 1x2 cell array of Nx2 refrence points cols:[x, y];
                'trainDataIdxs',{},...% A Nx2 array of indexs cols:[seqIdx, frameIdx]
                'testRefPts',{},...   % The refererence point pairs extracted from the test data: 1x2 cell array of Nx2 refrence points cols:[x, y];
                'testDataIdxs',{},... % A Nx2 array of indexs cols:[seqIdx, frameIdx]
                'NTrainFrames',{},... % The number of training reference frames used;
                'NTestFrames',{},...  % The number of testing reference frames used;
                'NTrainPts',{},...    % The number of training reference points from all train data frames;
                'NTestPts',{},...     % The number of testing reference points from all test data frames;
                'RMSE',{},...      %  RMSE - The root-mean-squared-error for the mapping as measured with the test data
                'maxError',{},...  %  maxError - The maximum displacement error observed over the test data set.           
                'mapFunctionPixels',{},...%  mapFunctionPixels - A function mapping from translates ch2->ch1 points arrays [Units: pixels].
                'mapFunctionMicrons',{}...%  mapFunctionMicrons - A function mapping from translates ch2->ch1 points arrays [Units:microns].
               );
    end

    properties
        Paths = struct('saveFile',[]); %All filenames are relative paths from the workingDir
    end

    properties (Hidden=true)
        version = 1; %For future file format version changes
    end

    properties (Dependent=true)
        ChannelROI; %[2x1] cell array of [xmin, xmax, ymin, ymax] bounding boxes
        NData;
        NMaps;
        SensorXSplit; %SensorSizeX/2
        SensorSizeX; 
        SensorSizeY;
        calibrated; % Boolean if the object is calibrated which means all of its parameters have been set and adding data leads to localization and mapping
        hasValidData; % Boolean: True if there is at least some valid data (data with successful reference localization)
        hasValidMap; % Boolean: True if there is some valid mapping (even the Null mapping)
    end
    
    properties (Access=protected, Transient=true)
        gainCalGuiFig;
        figHs; % cell array of figure handles created by the dip image gain calibration
    end    

    properties (Constant=true, Hidden=true)
        %Abstract properties inherited from Pickle
        saveFileExt = '.reganalysis'
        SaveableDataFormats = {'*.reganalysis', 'RegistrationAnalysis (.reganalysis)'};
        LoadableDataFormats = {'*.reganalysis;*.mat', 'All Loadable Formats (.reganalysis;.mat)';...
                               '*.reganalysis', 'RegistrationAnalysis (.reganalysis)';...
                               '*.mat', 'Mat-Files (.mat)'};
        LoadableRawDataFormats = {'*.mat', 'Mat-Files (.mat)'};
        
        MapAlgorithms = {'Null',... % Direct channel overlay with no use of control points (i.e. the baseline)
                         'GlobalAffine',... % Single global affine transform (uses fitgeotrans)
                         'LocalAffine',... % Locally weighted affine transform.
                         'SmoothAffine',... % Smoothed locally weighted affine transform (faster than local)
                         'LWM',... % Locally weighted mean (uses fitgeotrans)
                         'Polynomial',... % Ploynomial interpolation (uses fitgeotrans)
                         'NRS'... % Non-reflective similarity (uses fitgeotrans)
                        };
        MapAlgorithmDefaultParams = struct(...
                        'Null',[],...
                        'GlobalAffine',[],...
                        'LocalAffine',struct('NFitPoints',9,'NWeightPoints',4),...
                        'SmoothAffine',struct('NFitPoints',9,'WeightExponent',1),...
                        'LWM',struct('NFitPoints',12),...
                        'Polynomial',struct('Degree',3),...
                        'NRS',[]...
                        );
        NChannels=2; %Only 2 channels handled for now until we have a mircoscope that can do 4-channels...
        LocalizationFilterDescriptions = {'Boundary','Saturation','MinIntensity','MinSigma','MaxSigma','MaxPosSE','Distance'};
        DefaultPixelSaturationThreshold = 2^14-2^10; %Pixel value from raw image above which saturation is met.
        MinTrainPoints = 12; %Minimum number of training points needed to adequately train
    end



    methods
        function obj = RegistrationAnalysis(varargin)
            % RegistrationAnalysis
            %    Arguments are the same as for a the call to obj.initialize()
            %           
            if nargin>0
                obj.load(varargin{:});
            end
        end
        
        function load(obj, filename)
            %
            % Option 1:
            %  initialize - with RegistrationAnalysis object or saved registration analysis .reg file
            %
            %
            if iscell(filename) && numel(filename)==1
                filename = filename{1}; % Treat singleton cell-arrays the same as a single filename.
            end
            if ischar(filename)
                [~,base,ext] = fileparts(filename);
                obj.updateWaitbar(0,sprintf('Loading %s ...',[base ext]));
                switch ext
                    case obj.saveFileExt
                        obj.loadRA(filename); % Load from RegistrationAnalysis object or saved .reg file
                    case '.mat'
                        obj.loadMatData(filename);
                    otherwise
                        error('RegistrationAnalysis:load','Unknown file extension: %s',ext);
                end
            elseif iscell(filename)
                filelist = filename;
                if ~all(cellfun(@ischar, filelist))
                    error('RegistrationAnalysis:load','Loading of cell arrays only supported for cell arrays of filenames');
                end
                if ~all(cellfun(@(f) exist(f,'file'), filelist));
                    error('RegistrationAnalysis:load','Not all filenames exist');
                end
                if any(cellfun(@(f) ~isempty(strfind(f,obj.saveFileExt)),filelist))
                    error('RegistrationAnalysis:load','Cannot batch load .reganalysis files through load().');
                end
                N = numel(filelist);
                for n = 1:N
                    [~, base, ext] = fileparts(filelist{n});
                    obj.updateWaitbar((n/N)*.9,sprintf('Processing %s ...',[base ext]));
                    obj.loadMatData(filelist{n});
                end
            else
                error('RegistrationAnalysis:load','Bad arguments');
            end
            obj.reprocessMaps();
        end

        function loadMatData(obj,matFile)
            % Loads data from any saved file in matlab .mat format (regardless of actual file extension).
            % 
            % This is the only way to add new data.
            %
            % [in]
            %   matFile - String: a filename to open and load registration images sets from.
            S = load(matFile,'-mat');
            if isfield(S,'sequence')
                obj.loadSequence(matFile); %Load from generic sequence
            elseif isfield(S,'trn') && isfield(S,'tst')
                obj.loadHMMFiducial(matFile); % Load from HMM fiducial images
            elseif isfield(S,'Params')
                if isfield(S.Params,'CameraObj')
                    obj.loadSPTRegisterChannels(matFile); %This file was made by SPTRegisterChannels
                else
                    obj.loadLidkeFormat2(matFile); %This file is Lidke Format 2 (Circa 2016)
                end
            elseif isfield(S,'CRobj')
                obj.loadChannelRegistrationV1(matFile); % This is a ChannelRegistrationV1 save file
            else
                error('RegistrationAnalysis:load','Unknown .mat file format: %s',matFile);
            end
        end

        function [ccd_gain, ccd_background] = calibrateCCD(obj, beadsIm, bgIm)
            % Get a gain/offset calibration for a CCD.  Use dipimage calibration if given images.
            % Otherwise scalars can be given if already known and they will be check for validity.
            % [in]
            %   [Form1]: Give the calibration image sequences
            %       beadsIm - The gain calibration beads sequence, or a filename of a .mat containing such a
            %                 variable 'sequence' or 'Data'.
            %       bgIm - The gain calibration background sequence, or a filename of a .mat containing such a
            %                 variable 'sequence' or 'Data'.
            %   [Form2]: Give the CCD parameters
            %       ccd_gain - The gain in [e-/ADU]
            %       ccd_background - The background offset in [ADU]
            % [out]
            %   ccd_gain - CCD gain [e-/ADU]
            %   ccd_background - CCD background offset [ADU]
            if nargin<3; 
                bgIm = obj.CCDBackground;
                beadsIm = obj.CCDGain;
            end
        
            if isscalar(beadsIm) && isscalar(bgIm) %Check if user has given gain/offset values
                ccd_gain = beadsIm;
                ccd_background = bgIm;
                if ccd_background<0 || ccd_gain<=0
                    error('RegistrationAnalysis:calibrateCCD','Bad CCD Gain parameters: %.3g gain %.3g offset.',ccd_gain,ccd_background);
                end
            else
                if ischar(beadsIm)
                    Beads = load(beadsIm,'-mat');
                else
                    Beads.sequence = beadsIm;
                end
                if ischar(bgIm)
                    Bg = load(bgIm,'-mat');
                else
                    Bg.sequence = bgIm;
                end
                if isfield(Beads,'sequence')
                    beadsseq = Beads.sequence;
                elseif isfield(Beads,'Data')
                    beadsseq = Beads.Data;
                end

                if isfield(Bg,'sequence')
                    bgseq = Bg.sequence;
                elseif isfield(Bg,'Data')
                    bgseq = Bg.Data;
                end
                %Do cal_readnoise to get the calibrations and get rid of annoying figs it auto-pops-up
                fs = findobj('Type','figure','NumberTitle','on'); %record current figs
                cal_out = cal_readnoise(beadsseq, bgseq);
                obj.figHs = setdiff(findobj('Type','figure','NumberTitle','on'),fs); %Look for new figs

                ccd_background = cal_out(4); % background [ADU]
                ccd_gain = cal_out(2); % gain [e-/ADU]
            end
            obj.CCDGain = ccd_gain;
            obj.CCDBackground = ccd_background;
            obj.reprocessData();
        end

        function D = makeData(obj, seq, sourceName, dateTime)
            % Make a new structure for a squence data.
            if isa(seq,'dip_image')
                seq = single(seq);
            end
            sz = size(seq);
            if ~isempty(obj.frameSize) && ~all(obj.frameSize==sz(1:2))
                error('RegistrationAnalysis:BadFrameSize','Frame Size does not match');
            end
            if nargin>2
                D.sourceName = sourceName;
            else
                D.sourceName = [];
            end
            if nargin>3
                D.dateTime = dateTime;
            else
                D.dateTime = [];
            end
            D.Nframes = size(seq,3);
            D.rawSequence = seq; % convert from dipimage
            [Ridx,Cidx,Fidx] = ind2sub(sz,find(seq(:)>=obj.Params.PixelSaturationThreshold));
            D.saturatedPixels = [Ridx,Cidx,Fidx];
            D.NSaturatedPixels = arrayfun(@(n) sum(Fidx==n), (1:D.Nframes)');
            D.validFrames = true(D.Nframes,1);
%             if ~isempty(D.saturatedPixels)                
%                 D.validFrames(unique(D.saturatedPixels(:,3))) = false; % invalidate any frames with saturated pixels
%             end
            D.valid = any(D.validFrames);
            D.sequence = [];
            D.rawBoxCoords = [];
            D.rawLocalizations = [];
            D.localizationFilter = [];
            D.NLocalizations = 0;
            D.localizations = [];
            D.boxCoords = [];
            D.NChEmitters = zeros(D.Nframes,2); 
            D.NRefPts = zeros(D.Nframes,1);
            D.refPts = [];
            D.refPtsSE = [];
            D.refPtsIdx = [];
            D.RMSDisp = zeros(D.Nframes,1);
            D.MaxDisp = zeros(D.Nframes,1);
            D = obj.processData(D);
        end

        function D = processData(obj, D)
            % Actually process a loaded valid and calibrated data set to determine what (if any) usefull 
            % data exists.  Uses the obj.Params struct to control most of the processing parameters
            if ~obj.calibrated
                return;
            end
            D.sequence = obj.calibrateImage(D.rawSequence);
            if ~any(D.validFrames)
                return;
            end

            %Find boxes
            bxr = Boxxer2D(obj.frameSize, obj.Params.PSFSigma);
            [mx, mx_vals] = bxr.enumerateImageMaxima(bxr.filterDoG(D.sequence), obj.Params.NeighborhoodSize);
            fmx = bxr.filterMaximaAuto(mx, mx_vals, obj.Params.MinBoxes);
            D.rawBoxCoords = bxr.generateBoxCoords(fmx, D.Nframes, obj.Params.BoxSize);
            %Localize
            emitterImages = D.rawBoxCoords.makeROI(D.sequence); % make the emitter images
            g2d = Gauss2DsMAP(obj.Params.BoxSize, obj.Params.PSFSigma);
            [etheta, crlb, llh] = g2d.estimate(emitterImages, obj.Params.EstimationMethod);
            D.rawLocalizations = double([D.rawBoxCoords.covertToEmitterAbsoluteLocations(etheta); sqrt(crlb); llh'; single(D.rawBoxCoords.boxFrameIdx)]');
            %Filter localizations
            E = D.rawLocalizations; % E - The emitter matrix used for filtering columns are [x y I bg sigma, SEx, SEy, SEI, SEbg, SEsigma, FrameIdx] units are pixels.
            N = size(E,1);
            P = obj.Params;
            
            filter = zeros(N,1); %true if we will keep the localizations
            %Filter - Boundary (1)
            filterStepCode = 1; % increment this after each test.  This indicated which filter method was repsonsible for filtering point
            if ~isempty(P.BoundaryWidth) && P.BoundaryWidth>0
                b = P.BoundaryWidth;
                x=E(:,1);
                y=E(:,2);
                valid = x>b & abs(x-obj.SensorXSplit)>b & obj.frameSize(2)-x>b & ...
                        y>b & obj.frameSize(1)-y>b;
                filter(~filter & ~valid) = filterStepCode;
            end
            %Filter - Saturation (2)
            filterStepCode = filterStepCode+1;
            if ~isempty(P.PixelSaturationThreshold) && P.PixelSaturationThreshold>0
                rawEmitterImages = D.rawBoxCoords.makeROI(D.rawSequence);
                maxI = squeeze(max(max(rawEmitterImages,[],1),[],2));
                saturated = maxI >= P.PixelSaturationThreshold;
                filter(~filter & saturated) = filterStepCode;
            end
            %Filter - MinIntensity (3)
            filterStepCode = filterStepCode+1;
            if ~isempty(P.MinIntensity) && P.MinIntensity>0
                filter(~filter & E(:,3)<P.MinIntensity) = filterStepCode;
            end
            %Filter - MinSigma (4)
            filterStepCode = filterStepCode+1;
            if ~isempty(P.MinSigma) && P.MinSigma > 0
                filter(~filter & E(:,5)<P.MinSigma) = filterStepCode;
            end
            %Filter - MaxSigma (5)
            filterStepCode = filterStepCode+1;
            if ~isempty(P.MaxSigma) && P.MaxSigma > 0
                filter(~filter & E(:,5)>P.MaxSigma) = filterStepCode;
            end
            %Filter - MaxPositionSE (6)
            filterStepCode = filterStepCode+1;
            if ~isempty(P.MaxPositionSE) && P.MaxPositionSE > 0
                filter(~filter & (E(:,6)>P.MaxPositionSE | E(:,7)>P.MaxPositionSE)) = filterStepCode;
            end
            %Filter - MinLocalizationDist (7)
            filterStepCode = filterStepCode+1;
            if ~isempty(P.MinLocalizationDist) && P.MinLocalizationDist > 0              
                for frame = 1:D.Nframes
                    fidx = D.rawBoxCoords.frameIndexes{frame};
                    fidx = fidx(~filter(fidx)); % Select only non-filtered points.
                    frameE = E(fidx,:);
                    nFrameE = numel(fidx);
                    if nFrameE > 1
                        for n = 1:nFrameE-1
                            dists = sqrt(sum((repmat(frameE(n,[1,2]),nFrameE-n,1)-frameE(n+1:end,[1,2])).^2,2));                            
                            nearby = [n; n+find(dists <= P.MinLocalizationDist)];
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
            D.localizationFilter = filter;
            keepers = ~filter; % Keep anything that did not have a "code" for rejections
            D.localizations = D.rawLocalizations(keepers,:);
            D.NLocalizations = size(D.localizations,1);
            D.boxCoords = BoxCoords.filter(D.rawBoxCoords, keepers);
            %Identify reference points
            E = D.localizations;
            leftMask = E(:,1) < obj.SensorXSplit;
            approxShift = [obj.SensorXSplit, 0]; % approximate shift of ch2 relative to ch1.  This is subtracted from ch2 points to create the simplest mapping for identifiying the refPts
            D.NRefPts = zeros(D.Nframes,1);
            D.refPts = cell(D.Nframes, 1);
            D.refPtsSE = cell(D.Nframes, 1);
            D.refPtsIdx = cell(D.Nframes, 1);
            D.NRefPts = zeros(D.Nframes, 1);
            D.NChEmitters = zeros(D.Nframes, 2);
            D.RMSDisp = zeros(D.Nframes,1);
            D.MaxDisp = zeros(D.Nframes,1);
            for n=1:D.Nframes
                frameMask = E(:,end)==n;
                chIdxs = {find(leftMask & frameMask),find((~leftMask) & frameMask)}; %indexes of left side and right side emitters
                NChEmitters = cellfun(@numel, chIdxs);
                D.NChEmitters(n,:) = NChEmitters;
                if any(NChEmitters<1)
                    continue;
                end
                lPos = E(leftMask&frameMask,1:2);
                rPos = E((~leftMask)&frameMask,1:2) - repmat(approxShift,numel(chIdxs{2}),1); %Shift back by approxShift
                tempRefPtsIdx = zeros(NChEmitters(1), 1);
                displacements = zeros(NChEmitters(1), 1);
                for k=1:NChEmitters(1)
                    sqdist = sum((repmat(lPos(k,:),NChEmitters(2),1)-rPos).^2,2);
                    [~,minidx] = min(sqdist);
                    mindist = sqrt(sqdist(minidx));
                    if mindist <= P.MaxRefDisplacement
                        other = find(tempRefPtsIdx(1:k-1) == minidx,1,'first');
                        if ~isempty(other) 
                            if mindist>displacements(other)
                                continue;
                            else % Steal this point from other ch1 point. we are closer
                                displacements(other)=0;
                                tempRefPtsIdx(other)=0;
                            end
                        end
                        tempRefPtsIdx(k) = minidx;
                        displacements(k) = mindist;
                    end
                end
                if ~any(tempRefPtsIdx)
                    % Did not find any pairings
                    D.validFrames(n) = false;
                    continue;
                end                
                D.refPtsIdx{n} = [chIdxs{1}(tempRefPtsIdx>0), chIdxs{2}(tempRefPtsIdx(tempRefPtsIdx>0))];
                D.refPtsSE{n} = {E(D.refPtsIdx{n}(:,1),6:7),  E(D.refPtsIdx{n}(:,2),6:7)};
                D.refPts{n} = {E(D.refPtsIdx{n}(:,1),1:2),  E(D.refPtsIdx{n}(:,2),1:2)};
                D.NRefPts(n) = size(D.refPts{n}{1},1);
                assert(numel(unique(D.refPtsIdx{n}(:,2)))==size(D.refPtsIdx{n},1));
                displacements = displacements(tempRefPtsIdx>0);
                D.RMSDisp(n) = sqrt(mean(displacements.^2));
                D.MaxDisp(n) = max(displacements);
            end
            D.valid = any(D.validFrames);
        end

        function appendData(obj, D)
            if isempty(obj.data)
                obj.data=D;
            else
                obj.data(end+1,1)=D;
            end
        end
        
        function clearData(obj)
            %Clear all data from object
            obj.data(:) = []; 
            obj.frameSize = [];
            obj.initialized = false;
            obj.reprocessMaps();
        end

        function deleteData(obj, frameIdxs)
            obj.data(frameIdxs) = []; 
            if obj.NData<1
               obj.frameSize = [];
               obj.initialized = false;
            end
            obj.reprocessMaps();
        end

        function invalidateFrame(obj, idx)
            % idx - a 2x1 array of [seqIdx, frameIdx].
            obj.data(idx(1)).validFrames(idx(2)) = false;
            obj.reprocessMaps();
        end

        function data = getDataCellArray(obj)
            if obj.NData<1
                data = [];
                return;
            end
            seqIdx =   cellcat(cellmap(@(n) repmat({n},obj.data(n).Nframes,1), (1:obj.NData)'));
            frameIdx = num2cell(cellcat(cellmap(@(n) (1:n)', {obj.data.Nframes}')));
            valid =   num2cell(vertcat(obj.data.validFrames));
            dateTime = cellcat(cellmap(@(dt,n) repmat({datestr(dt,'[YYYY-MM-DD] HH:mm:ss')},n,1), {obj.data.dateTime}', {obj.data.Nframes}'));

            NChEmitters = num2cell(vertcat(obj.data.NChEmitters));
            NRefPts = num2cell(vertcat(obj.data.NRefPts));
            NSaturated = num2cell(vertcat(obj.data.NSaturatedPixels));
            RMSD = num2cell(vertcat(obj.data.RMSDisp));
            MaxDisp = num2cell(vertcat(obj.data.MaxDisp));
            SourceName = cellcat(cellmap(@(dt,n) repmat({dt},n,1), {obj.data.sourceName}', {obj.data.Nframes}'));
            data = [seqIdx, frameIdx, dateTime,valid, NRefPts, NChEmitters(:,1), NChEmitters(:,2),...
                    NSaturated,RMSD,MaxDisp,SourceName]; 
        end

        function M = computeMap(obj, algorithm, params, trainDataIdxs, testDataIdxs)
            allAlgs = fieldnames(obj.MapAlgorithmDefaultParams);
            algorithm = allAlgs{find(strcmpi(algorithm, allAlgs),1,'first')};
            M.name = algorithm;
            M.algorithm = algorithm;
            if nargin<3 || isempty(params) 
                M.params = obj.MapAlgorithmDefaultParams.(algorithm);
            else
                M.params  = params;
            end
            if nargin<4
                [M.trainDataIdxs, M.testDataIdxs] = obj.getTrainTestDataIdxs();
            else
                M.trainDataIdxs = trainDataIdxs;
                M.testDataIdxs = testDataIdxs;
            end
            
            
            M.trainRefPts = obj.collectRefPts(M.trainDataIdxs);
            M.testRefPts = obj.collectRefPts(M.testDataIdxs);
            M.NTrainFrames = size(M.trainDataIdxs,1);
            M.NTestFrames = size(M.testDataIdxs,1);
            M.NTestPts = size(M.testRefPts{1},1);
            M.NTrainPts = size(M.trainRefPts{1},1);
            M.mapFunctionPixels = [];
            M.mapFunctionMicrons = [];
            M.RMSE = []; 
            M.maxError = []; 
            M.valid = false;
            if M.NTrainPts>=obj.MinTrainPoints
                trainRefPtsMicrons = cellmap(@(P) P.*repmat(obj.pixelSize,size(P,1),1), M.trainRefPts);
                try
                    switch algorithm
                        case 'Null'
                            M.mapFunctionPixels = obj.generateChannelMap_null([obj.SensorXSplit,0]);
                            M.mapFunctionMicrons = obj.generateChannelMap_null([obj.SensorXSplit.*obj.pixelSize,0]);
                        case 'GlobalAffine'
                            M.mapFunctionPixels = obj.generateChannelMap_fitgeotrans(M.trainRefPts, 'affine');
                            M.mapFunctionMicrons = obj.generateChannelMap_fitgeotrans(trainRefPtsMicrons, 'affine');
                        case 'LocalAffine'
                            M.mapFunctionPixels = obj.generateChannelMap_localAffine(M.trainRefPts, M.params.NFitPoints, M.params.NWeightPoints);
                            M.mapFunctionMicrons = obj.generateChannelMap_localAffine(trainRefPtsMicrons, M.params.NFitPoints, M.params.NWeightPoints);
                        case 'SmoothAffine'
                            M.mapFunctionPixels = obj.generateChannelMap_smoothAffine(M.trainRefPts, M.params.NFitPoints, M.params.WeightExponent);
                            M.mapFunctionMicrons = obj.generateChannelMap_smoothAffine(trainRefPtsMicrons, M.params.NFitPoints, M.params.WeightExponent);
                        case 'LWM'
                            M.mapFunctionPixels = obj.generateChannelMap_fitgeotrans(M.trainRefPts, 'lwm', M.params.NFitPoints);
                            M.mapFunctionMicrons = obj.generateChannelMap_fitgeotrans(trainRefPtsMicrons, 'lwm', M.params.NFitPoints);
                        case 'Polynomial'
                            M.mapFunctionPixels = obj.generateChannelMap_fitgeotrans(M.trainRefPts, 'polynomial', M.params.Degree);
                            M.mapFunctionMicrons = obj.generateChannelMap_fitgeotrans(trainRefPtsMicrons, 'polynomial', M.params.Degree);
                        case 'NRS'
                            M.mapFunctionPixels = obj.generateChannelMap_fitgeotrans(M.trainRefPts, 'nonreflectivesimilarity');
                            M.mapFunctionMicrons = obj.generateChannelMap_fitgeotrans(trainRefPtsMicrons, 'nonreflectivesimilarity');
                    end
                catch Err
                    switch Err.identifier
                        case 'images:geotrans:requiredNonCollinearPoints'
                            fprintf('Failure to run fitgeotrans for: %s\n%s\n%s\n',algorithm,Err.identifier, Err.message);
                            return;
                        otherwise
                            Err.rethrow
                    end
                end
                [M.RMSE, M.maxError] = obj.testChannelMap(M.mapFunctionPixels, M.testRefPts);
                M.valid = true;
            end
        end

        function appendMap(obj, M)
            if isempty(obj.maps)
                obj.maps=M;
            else
                obj.maps(end+1,1)=M;
            end
        end
        
        function clearMaps(obj)
            %Clear all maps from object
            obj.maps(:) = []; 
        end

        function deleteMap(obj, idxs)
            obj.maps(idxs) = []; 
        end
        
        function data = getMapsCellArray(obj)
            if obj.NMaps<1
                data = [];
                return;
            end
            names = {obj.maps.name};
            algs = {obj.maps.algorithm};
            valid = {obj.maps.valid};
            rmse = {obj.maps.RMSE};
            maxError = {obj.maps.maxError};
            NTrainPts = {obj.maps.NTrainPts};
            NTestPts = {obj.maps.NTestPts};
            NTrainFrames = {obj.maps.NTrainFrames};
            NTestFrames = {obj.maps.NTestFrames};
            trainDataIdxs = cellmap(@(d) strjoin(cellmap(@mat2str,num2cell(d,2)),'; '),{obj.maps.trainDataIdxs});
            testDataIdxs = cellmap(@(d) strjoin(cellmap(@mat2str,num2cell(d,2)),'; '),{obj.maps.testDataIdxs});
            data = vertcat(names, algs, valid, rmse, maxError, NTrainPts, NTestPts, NTrainFrames, NTestFrames, ...
                            trainDataIdxs, testDataIdxs)'; 
        end
        
        function [bestMapFunc, bestMapStruct] = getOptimalMapPixels(obj)
            % Get the optimal mapping in pixels
            valid = [obj.maps(:).valid];
            if ~any(valid)
                error('RegistrationAnalysis:getOptimalMapPixels','No valid mappings.  Check data sources.');
            end
            [~,mapIdx] = min([obj.maps.RMSE]);
            bestMapStruct = obj.maps(mapIdx);
            bestMapFunc = bestMapStruct.mapFunctionPixels;
        end
        
        function [bestMapFunc, bestMapStruct] = getOptimalMapMicrons(obj)
            [~,mapIdx] = min([obj.maps.RMSE]);
            bestMapStruct = obj.maps(mapIdx);
            bestMapFunc = bestMapStruct.mapFunctionMicrons;
        end
        
        function Ts = shiftRPTTracksMicrons(obj,rpt)
            % Use the optimal mapping to shift tracks from an rpt object and return them in microns.
            % We assume the tracked dataset is on the right side.  Returned tracks are in units of microns
            % with respect to the left side of the image
            % [in]
            %   rpt - A RPT object fully tracked for the right side in a region correpsonding to this
            %         RegistrationAnalysis object
            % [out]
            %   Ts - A cellarray of RPT tracks in units of microns.
            map = obj.getOptimalMapMicrons();
            Ts = rpt.getTracks();
            for n=1:numel(Ts)
                Ts{n}(:,[2,3]) = map(Ts{n}(:,[2,3]));
            end
        end
        
        function Ts = shiftRPTTracksPixels(obj,rpt)
            % Use the optimal mapping to shift tracks from an rpt object and return them in pixels.
            % We assume the tracked dataset is on the right side.  Returned tracks are in units of pixels
            % with respect to the left side of the image
            % [in]
            %   rpt - A RPT object fully tracked for the right side in a region correpsonding to this
            %         RegistrationAnalysis object
            % [out]
            %   Ts - A cellarray of RPT tracks in units of pixels.
            map = obj.getOptimalMapPixels();
            Ts = rpt.getTracksPixels();
            for n=1:numel(Ts)
                Ts{n}(:,[2,3]) = map(Ts{n}(:,[2,3]));
            end
        end

        function Ts = shiftTracksPixels(obj,Ts)
            % Use the optimal mapping to shift tracks given in units of pixels.
            % We assume the tracked dataset is on the right side.  Returned tracks are in units of pixels
            % with respect to the left side of the image
            % [in]
            %   Ts - A cell array of tracks (in the RPT format) in units of
            %   pixels
            % [out]
            %   Ts - A cellarray of RPT tracks in units of pixels.
            map = obj.getOptimalMapPixels();
            for n=1:numel(Ts)
                Ts{n}(:,[2,3]) = map(Ts{n}(:,[2,3]));
            end
        end
        
        %% Visualization and GUI 
        
        function f = viewSequence(obj, seqIdx)
            if ~obj.calibrated
                error('RegistrationAnalysis:calibration','Object is not calibrated.  Check pixelSize and gain/background parameters are set.');
            end
            f = obj.viewMaximizedDipFig(obj.data(seqIdx).sequence);
        end

        function f = viewRawSequence(obj, seqIdx)
            f = obj.viewMaximizedDipFig(obj.data(seqIdx).rawSequence);
        end

        function f = viewSequenceSaturation(obj, seqIdx)
            f = obj.viewMaximizedDipFig(obj.data(seqIdx).rawSequence);
            dipmapping(f,'saturation');
        end

        function f = viewSelectedBoxes(obj, seqIdx)
            box_im = obj.data(seqIdx).rawBoxCoords.plotBoxes(obj.data(seqIdx).sequence,obj.data(seqIdx).localizationFilter, ['Accepted', obj.LocalizationFilterDescriptions]);
            f = obj.viewMaximizedDipFig(box_im);
        end
        
        function plotEmitterIntensity(obj, axH, seqIdxs)
            if nargin<3 || isempty(seqIdxs)
                seqIdxs = 1:obj.NData;
            end
            seqIdxs( [obj.data.NLocalizations]<=0 ) = []; % Delete any sequences with no localizations
            if isempty(seqIdxs);
                return;
            end
            if nargin == 1 || isempty(axH)
                f=figure();
                whitebg(f);
                axH = axes();
            end
            axes(axH);
            Is = cellmap(@(i) obj.data(i).localizations(:,3), seqIdxs);
            names = cellmap(@(i) sprintf('[Seq:%i] %s',i,obj.data(i).sourceName),seqIdxs);
            BaseRPT.plot_distributions(axH,Is, names, true);
            xlabel('Intensity (photons)');
            title('Emitter Intensity Distribution');
        end

        function plotEmitterSigma(obj, axH, seqIdxs)
            if nargin == 1 || isempty(axH)
                f=figure();
                whitebg(f);
                axH = axes();
            end
            if nargin<3 || isempty(seqIdxs)
                seqIdxs = 1:obj.NData;
            end
            axes(axH);
            Sigmas = cellmap(@(i) obj.data(i).localizations(:,5), seqIdxs);
            names = cellmap(@(i) sprintf('[Seq:%i] %s',i,obj.data(i).sourceName),seqIdxs);
            BaseRPT.plot_distributions(axH,Sigmas, names, true);
            xlabel('Sigma (pixels)');
            title('Emitter Sigma Distribution');
        end

        function plotEmitterLLH(obj, axH, seqIdxs)
            if nargin == 1 || isempty(axH)
                f=figure();
                whitebg(f);
                axH = axes();
            end
            if nargin<3 || isempty(seqIdxs)
                seqIdxs = 1:obj.NData;
            end
            axes(axH);
            posSE = cellmap(@(i) obj.data(i).localizations(:,end-1), seqIdxs);
            names = cellmap(@(i) sprintf('[Seq:%i] %s',i,obj.data(i).sourceName),seqIdxs);
            BaseRPT.plot_distributions(axH,posSE, names, true);
            xlabel('Log-Likelihood');
            title('Emitter Log-Likelihood Distribution');
        end

        function plotEmitterPosition(obj, axH, seqIdxs)
            if nargin == 1 || isempty(axH)
                f=figure();
                whitebg(f);
                axH = axes();
            end
            if nargin<3 || isempty(seqIdxs)
                seqIdxs = 1:obj.NData;
            end
            axes(axH);
            pos = cellmap(@(i) obj.data(i).boxCoords.covertToEmitterBoxCenterRelativeLocations(obj.data(i).localizations(:,[1,2])), seqIdxs);
            posX = cellmap(@(p) p(:,1), pos);
            names = cellmap(@(i) sprintf('X: [Seq:%i] %s',i,obj.data(i).sourceName),seqIdxs);
            BaseRPT.plot_distributions(axH,posX, names, true);
            hold('on');
            posY = cellmap(@(p) p(:,2), pos);
            names = cellmap(@(i) sprintf('Y: [Seq:%i] %s',i,obj.data(i).sourceName),seqIdxs);
            Hs = BaseRPT.plot_distributions(axH,posY, names, true);
            for h=Hs
                h.LineStyle='--';
            end
            bs = max(obj.Params.BoxSize);
            xlim([-floor(bs/2), ceil(bs/2)]);
            xlabel('Relative Box Position (pixels)');
            title('Emitter Poistion Error Distribution');
        end

        function plotEmitterPositionSE(obj, axH, seqIdxs)
            if nargin == 1 || isempty(axH)
                f=figure();
                whitebg(f);
                axH = axes();
            end
            if nargin<3 || isempty(seqIdxs)
                seqIdxs = 1:obj.NData;
            end
            axes(axH);
            posSE = cellmap(@(i) mean(obj.data(i).localizations(:,[6,7]),2), seqIdxs);
            names = cellmap(@(i) sprintf('[Seq:%i] %s',i,obj.data(i).sourceName),seqIdxs);
            BaseRPT.plot_distributions(axH,posSE, names, true);
            xlabel('Position Uncertainty (pixels)');
            title('Emitter Poistion Error Distribution');
        end

        function f=viewFramePairing(obj, seqIdx, frameIdx)
            if nargin<3
                frameIdx=1;
            end
            f = figure('Position',[10,10,1000,800]);
            whitebg(f);
            conColor = [0,1,0.3];
            cm = [hot(65);conColor];
            colormap(cm); 
            ax=axes();
            ax.SortMethod='ChildOrder';
            ax.YDir='reverse';
            ax.TickDir='out';
            ax.XMinorTick='on';
            ax.YMinorTick='on';
            ax.ZMinorTick='off';
            ax.ZTick=[];
            ax.Box='on';
            ax.BoxStyle='full';
            ax.Projection='Orthographic';
            grid(ax,'on');
            zSep = 30;
            localizationSigmaDrawRadius = 3; % Number of sigmas out to draw localization.
            S = BaseRPT.makeLocalizationSurface();
            D = obj.data(seqIdx);
            frames = D.sequence(:,:,frameIdx);
            overlayMask= D.boxCoords.overlayMask(); %This corresponds to the box outlines
            overlayMask = overlayMask(:,:,frameIdx);
            leftMask = overlayMask(:,1:obj.SensorXSplit);
            rightMask = overlayMask(:,obj.SensorXSplit+1:end);
            leftIm = fix(64*mat2gray(frames(:,1:obj.SensorXSplit)));
            rightIm = fix(64*mat2gray(frames(:,obj.SensorXSplit+1:end)));
            leftAlpha = 0.5*double(leftMask>0);
            rightAlpha = 0.5*double(rightMask>0);
            sz=size(leftIm);
            axis([0, obj.SensorXSplit, 0, sz(1), 0, zSep]);
            pbaspect([sz(2)/sz(1), 1, zSep/sz(2)]);
            view(0,90);
            BX = repmat([0, sz(2)],2,1);
            BY = repmat([0; sz(1)],1,2);
            BZ = zeros(2,2);
            hold('on');
            ax.ALim = [0,1]; % Set the alpha axis to 0 to 1.
            surface(BX,BY,BZ,leftIm,'AlphaData',leftAlpha,'FaceAlpha','texturemap','FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            surface(BX,BY,BZ+zSep,rightIm,'AlphaData',rightAlpha,'FaceAlpha','texturemap','FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            transF = obj.generateChannelMap_null();
            lpts = D.refPts{frameIdx}{1};
            rpts = transF(D.refPts{frameIdx}{2});
            se = D.refPtsSE{frameIdx};
            
            for i=1:size(lpts,1)
                plot3([lpts(i,1), rpts(i,1)],[lpts(i,2), rpts(i,2)],[0,zSep],'LineStyle','-','Color',conColor,...
                      'MarkerSize',3,'MarkerEdgeColor','none');
                %draw localizations
                BaseRPT.drawLocalizationSurface(lpts(i,:), localizationSigmaDrawRadius*se{1}(i,:), 0, 66, S);
                BaseRPT.drawLocalizationSurface(rpts(i,:), localizationSigmaDrawRadius*se{2}(i,:), zSep, 66, S);
            end
            xlabel('X (px)');
            ylabel('Y (px)');
            title('Registration Map (Ch. 1 Coords)');
        end

        function f=viewMapDisplacement(obj, mapIdxs, testDataIdx)
            % mapList - cell array of mapping names to consider.  defaults to all mappings
            if nargin<2 || isempty(mapIdxs)
                mapIdxs = 1:obj.NMaps;
            else
                mapIdxs = mapIdxs(:)';
            end
            if nargin<3
                testDataIdx = obj.maps(mapIdxs(1)).testDataIdxs;
            end
            testData = obj.collectRefPts(testDataIdx);
            f=figure();
            whitebg(f);
            ms=7;
            ax=axes();
            hold(ax);
            xlim([0,obj.SensorXSplit]);
            ylim([0,obj.SensorSizeY]);
            rectangle('Position',[0,0,obj.SensorXSplit, obj.SensorSizeY],'EdgeColor',[0 0 0]);
            axis('equal');
            ax.YDir='reverse';
            ax.Box='on';
            ax.BoxStyle='full';
            grid('on');
            grid('minor');
            plot(testData{1}(:,1), testData{1}(:,2),'ro','MarkerSize',ms,'DisplayName','Channel1')
            styles = {'r+','g+','w+','b+','mx','gx','rx','bx'};
            for n=mapIdxs
                M = obj.maps(n);
                ch2 = M.mapFunctionPixels(testData{2});
                desc = sprintf('Channel2[%s]: RMSE:%.3g (nm) MaxErr:%.3g (nm)',M.name, M.RMSE, M.maxError);
                plot(ch2(:,1), ch2(:,2), styles{n}, 'MarkerSize',ms,'DisplayName', desc);
            end            
            xlabel('X position (Channel1) [px]');
            ylabel('Y position (Channel1) [px]');
            H=legend('location','best');
            H.FontSize=11;
            title('Registration Map Comparison');
        end

        function plotMapError(obj, axH, mapIdx, testDataIdx)
            if nargin == 1 || isempty(axH)
                f=figure();
                whitebg(f);
                f.Color=[0 0 0];
                axH = axes();
            end
            axes(axH);
            if nargin<3
                [~,mapIdx] = min([obj.maps.RMSE]);
            end
            if nargin<4
                testDataIdx = obj.maps(mapIdx).testDataIdxs;
            end
            testData = obj.collectRefPts(testDataIdx);
            ch1 = testData{1};
            interpMethod = 'natural';
 
            M = obj.maps(mapIdx);
            err = obj.measureMapError(M, testData);
            rmse = sqrt(mean(err(:).^2));
            max_err =max(err(:));
            
%             all_max_err = max([obj.maps(~strcmp({obj.maps(:).algorithm},'Null')).maxError]);
            all_max_err = max_err;
            hold('on');
            axH.YDir='reverse';
            axH.Box='on';
            axH.BoxStyle='full';
            grid('on');
            grid('minor');
            [Xq,Yq] = meshgrid(0:obj.SensorXSplit,0:obj.SensorSizeY);
            Vq = griddata(ch1(:,1), ch1(:,2), err, Xq, Yq, interpMethod);
            surface(Xq,Yq,Vq,'EdgeColor','none');
            xlim([0,obj.SensorXSplit]);
            ylim([0,obj.SensorSizeY]);
            zlim([0, all_max_err]);
            caxis([0, all_max_err]);
            aspectRatio = obj.SensorXSplit/obj.SensorSizeY;
            pbaspect(axH,[aspectRatio, 1, 0.4*min(1,aspectRatio)]);
            view(75,30);
            colormap('jet');
            H = colorbar();
            H.Label.String = sprintf('Error (nm)');
            axH.Projection='Orthographic';
            axH.TickDir='out';    
            xlabel('X position (Channel1) [px]');
            ylabel('Y position (Channel1) [px]');
            zlabel(sprintf('Error[%s] (nm)',M.name));
            title(sprintf('Map - %s - RMSE: %.3g(nm) MaxError: %.3g(nm)',M.name,rmse,max_err));
        end

        function f = viewMapVectors(obj, mapIdx, testDataIdx)
            vecScale=10;
            nullM = obj.computeMap('Null');
            M = obj.maps(mapIdx);
            if nargin<3
                testDataIdx = M.testDataIdxs;
            end
            testData = obj.collectRefPts(testDataIdx);
            ch1 = testData{1};
            ch2_null = nullM.mapFunctionPixels(testData{2});
            ch2_map = M.mapFunctionPixels(testData{2});
            f = figure();
            whitebg(f);
            ax=axes();
            hold(ax);
            xlim([0,obj.SensorXSplit]);
            ylim([0,obj.SensorSizeY]);
            rectangle('Position',[0,0,obj.SensorXSplit, obj.SensorSizeY],'EdgeColor',[0 0 0]);
            axis('equal');
            ax.YDir='reverse';
            ax.Box='on';
            ax.BoxStyle='full';
            grid('on');
            grid('minor');
            xlabel('X position (Ch.1) [px]');
            ylabel('Y position (Ch.1) [px]');
            
            ms=4;
            vec_null = vecScale*(ch2_null-ch1);
            vec_map = vecScale*(ch2_map-ch1);
            plot(ch1(:,1),ch1(:,2),'rs','MarkerSize',ms,'DisplayName','Ch.1');
            desc = sprintf('Ch.2[%s]: RMSE:%.3g (nm) MaxErr:%.3g (nm)', nullM.name, nullM.RMSE, nullM.maxError);
%             plot(ch2_null(:,1),ch2_null(:,2),'g+','MarkerSize',ms,'DisplayName',desc);
            desc = sprintf('Ch.2[%s]: RMSE:%.3g (nm) MaxErr:%.3g (nm)',M.name, M.RMSE, M.maxError);
            plot(ch2_map(:,1),ch2_map(:,2),'+','MarkerFaceColor',[0.4, 0.4, 1],'MarkerSize',ms,'DisplayName',desc);
            desc=sprintf('Ch.2[%s] Displacement[magnification: x%i]','Null',vecScale);
%             quiver(ch1(:,1), ch1(:,2), vec_null(:,1), vec_null(:,2),'DisplayName',desc,'Color','g','AutoScale','off','MaxHeadSize',0.05);
            desc=sprintf('Ch.2[%s] Displacement[magnification: x%i]',M.name,vecScale);
            quiver(ch1(:,1), ch1(:,2), vec_map(:,1), vec_map(:,2),'DisplayName',desc,'Color',[0.4, 0.4, 1],'AutoScale','off');
            H=legend('location','best');
            H.FontSize=11;
            title('Channel Displacement Vectors');            
        end
        



        function transFunc = generateChannelMap_null(obj, translation)
            % Generate the default null channel mapping.
            if nargin==1
                translation = [obj.SensorXSplit,0];
            end
            transFunc = @(p) p - repmat(translation,size(p,1),1);
        end
        
        function refData = collectRefPts(obj, dataIdxs)
            % [in]
            %   dataIdxs - a cell array of indexes.  Each index is size 1x2 and specifies the sequenceIdx and
            %   frameIdx of the reference points to include in the reference data set.
            % [out] refData - 2x1 cell array of reference points for each channel.  Each cell is the reference
            % points for that channel as Nx2 arrays
            refDataCells = cellmap(@(idx) obj.data(idx(1)).refPts{idx(2)}, num2cell(dataIdxs,2));
            if any(cellfun(@isempty,refDataCells))
                error('RegistrationAnalysis:collectRefPts','Got dataIdxs to non-existant data');
            end
            refData{1} = cell2mat(cellmap(@(data) data{1}, refDataCells));
            refData{2} = cell2mat(cellmap(@(data) data{2}, refDataCells));
        end

        function error = measureMapError(obj, map, testData)
            % Return errors in nm
            ch1 = testData{1};
            ch2 = map.mapFunctionPixels(testData{2});
            error = 1e3*sqrt(sum(((ch1-ch2).*repmat(obj.pixelSize,size(testData{2},1),1)).^2,2));
        end

    end

    methods %Dependent property definitions
        function n = get.NData(obj)
            % The number of fiducial sequences.
            n = numel(obj.data);
        end

        function n = get.NMaps(obj)
            % The number of maps being maintained
            n = numel(obj.maps);
        end
        
        function x = get.SensorXSplit(obj)
            % The X position at which the sensor is split.  This is taken to be 1/2 the SensorSizeX.
            if isempty(obj.frameSize)
                x=[];
            else
                x = obj.frameSize(2)/2; %SensorSizeX is the #columns in a frame which is second index in frameSize
            end
        end

        function x = get.SensorSizeX(obj)
            % The size of the sensor in the X direction (the readout direction).  This is the same as the
            % number of columns in the frames so this corresponds to obj.frameSize(2)
            if isempty(obj.frameSize)
                x=[];
            else
                x = obj.frameSize(2);
            end
        end

        function y = get.SensorSizeY(obj)
            % The size of the sensor in the Y direction (the non-readout direction).  This is the same as the
            % number of rows in the frames so this corresponds to obj.frameSize(1)
            if isempty(obj.frameSize)
                y=[];
            else
                y = obj.frameSize(1);
            end
        end

        function roi = get.ChannelROI(obj)
            if isempty(obj.frameSize)
                roi = [];
            else
                roi = {[1 obj.SensorXSplit 1 obj.frameSize(1)], [obj.SensorXSplit+1 obj.frameSize(2) 1 obj.frameSize(1)]};
            end
        end
        
        function val = get.calibrated(obj)
            val = logical(~isempty(obj.pixelSize) && ~isempty(obj.CCDGain) && ~isempty(obj.CCDBackground));
        end
        function val = get.hasValidData(obj)
            val = ~isempty(obj.data) && any(obj.data(:).valid);
        end
        function val = get.hasValidMap(obj)
            val = ~isempty(obj.maps) && any(obj.maps(:).valid);
        end
 
    end %Dependent property definitions

    methods (Static=true)
        function tform_func = generateChannelMap_fitgeotrans(refPts, method, varargin)
            % This uses matlab's fitgeotrans image coordinate registration to create an
            % channel2->channel1 coordinate mapping using reference point pairs.  Note we always work in
            % absolute coordinates.
            % 
            % This function returns another function that performs the actual mapping.
            % 
            % [in]
            %   refPts - 2x1 cell array of paired reference points for the two channels. (absolute coords)
            %   method - the method for fitgeotrans and associated arguments;
            % [out]
            %   tform_func - A function that takes in ch2Points in a Nx2 matrix and returns the translated
            %                positions in ch1 coords
            if nargin==1
                method='lwm';
                varargin={16};
            end
            tform = fitgeotrans(double(refPts{1}),double(refPts{2}),method, varargin{:});
            tform_func  = @(ch2Points) tform.transformPointsInverse(ch2Points);            
        end
        
        function tform_func = generateChannelMap_localAffine(refPts, nFitPoints, nWeightPoints)
            % This uses the localy weighted affine mapping coordinate registration to create an
            % channel2->channel1 coordinate mapping using reference point pairs.  Note we always work in
            % absolute coordinates.  And channel 2 is always on the right half of the sensor image.
            % 
            % This function returns another function that performs the actual mapping.
            % 
            % [in]
            %   refPts - 2x1 cell array of paired reference points for the two channels. (absolute coords)
            %   nFitPoints - [Default=9] The number of fitting points to use in estimating the local affine tforms.
            %   nWeightPoints - [Default=4] The number of local weighting points in use when combining local affine maps for
            %                    final weighting of the returned mapping
            % [out]
            %   tform_func - A function that takes in ch2Points in a Nx2 matrix and returns the translated
            %                positions in ch1 coords
            if nargin<3
                nWeightPoints=4;
            end
            if nargin<2
                nFitPoints=9;
            end
            affineTFs = RegistrationAnalysis.findAffineTform(refPts, nFitPoints);
            tform_func  = @(ch2Points) RegistrationAnalysis.affineLocalTform(ch2Points, affineTFs, refPts{2}, nWeightPoints);
        end

        function tform_func = generateChannelMap_smoothAffine(refPts, nFitPoints, distanceWeightExponent)
            % This uses the smoothed localy weighted affine mapping coordinate registration to create an
            % channel2->channel1 coordinate mapping using reference point pairs.  Note we always work in
            % absolute coordinates.  And channel 2 is always on the right half of the sensor image.
            % 
            % This function returns another function that performs the actual mapping.
            % 
            % [in]
            %   refPts - 2x1 cell array of paired reference points for the two channels. (absolute coords)
            %   nFitPoints - [Defulat=9] The number of fitting points to use in estimating the local affine tforms.
            %   distanceWeightExponent - [Default=2] The exponent to weight the contribution of the mappings with
            %                                      based on distance.  Larger values correspond to more
            %                                      concentrated local weighting
            % [out]
            %   tform_func - A function that takes in ch2Points in a Nx2 matrix and returns the translated
            %                positions in ch1 coords
            if nargin<2
                nFitPoints = 9;
            end
            if nargin<3
                distanceWeightExponent = 2;
            end
            affineTFs = RegistrationAnalysis.findAffineTform(refPts, nFitPoints);
            tform_func  = @(ch2Points) RegistrationAnalysis.affineSmoothTform(ch2Points, affineTFs, refPts{2}, distanceWeightExponent);
        end
        
        function t_points = affineLocalTform(points,tForms,ch2RefPoints,nWeightPoints)
            % This uses a weighted average of affine transforms from a fixed number of local points.  
            % There are Nref, reference points each of which has a computed 3x3 affine transform.
            % For each point to transform we select the nWeightPoints closes reference points and weight thier
            % tranformations by 1/distance.
            % 
            % [in] points - size:[Npoints,2] vector of points from channel2 to transform.  Column1=x Column2=y.
            % [in] tForms - size:[3,3,NRefPoints] array of 3x3 transform matricies.
            %              where the M matrix has the structure:
            %
            %               [a b c][x] [xp]
            %               [d e f]|y]=[yp]
            %               [0 0 1][1] [1 ]
            %              Each 3x3 matrix corresponds to a reference point, and is assumed to be locally accurate.
            % [in] ch2RefPoints - size:[2,NRefPoints] The locations in channel2 corresponding to each affine 
            %                                      transform
            % [in] nWeightPoints - integer>=1 - The number of nearby points to include in the weighted average.
            %                                 the nWeightPoints nearest reference points to each point will be used
            %                                 to compute the transform. Using a 1/distance weighting. [Default=4]
            % [out] t_points - size:[Npoints,2] - The transformed point locations.
            if nargin<4
                nWeightPoints=4;
            end
            Npoints = size(points,1);
            X = [points,ones(Npoints,1)]'; %points extended with 1's for general affine tForms
            t_points = zeros(3,Npoints);
            for n=1:Npoints
                dists = sqrt((points(n,1) - ch2RefPoints(:,1)).^2 + (points(n,2) - ch2RefPoints(:,2)).^2);
                [~, sidx] = sort(dists);
                if dists(sidx(1))==0 %This is exactly at a reference points so we can't weight by distance.
                    t_points(:,n) =tForms(:,:,sidx(1)) * X(:,n);
                else
                    nearby_ref_idxs = sidx(1:nWeightPoints)'; % nearby_ref_idxs is the indexs of the nWeightPoints closest reference points
                    weight_norm = 1./sum(1./dists(nearby_ref_idxs));
                    for i = nearby_ref_idxs
                        t_points(:,n) = t_points(:,n) + (weight_norm * (1/dists(i)) .* (tForms(:,:,i) * X(:,n)));
                    end
                end
            end
            t_points = t_points(1:2,:)';
        end
        
        function t_points = affineSmoothTform(points,tForms,refPoints,distanceWeightExponent)
            %AffineTform Affine Transform of a set of points
            % This performs a locally weighted affine transform that is globablly smooth, in that we use all
            % reference points in the weighting of the transform.  This eliminates discontinuities in the effective
            % tranformation.  We weight by 1/distance^k.  Where we call k the 'distanceWeightExponent' which has a
            % defult value of 2, but could be any number >=1.  Larger numbers will give more weight to closer
            % reference points.
            %
            % [in] points - size:[Npoints,2] vector of points from channel2 to transform.  Column1=x Column2=y.
            % [in] tForms - size:[3,3,NRefPoints] array of 3x3 transform matricies.
            %              where the M matrix has the structure:
            %
            %               [a b c][x] [xp]
            %               [d e f]|y]=[yp]
            %               [0 0 1][1] [1 ]
            %              Each 3x3 matrix corresponds to a reference point, and is assumed to be locally accurate.
            % [in] refPoints - size:[2,NRefPoints] The locations that each affine transform was caluclated at and
            %                                      is assumed to be locally accurate at.
            % [in] distanceWeightExponent - scalar>0 - The exponent to raise 1/distance^distanceWeightExponent where
            %                                          larger numbers will lead to more locally-weighted transformations.
            %                                          [Default=2]
            % [out] t_points - size:[Npoints,2] - The transformed point locations.
            if nargin<4
                distanceWeightExponent=2;
            end
            Npoints = size(points,1);
            NrefPoints = size(refPoints,1);
            X = [points,ones(Npoints,1)]'; %points extended with 1's for general affine tForms
            t_points = zeros(3,Npoints);
            for n=1:Npoints
                dists = sqrt((points(n,1) - refPoints(:,1)).^2 + (points(n,2) - refPoints(:,2)).^2);
                zidx = find(dists==0);
                if(numel(zidx)>1)
                    warning('Reference points not unique');
                end
                if numel(zidx)==1 %This is exactly at a reference points so we can't weight by distance.
                    t_points(:,n) = tForms(:,:,zidx) * X(:,n);
                else
                    weights = dists.^(-distanceWeightExponent);
                    weight_norm = 1/sum(weights);
                    for i = 1:NrefPoints
                        t_points(:,n) = t_points(:,n) + (weight_norm * weights(i) * (tForms(:,:,i) * X(:,n)));
                    end
                end
            end
            t_points = t_points(1:2,:)';
        end
        
        function M = findAffineTform(refPts, nFitPoints)
            %FindAffineTform Find Affine Transformation Matrix
            %   FindAffineTform  uses the closet nRefPoints to each point in finds the Affine transform of X to Xp
            %   by M*ch2_pts = ch1_pts 
            %   [in]
            %       refPts - 2x1 cell array of paired reference points for the two channels. (absolute coords)
            %       nFitPoints - The number of points to use in the local estimation of the local affine map [default=9]     
            if nargin<2
                nFitPoints = 9;
            end            
            nPoints = size(refPts{1},1);
            M = zeros(3,3,nPoints);
            for i = 1:nPoints
                dist_sq = sum((refPts{2} - repmat(refPts{2}(i,:),nPoints,1)).^2,2);
                [~, sidx] = sort(dist_sq);
                close_idxs = sidx(1:nFitPoints);
                invA = pinv([refPts{1}(close_idxs,:), ones(nFitPoints,1)]);
                M(:,:,i)= pinv([(invA*refPts{2}(close_idxs,:))'; 0, 0, 1]);
            end            
        end
        
    end %public static methods

    methods (Access=protected)
        function [trainDataIdxs, testDataIdxs] = getTrainTestDataIdxs(obj)
            trainDataIdxs = [];
            testDataIdxs = [];     
            if obj.NData<1 || ~obj.calibrated
                return;
            end                
            % Look for data with either 1 or 2 good frames with valid reference points.
            firstSeqIdxs=[];
            secondSeqIdxs=[];
            for k=1:obj.NData
                [N,I] = sort(obj.data(k).NRefPts,1,'descend');
                if N(1)>0
                    % These will be [dataIdx, frameIdx] for the frame with the most valid datapoints from
                    % each dataset.  Datasets with no valid reference points will be excluded;
                    firstSeqIdxs=[firstSeqIdxs; k, I(1)]; %#ok<AGROW>
                end
                if numel(N)>1 && N(2)>0
                    % These will be [dataIdx, frameIdx] for the frame with the *second* most valid datapoints from
                    % each dataset.  Datasets with no valid reference points will be excluded;
                    secondSeqIdxs=[secondSeqIdxs; k, I(2)]; %#ok<AGROW>
                end
            end
            if isempty(firstSeqIdxs)
                % No useful frames of data
                return          
            elseif isempty(secondSeqIdxs)
                NSingleFrameSeq = size(firstSeqIdxs,1);
                if NSingleFrameSeq==1
                    %only one set of reference points.  This is not optimal.  We must train and test
                    %with the same data.
                    trainDataIdxs = firstSeqIdxs;
                    testDataIdxs = firstSeqIdxs;
                else
                    %split data into train and test datasets
                    Ntrain = ceil(NSingleFrameSeq/2);
                    trainDataIdxs = firstSeqIdxs(1:Ntrain,:);
                    testDataIdxs = firstSeqIdxs(Ntrain+1:end,:);
                end
            else
                %Use the first frame of any dataset to train, use any datasets with second frames to test.
                trainDataIdxs = firstSeqIdxs;
                testDataIdxs = secondSeqIdxs;
            end
        end

        function reprocessData(obj, dataIdxs)
            if obj.NData < 1 || ~obj.calibrated
                return;
            end
            if nargin<2
                dataIdxs = 1:obj.NData;
            end
            obj.updateWaitbar(0,'Reprocessing Data...');
            N = numel(dataIdxs);
            for n=dataIdxs
                obj.updateWaitbar((n/N)*.9,'Reprocessing Data ...');
                obj.data(n) = obj.processData(obj.data(n));
            end
            obj.reprocessMaps();
        end

        function reprocessMaps(obj, idxs)
            obj.updateWaitbar(0.9,'Reprocessing Maps...');
            if nargin < 2
                idxs = 1:obj.NMaps;
            else
                idxs = idxs(:)';
            end
            if isempty(idxs)
                for n=1:numel(obj.MapAlgorithms)
                    obj.appendMap(obj.computeMap(obj.MapAlgorithms{n}));
                end
            else
                for n=idxs
                    obj.maps(n) = obj.computeMap(obj.maps(n).algorithm, obj.maps(n).params);
                end
            end            
            obj.updateWaitbar(1);
        end

        function generateChannelMaps(obj)
            % Select training and testing datasets and generate all mappings with default parameter values
            trainDataIdxs = cellmap(@(seqIdx) [seqIdx, 1], 1:obj.NData);
            testSeqs = find(obj.NFrames>1); % will test on second frame of any sequences with more than 1 frame
            if isempty(testSeqs)
                %couldn't find any sequences with a second frame must use last train sequence for testing
                testDataIdxs = trainDataIdxs(end);
                trainDataIdxs(end) = [];
            else
                testDataIdxs = cellmap(@(seqIdx) [seqIdx, 2], testSeqs);
            end
            for alg = obj.MapAlgorithms
                obj.updateMapAlgorithm(alg{1}, alg{1}, [], trainDataIdxs, testDataIdxs);
            end
        end

        function [RMSE, maxError] = testChannelMap(obj, mapF, testData)
            % Return errors in nm
            ch1 = testData{1};
            ch2 = mapF(testData{2});
            nPts = size(ch1,1);
            error = sqrt(sum((repmat(obj.pixelSize,nPts,1).*(ch1-ch2)).^2,2));
            RMSE = 1e3*sqrt(mean(error.^2));
            maxError = 1e3*max(error);
        end

        function cal_im = calibrateImage(obj,im)
            %Gain calibrate a raw CCD image, using CCDGain CCDBackground properties
            % [in]
            %   im - 2d or 3d image to be calibrates
            % [out]
            %   cal_im - calibrated image normalized to be positive
            cal_im = single(obj.CCDGain*(im - obj.CCDBackground));
            cal_im(:) = max(0, cal_im(:));            
        end

        function checkCalibrated(obj)
            % Check the object is calibrated.  When the object is calibrated it has all the necessary 
            % information to begin processing the data.
            if isempty(obj.pixelSize)
                error('RegistrationAnalysis:NotCalibrated', 'Pixel Size must be set to complete calibration.');
            end
            if ~obj.calibrated
                error('RegistrationAnalysis:NotCalibrated', 'Data is not calibrated.');
            end
        end

        function checkValidData(obj)
            % Check the the object is calibrated and has valid data avalible
            obj.checkCalibrated();
            if obj.NData<1
                error('RegistrationAnalysis:NoData', 'No valid data has been loaded.');
            end
        end

        function setWorkingDir(obj, dataPath)
            if ~obj.initialized
                [wdir, baseName, ~] = fileparts(dataPath);
                obj.workingDir = wdir;
                obj.Paths.saveFile = [baseName obj.saveFileExt];
            end
        end

        function loadRA(obj, raFilePath)
            [dirpath,~,~] = fileparts(raFilePath);
            S=load(raFilePath,'-mat');
            if obj.initialized
                obj.savePreservedProperties();
                pp=obj.preservedProperties;
                obj.resetObject();
                obj.preservedProperties=pp;
            end
            obj.copyobj(S.obj);
            obj.workingDir = collapsepath(dirpath);
            obj.initialized = true;
            obj.dirty = false;
        end

        function loadHMMFiducial(obj,dataPath)
            S = load(dataPath, '-mat');
            trn = single(S.trn);
            tst = single(S.tst);
            sz = size(trn);
            if obj.NData<1 || isempty(obj.frameSize)
                obj.frameSize = sz(1:2);
            end
            [~,baseName,~] = fileparts(dataPath);
            dateTime = obj.parseFileNameDateTime(dataPath);
            obj.appendData(obj.makeData(trn, ['trn#' baseName], dateTime));
            obj.appendData(obj.makeData(tst, ['tst#' baseName], dateTime));
            obj.setWorkingDir(dataPath);
            obj.initialized=true;
        end

        function loadSequence(obj,dataPath)
            % Load a .mat file with a sequence variable
            S = load(dataPath, '-mat');
            seq = single(S.sequence);
            sz = size(seq);
            if obj.NData<1 || isempty(obj.frameSize)
                obj.frameSize = sz(1:2); % dip_image sizes are reported backwards  set frameSize to [sizeY, sizeX]
            end
            [~,baseName,~] = fileparts(dataPath);
            dateTime = obj.parseFileNameDateTime(dataPath);
            obj.appendData(obj.makeData(seq, baseName, dateTime));
            obj.setWorkingDir(dataPath);
            obj.initialized=true;
        end

        function loadLidkeFormat2(obj,dataPath)
            %This is yet anouther new .mat format that is being produced now (Circa 2016).  Need to ask Keith about
            %this.  Argh...
            S = load(dataPath, '-mat');
%             P = S.Params;
            warning('RegistrationAnalyis:BadFormat','LidkeFormatV2 is not supported.');
%             seq = permute(single(S.Data),[2 1 3]); % Not sure if I should flip it or not here. 
            seq = single(S.Data); % Not sure if I should flip it or not here.  
            if obj.NData<1 || isempty(obj.frameSize)
                sz=size(seq);
                obj.frameSize = sz(1:2);
            end
            [~,baseName,~] = fileparts(dataPath);
            dateTime = obj.parseFileNameDateTime(dataPath);
            obj.appendData(obj.makeData(seq, baseName, dateTime));
            obj.setWorkingDir(dataPath);
            obj.initialized=true;
        end

        function loadSPTRegisterChannels(obj,dataPath)
            % Load a .mat file with a sequence variable
            S = load(dataPath, '-mat');
            P = S.Params;
            if obj.NData<1 || isempty(obj.pixelSize)
                obj.pixelSize = double([P.PixelSzX, P.PixelSzY]);
            end
            if obj.NData<1 || isempty(obj.frameSize)
                obj.frameSize = double([P.YPixels, P.XPixels]); % Not sure about the orientation here with respect to X&Y
            end
            if obj.NData<1 || isempty(obj.CCDBackground)
                obj.CCDBackground = double(P.CCDoffset); % This is normallay a fixed default value and very possibly is wrong
            end
            if obj.NData<1 || isempty(obj.CCDGain)
                obj.CCDGain = 1/double(P.Gain); % Convert units. This is normallay a fixed default value and very possibly is wrong
            end
            seq = permute(single(P.Data),[2 1 3]);
            [~,baseName,~] = fileparts(dataPath);
            dateTime = obj.parseFileNameDateTime(dataPath);
            obj.appendData(obj.makeData(seq, baseName, dateTime));
            obj.setWorkingDir(dataPath);
            obj.initialized=true;
        end

        function loadChannelRegistrationV1(obj,dataPath)
            S = load(dataPath, '-mat');
            reg = S.CRobj;
            %This ChReg data type did not save the raw images which RegistrationAnalysis works with
            %So we attempt to reconstruct a fiducial image from the saved boxes...
            im = RegistrationAnalysis.reconstructChannelRegistrationV1Image(reg);
            sz = size(im);
            if obj.NData<1 || isempty(obj.frameSize)
                obj.frameSize = sz(1:2);
            end
            if obj.NData<1 || isempty(obj.pixelSize)
                obj.pixelSize = [reg.PixelSize, reg.PixelSize];
            end
            if obj.NData<1 || isempty(obj.CCDGain)
                obj.CCDGain = 1./reg.CCDgain;
            end
            if obj.NData<1 || isempty(obj.CCDBackground)
                obj.CCDBackground = reg.CCDoffset;
            end
            im = im/obj.CCDGain + obj.CCDBackground; % Purposely re-encode the RAW image so everything works nicely.
            [~,baseName,~] = fileparts(dataPath);
            dateTime = obj.parseFileNameDateTime(dataPath);
            obj.appendData(obj.makeData(im, ['CRV1#' baseName], dateTime));
            obj.setWorkingDir(dataPath);
            obj.initialized = true;
        end
    end % protected methods

    methods (Access=protected, Static=true)
        function dateTime = parseFileNameDateTime(filePath)
            %Try to parse date from file name
            [~,baseName,~] = fileparts(filePath);
            parts=strsplit(baseName,'-');
            if numel(parts)==6
                dateTime = datetime([cellfun(@str2double,parts(2:6)), 0]);
            elseif numel(parts)>6
                dateTime = datetime(cellfun(@str2double,parts(end-5:end)));
            else
                dateTime = [];
            end
        end

        function im = reconstructChannelRegistrationV1Image(CRobj)
           im = zeros(flip(CRobj.ImageSize));
           boxIdx = (1:CRobj.BoxSize)-1;
           for ii = 1:size(CRobj.BoxStartPixel,1)
                for jj = 1:size(CRobj.BoxStartPixel,2)
                    im(CRobj.BoxStartPixel(ii,jj,2)+boxIdx,...
                        CRobj.BoxStartPixel(ii,jj,1)+boxIdx) = CRobj.ImageStack(:,:,ii-1,jj-1);
                end
            end
        end
    end %protected static methods

    methods (Access=protected) % Abstract methods inherited from Pickle            
       function val = getProtectedProperty(obj, name)
           %This is necessary for Pickle functionality to be able to access subclass protected variables
           val = obj.(name);
       end
    
       function modifyProtectedProperty(obj, name, newval)
           %This is necessary for Pickel functionality to be able to change subclass protected variables
           obj.(name)=newval;
       end
    end % Abstract methods inherited from Pickle
end

