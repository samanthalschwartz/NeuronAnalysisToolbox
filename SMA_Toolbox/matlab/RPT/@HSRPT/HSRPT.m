% HSRPT.m
% Mark J. Olah (mjo@cs.unm.edu)
% 02/15
%
% Hyper-Spectral Robust Particle Tracking.

classdef HSRPT < BaseRPT
   properties (Constant=true, Hidden=true)
        %Abstract properties inherited from Pickle
        saveFileExt = '.hsrpt'
        SaveableDataFormats = {'*.hsrpt', 'HSRPT (.hsrpt)'};
        LoadableDataFormats = {'*.hsdata;*.hsrpt','All Loadable Sources (.hssdata,.hsrpt)';...
                               '*.hsdata','HSData file (.hsrpt)'; '*.hsrpt', 'HSRPT file (.hsrpt)'};

        %Abstract properties inherited from BaseRPT

        %These methods allow the generic methods of BaseRPT to do alot of common work
        DataClass = 'HSData';
        DataFileExt = HSData.saveFileExt;

        % Emitters are an format mainly for internal use.  
        % The information is the similar to the Localizations, but we leave more information in pixel
        % coordinates and retain mapping to boxCoords index in boxIdx column.  This allows association
        % between the emitter information and the box information
        NEmitterColumns = 19;
        EmitterColumnNames = {'x', 'y', 'L', 'I', 'bg', 'sigmaX', 'sigmaY', 'sigmaL',...
                              'SE_x', 'SE_y', 'SE_L', 'SE_I', 'SE_bg', 'SE_sigmaX', 'SE_sigmaY', 'SE_sigmaL',...
                              'LLH','boxIdx','frameIdx'};
        EmitterColumnUnits = {'px', 'px', 'px', 'photons', 'photons/px', 'px','px', 'px', 'px','px','px', 'photons', 'photons/px', 'px','px','px','','index','index'};



        % Table/Matrix columns format for Localization results
        NLocalizationColumns = 17;
        LocalizationColumnNames = {'t', 'x', 'y', 'L', 'I', 'bg', 'sigmaX', 'sigmaY', 'sigmaL',... %1-8
                                   'SE_x', 'SE_y', 'SE_L', 'SE_I', 'SE_bg', 'SE_sigmaX', 'SE_sigmaY', 'SE_sigmaL','frame'}; %9-17
        LocalizationColumnUnits = {'s', 'um', 'um', 'nm', 'photons', 'photons/px', 'um', 'um', 'nm',...
                                   'um', 'um', 'nm', 'photons', 'photons/px', 'um', 'um', 'nm','index'};
        LocalizationColumnDescriptions = {'Time','x-Position estimate', 'y-Position estimate','Wavelength estimate',...
                                   'Intensity estimate','Mean background intensity per pixel estimate',...
                                   'Apparent gaussian sigma X estimate','Apparent gaussian sigma Y estimate','Apparent gaussian sigma L estimate',...
                                   'Standard error of x-position estimate','Standard error of y-position estimate','Standard error of Wavelength estimate',...
                                   'Standard error of intensity estimate', 'Standard error of background intensity estimate',...
                                   'Standard error of apparent gaussian sigma X estimate','Standard error of apparent gaussian sigma Y estimate',...
                                   'Standard error of apparent gaussian sigma L estimate',...
                                   'Frame index'};

        % Table/Matrix columns format for Track results
        NTrackColumns = 18;
        TrackColumnNames = {'t', 'x', 'y', 'L', 'I', 'bg', 'sigmaX', 'sigmaY', 'sigmaL',...
                            'SE_x', 'SE_y', 'SE_L', 'SE_I', 'SE_bg', 'SE_sigmaX', 'SE_sigmaY', 'SE_sigmaL','frame'};
        TrackColumnUnits = {'s', 'um', 'um', 'nm', 'photons', 'photons/px', 'um', 'um', 'nm',...
                                 'um', 'um', 'nm', 'photons', 'photons/px', 'um', 'um', 'nm', 'index'};        
        TrackColumnDescriptions = {'Time', 'x-Position estimate', 'y-Position estimate', 'Wavelength estimate',...
                                   'Intensity estimate', 'Mean background intensity per pixel estimate',...
                                   'Apparent gaussian sigma X estimate','Apparent gaussian sigma Y estimate','Apparent gaussian sigma L estimate',...
                                   'Standard error of x-position estimate','Standard error of y-position estimate','Standard error of Wavelength estimate',...
                                   'Standard error of intensity estimate', 'Standard error of background intensity estimate',...
                                   'Standard error of apparent gaussian sigma X estimate','Standard error of apparent gaussian sigma Y estimate',...
                                   'Standard error of apparent gaussian sigma L estimate',...
                                   'Frame index'};
        DefaultWindowSize=[800 750];
    end
    
    properties
        %Parameters controlling [phase=3]: The identification of fetures representing possible point
        % emitters, buy looking for local maxima in a filtered image.
        ParamsFindMaxima=struct(... %phaseIdx=3
            'method','DoG',... %Options ['LoG'=Laplacian of Gaussian, 'DoG'=Difference of Gaussian'
            'filterSigmas',[ 3.0  4.8 5.8;... %L
                            1.14 1.3 2.5;... %Y
                            0.9  1.1 2.0],...%X
            'maximaNeighborhoodSize',7, ... %Options odd number >=3
            'scaleNeighborhoodSize',7 ... %Options odd number >=3
            );

        %Parameters controlling [phase=4]: The filtering of maxima to remove the background peaks and
        % form a set of candidate imitter boxes or sub-images that are to be fit in next stage.
        ParamsFilterMaxima=struct(...%phaseIdx=4
            'maximaThreshold',-1, ... %If empty or negative we will estimate this
            'optimalBoxSize',[10 12 14; 8 9 9; 7 8 9]...
            );

        %Parameters controlling [phase=5]: The localization (fitting) of the
        % identified candidate emitters from the boxes identified in phase 4
        ParamsLocalizeEmitters=struct(...%phaseIdx=5
            'model','GaussHSsMAP',... %The class name of the emitter model to use
            'estimator','Newton'...   %The estimation technique to use
        );

        %Parameters controlling [phase=6]: The filtration of the fitted
        %emitters for quality.
        ParamsFilterEmitters=struct(...%phaseIdx=6
            'minIntensity',300, ...
            'minSigma',[0.50, 0.8, 1.0],... %Sigma is allowed to go as small as 0.5 pixels.
            'maxSigma',[3.0, 5.0, 8.0],...
            'maxPositionSE',0.3,... %Maximum SE in either of the position arguments as computed by sqrt(crlb)
            'maxLambdaSE', 0.6,... %Maximum SE in either of the position arguments as computed by sqrt(crlb)
            'certVsUniformModel', -1,... % <1  0.95 = 95% certainty the emiitter model is correct;
            'certVsNoiseModel', -1, ...
            'overlapDistance',2.5, ... %Overlap 
            'overlapLambdaDistance',3.5 ... %Overlap 
            );
        
        %Parameters controlling [phase=7]: The tracking of the fitted
        %emitters.
        ParamsTrack=struct(...%phaseIdx=7
            'D',1.3,... % pixels^2/frame
            'Kon',0.1,... % 1/frame
            'Koff',0.1,... % 1/frame
            'LambdaR',0.3, ... %px
            'MaxSpeed',2.5,... % px/frame
            'MaxGapCloseFrames',20, ... %frames
            'MinGapCloseTrackLength',1,... 
            'MinFinalTrackLength',5 ...
            );  
    end

    properties (Dependent=true)       
        frameSize; %size [L Y X] of an idivudal frame from the ROI
        nFrames;
        ROIOrigin; %The 
        ROIPhysical;
        ROIPhysicalOrigin;
        ROIColorMap;
    end

    properties (Hidden=true)
        version=1; %For future file format version changes
    end

    methods
        function obj=HSRPT( varargin )
            % Input (options):
            %  (1) <empty>
            %  (2) hsrpt_filepath
            %  (3) hsdata_filepath, roi
            %  (4) hsdata_obj, roi
            %  
            % roi: (optional) allows the choice of ROI to be saved selected.   The roi variable if given can be:
            %   (I) integer index into data.ROI
            %   (II) 1x6 integer array with [xmin, xmax, ymin, ymax, Lmin, Lmax]
            %   (III) 1x8 integer array with [xmin, xmax, ymin, ymax, Lmin, Lmax, tmin, tmax]
            %
            % Option (1) makes an empty HSRPT with no .hsrpt filename
            % Option (2) Opens a saved HSRPT object from an .HSRPT file
            % Option (3) Makes a new HSRPT object using the data stored in the
            %           file path.  The new HSRPT will not have a filename until it is calling save() or saveas().
            %           The roi is optional and defaults to the whole frame.
            % Option (4) The same as option (4), except the data is given as an
            %            object instead of a filename.
            %
            % These argurments are the same as the load method which can be used
            % to reload a different HSRPT into an already created object.
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
               stats.FrameSizePixels = obj.frameSize([3,2]); %X Y
               stats.FrameSizeMicron = obj.frameSize([3,2])*obj.data.pixelSize; % X Y
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
            fprintf('Phase: Find Maxima\n');
            tic;

            obj.boxxer = Boxxer3D(obj.frameSize, P.filterSigmas);
            frames=obj.getFrames();
            R.filteredSumImage = HSData.makeRGBSumImage(obj.getFilteredFrames());
            
            [R.rawMaxima, R.rawMaximaVals] = obj.boxxer.scaleSpaceDoGMaxima(frames,...
                                                P.maximaNeighborhoodSize, P.scaleNeighborhoodSize); %[L Y X]           
            R.rawMaximaImage = obj.computeMaximaImage(R.rawMaxima, R.rawMaximaVals);
            obj.ResultsFindMaxima = R;

            obj.updateWaitbar(1);
            obj.times.findMaxima = toc;
            fprintf(' Time: %.3fs\n',obj.times.findMaxima);
            obj.setPhase(3);
        end

        function filterMaxima(obj)
            obj.checkPhase(3); 
            P = obj.checkParamsFilterMaxima();
            fprintf('Phase: Filter Maxima\n');
            obj.updateWaitbar(0,'Phase: Filter Maxima');
            tic;

            Maxima = obj.ResultsFindMaxima;
            if P.maximaThreshold<0;
                smax = sort(Maxima.rawMaximaVals,1,'descend');
                [~,R.maximaThreshold]=triThres(smax);
            else
                R.maximaThreshold=P.maximaThreshold;
            end
            obj.updateWaitbar(0.2);

            R.filter = Maxima.rawMaximaVals>=R.maximaThreshold;
            R.maxima = Maxima.rawMaxima(:,R.filter);
            R.maximaVals = Maxima.rawMaximaVals(R.filter);
            R.maximaImage = obj.computeMaximaImage(R.maxima, R.maximaVals);
            R.boxCoords = obj.boxxer.generateBoxCoords(R.maxima([1,2,3,5],:), obj.nFrames, P.optimalBoxSize);
            [R.emitterImages, R.emitterFrameIdx] = R.boxCoords.makeROI(obj.getFrames());
            obj.ResultsFilterMaxima = R;

            obj.updateWaitbar(1);
            obj.times.filterMaxima = toc;
            fprintf(' Time: %.3fs\n',obj.times.filterMaxima);
            obj.setPhase(4);
        end
        
        function localizeEmitters(obj)
            obj.checkPhase(4); 
            P = obj.checkParamsLocalizeEmitters();
            obj.updateWaitbar(0,'Phase: Localize Emitters');
            fprintf('Phase: Localize Emitters');
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
                positions = 0.5+double(boxCoords.boxCenter(:,idxs) - boxCoords.boxOrigin(:,idxs));
                thetaInit{n} = double([positions; zeros(2,numel(idxs)); scales]);
                if obj.emitterModel{n}.nParams==4
                    thetaInit{n} = thetaInit{n}(1:4,:); %For the 4-parameter models don't include the sigma-scales.
                end
%                 [etheta{n}, crlb{n}, llh{n}] = obj.emitterModel{n}.estimate(im_list{n}, P.estimator, thetaInit{n});
                [etheta{n}, crlb{n}, llh{n}] = obj.emitterModel{n}.estimate(im_list{n}, P.estimator);
                obj.updateWaitbar(0.1+0.8*n/N);
            end
            obj.ResultsLocalizeEmitters.rawTheta = etheta; %rawTheta is in cell-based format for easy re-use with the emitter Model
            obj.ResultsLocalizeEmitters.thetaInit = thetaInit; %save the theta init for later debugging.

            %Form localizations (expand sigma to SigmaX and SigmaY)
            sigmaX=obj.data.psf(1);
            sigmaY=obj.data.psf(2);
            etheta = [etheta{:}];
            etheta = [etheta(1:5,:); etheta(6,:) * sigmaX; etheta(6,:) * sigmaY; etheta(7,:)];
            crlb = sqrt([crlb{:}]);
            crlb = [crlb(1:5,:); crlb(6,:) * sigmaX; crlb(6,:) * sigmaY; crlb(7,:)];
            llh = vertcat(llh{:})';
            boxIdx = double([boxIdx{:}]);
            frameIdx = double([frameIdx{:}]) + obj.ROI(7)-1; 
            E = [etheta; crlb; llh; boxIdx; frameIdx]'; %internal emitter format
            shift = double(boxCoords.boxOrigin([3,2,1],boxIdx)')-1; %shift switches x/y since boxxer deals in row/col and locs are in x/y
            E(:,1:3) = E(:,1:3) + shift + repmat(obj.ROIOrigin,size(E,1),1); %correct for box coords

            %Sorted emitters by frame Idx to maintain relationship with localizations
            [~,sidx] = sort(E(:,end));
            E = E(sidx,:);

            obj.ResultsLocalizeEmitters.rawEmitters = E;            


%             P = obj.checkParamsLocalizeEmitters();
%             ims = obj.ResultsFilterMaxima.emitterImages;
%             frameIdx = double(obj.ResultsFilterMaxima.emitterFrameIdx{1});% + obj.ROI(7) -1;
%             
%             obj.initializeEmitterModel();
%             [etheta, crlb, llh] = obj.emitterModel.estimateMAP(ims, P.estimator);
%             sigmaX=obj.ParamsFindMaxima.filterSigma(3,1);
%             sigmaY=obj.ParamsFindMaxima.filterSigma(2,1);
%             E = [etheta(1:3,:) + double(obj.ResultsFilterMaxima.boxCoords.origin([3,2,1],:)-1);...
%                  etheta(4:5,:); etheta(6,:) * sigmaX; etheta(6,:) * sigmaY; etheta(7,:);...
%                  sqrt(crlb(1:5,:)); sqrt(crlb(6,:))*sigmaX; sqrt(crlb(6,:))*sigmaY; sqrt(crlb(7,:));...
%                  frameIdx];
%             obj.ResultsLocalizeEmitters.rawEmitters=E;
%             obj.ResultsLocalizeEmitters.rawTheta = etheta;
%             obj.ResultsLocalizeEmitters.LLH = llh;
            
            obj.updateWaitbar(1);
            obj.times.localizeEmitters = toc;
            fprintf(' Time: %.3fs\n',obj.times.localizeEmitters);
            obj.setPhase(5);
        end
        
        function initializeEmitterModel(obj)
            % make sure the emitterModel transient property is initialized
            boxCoords = obj.getBoxCoords();
            P = obj.ParamsLocalizeEmitters;
            emitterConstructor = str2func(P.model);
            size_cats = boxCoords.sizeCategories;
            Nsize_cats = boxCoords.NsizeCategories;
            obj.emitterModel = cellmap(@(i) emitterConstructor(size_cats([3,2,1],i), [obj.data.psf, 1]), 1:Nsize_cats);
        end
        
        function filterEmitters(obj)
            obj.checkPhase(5); 

            fprintf('Phase: Filter Emitters\n');
            tic;

            P = obj.checkParamsFilterEmitters();
            E = obj.ResultsLocalizeEmitters.rawEmitters;
            N = size(E,1); % number of raw localziations
            boxes = obj.getBoxCoords(); % the box coords object
            filter = true(N,1); %true if we will keep the localizations
            if ~isempty(P.minIntensity)
                filter = filter & E(:,4)>=P.minIntensity;
            end
            %dist filter
            min_dist_radius=2;
            min_Ldist_radius=3;
            [~,sidx]=sort(E(:,4),2,'descend');
            for i=sidx(1:N)
                if(filter(i))
                    dist = sqrt((E(:,1)-E(i,1)).^2+(E(:,2)-E(i,2)).^2);
                    Ldist = abs(E(:,3)-E(i,3));
                    dist_filter= ~((dist<min_dist_radius) & (Ldist<min_Ldist_radius) & (E(:,end)==E(i,end)));
                    dist_filter(i)=true;
                    filter =filter & dist_filter;
                    nFilt = sum(~dist_filter);
                    if nFilt>0
                        fprintf('Eliminating: %i close candidates.  Now %i/%i Left\n',nFilt, sum(filter),length(filter));
                    end
                end
            end

            R.filter = filter;
            R.emitters = E(filter,:);
            R.theta = cellmap(@(k) obj.ResultsLocalizeEmitters.rawTheta{k}(:,filter(boxes.scaleIndexes{k})), 1:1);
            obj.ResultsFilterEmitters = R;

            obj.times.filterEmitters = toc;
            fprintf(' Time: %.3fs\n',obj.times.filterEmitters);
            obj.setPhase(6);
        end
        
        function trackEmitters(obj)
            obj.checkPhase(6);
 
            fprintf('Phase: Track Emitters\n');
            tic;

            P = obj.checkParamsTrack();
            E = obj.ResultsFilterEmitters.emitters;
            L = obj.makeLocalizations(E);
            roi = [obj.ROI(1:4), obj.ROI(7:8)];
            frameIdx = E(:,end);
            pos = E(:,1:2);
            SE_pos = E(:,9:10);
            features =  E(:,3);
            SE_features = E(:,11);
            obj.tracker = BaseTrack(roi, frameIdx, pos, SE_pos, features, SE_features);
            obj.tracker.D = P.D;
            obj.tracker.kon = P.Kon;
            obj.tracker.koff = P.Koff;
            obj.tracker.featureTol = P.LambdaR;
            obj.tracker.maxSpeed = P.MaxSpeed;
            obj.tracker.maxGapCloseFrames = P.MaxGapCloseFrames;
            obj.tracker.minGapCloseTrackLength = P.MinGapCloseTrackLength;
            obj.tracker.minFinalTrackLength = P.MinFinalTrackLength;

            obj.tracker.doLAP();

            R.tracks = obj.tracker.makeTracksArray(L);
            R.trackLengths = cellfun(@(t) size(t,1), R.tracks);

            %Sort the tracks
            [R.trackLengths, sidx] = sort(R.trackLengths,2,'descend');
            R.tracks = R.tracks(sidx);
            
%             dT = obj.data.frameT;
%             for ti=1:length(R.tracks)
%                 Ts=R.tracks{ti};
%                 times = (Ts(:,end)-1)*dT;
%                 R.tracks{ti} = [times, Ts];
%             end
            obj.ResultsTrack = R;

            obj.times.trackEmitters = toc;
            fprintf(' Time: %.3fs\n',obj.times.trackEmitters);
            obj.setPhase(7);
        end

        %% Data Retrieval Methods
        function L=makeLocalizations(obj, E)
            pS = obj.data.pixelSize;
            lambda = obj.data.lambdaPixelsToWavelength(E(:,3));
            Slambda = obj.data.lambdaPixelsToWavelength(E(:,3)-E(:,8))-lambda;
            SE_lambda = obj.data.lambdaPixelsToWavelength(E(:,3)-E(:,11))-lambda;
            SE_Slambda = obj.data.lambdaPixelsToWavelength(E(:,3)-E(:,16))-lambda;
            Fs = E(:,end)-1+obj.ROI(7);
%             ROIx=obj.ROI(1)-1;
%             ROIy=obj.ROI(3)-1;
            Ts = (Fs-1)*obj.data.frameT;
            L = [Ts, E(:,1:2)*pS, lambda, E(:, 4:5), E(:, 6:7)*pS, Slambda,...
                 E(:, 9:10)*pS, SE_lambda, E(:,12:13), E(:,14:15)*pS, SE_Slambda, Fs];
        end

        function [L,col_names] = getLocalizations(obj)
            % Localizations are sorted by frame
            % L - matrix format with column names given and each row a localization.  
            % col_names - Cell array of descriptions for each column of L.
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
            col_names = obj.LocalizationColumnNames;
        end

        %% Data access
        function RGB=getFramesRGB(obj)
            obj.checkPhase(2);
            RGB=HSData.makeRGB(obj.getFrames(), obj.ROIColorMap);
        end


        %% Plotting and Visualization
        % All methods make a new figure if called from without axes
        function plotSumImage(obj,axH)
            obj.checkPhase(2);
            if nargin==1
                figure('Units','pixels','Position',[10 10 obj.DefaultWindowSize]);
                axH = axes('Units','pixels');
            end
            BaseRPT.plot_imagesc(axH,obj.sumImage,obj.ROIOrigin)
            if nargin==1
                obj.positionImageAxes(axH, obj.data.physicalSize, [0 0 obj.DefaultWindowSize],[10 10 10 10]);
            end
        end

        function plotSpectra(obj,axH)
            obj.checkPhase(2);
            if nargin==1
                f=figure('Units','pixels','Position',[10 10 obj.DefaultWindowSize]);
                axH = axes('Units','pixels');
                whitebg(f);
            end
            obj.data.plotSpectralSeries(axH,obj.ROI);
           
        end

        function f=viewFramesRGB(obj)
            RGB = obj.getFramesRGB();
            dipframes=joinchannels('RGB',RGB);
            f = obj.viewMaximizedDipFig(dipframes);
            f.Name = 'RGB Image';
        end

        function f=viewFilteredFramesRGB(obj,varargin)
            obj.checkPhase(2);
            frames=obj.getFilteredFrames(varargin{:});
            RGB=HSData.makeRGB(frames, obj.ROIColorMap );
            f = obj.viewMaximizedDipFig(joinchannels('RGB',RGB));
            f.Name = 'RGB Filtered Image';
        end


        function plotRawMaximaImageRGB(obj, axH)
            obj.checkPhase(3);
            if nargin==1
                f = figure('Units','pixels','Position',[10 10 obj.DefaultWindowSize]);
                axH = axes('Units','pixels');
            end
            im = obj.ResultsFindMaxima.rawMaximaImage;
            RGB = HSData.makeRGB(im);
            BaseRPT.plot_imagesc(axH,RGB,obj.ROIOrigin)
            if nargin==1
                obj.positionImageAxes(axH, obj.data.physicalSize, [0 0 obj.DefaultWindowSize],[10 10 10 10]);
                f.Name = 'Raw Maxima RGB';
            end
        end

        function plotMaximaImageRGB(obj, axH)
            obj.checkPhase(4);
            if nargin==1
                f = figure('Units','pixels','Position',[10 10 obj.DefaultWindowSize]);
                axH = axes('Units','pixels');
            end
            im = obj.ResultsFilterMaxima.maximaImage;
            RGB = HSData.makeRGB(im);
            BaseRPT.plot_imagesc(axH,RGB,obj.ROIOrigin)
            if nargin==1
                obj.positionImageAxes(axH, obj.data.physicalSize, [0 0 obj.DefaultWindowSize],[10 10 10 10]);
                f.Name = 'Raw Maxima RGB';
            end
        end

        function f = viewBoxesMovieRGB(obj)
            obj.checkPhase(4);
            obj.updateWaitbar(0.1,'Generating Boxes 2D');
            im = single(obj.boxxer.plotRGBBoxCoordsDIP(obj.getFrames(), obj.getBoxCoords(), obj.ROIColorMap));
            im = joinchannels('RGB',im(:,:,:,1), im(:,:,:,2),im(:,:,:,3));
            f = obj.viewMaximizedDipFig(im);
            obj.updateWaitbar(1);
        end

       function plot3DTrackSequence(obj, opts)
            %  opts [optional] - If provided this must be 2rd argument, so you must provide trackIds also.  
            %     It is a strcut which allow for setting of several options.
            %  
            Ts = obj.getTracks();            
            Cs = cellmap(@(i) i*ones(size(Ts{i},1),1), 1:length(Ts));
            opts.cLabel = 'TrackID';
            opts.trackColorMap = @prism;
            opts.trackColorRange = length(Ts);
            RPT.plotTracks3D(Ts, Cs, obj.ROIPhysical, obj.sumImage, opts);
       end

       function plot3DTrackWavelength(obj, opts)
            %  opts [optional] - If provided this must be 2rd argument, so you must provide trackIds also.  
            %     It is a strcut which allow for setting of several options.
            %  
            Ts = obj.getTracks();            
            Cs = cellmap(@(T) T(:,4), Ts);
            opts.cLabel = 'Wavelength';
            opts.trackColorMap = @(~) flipud(obj.ROIColorMap);
            opts.trackColorRange = size(obj.ROIColorMap,1);
%             [im, imColormap] = rgb2ind(obj.sumImage, obj.ROIColorMap);
%             opts.imageColorMap = @(~) imColormap;
%             opts.imageColorRange =  size(imColormap,1);
            HSRPT.plotTracks3D(Ts, Cs, obj.ROIPhysical, obj.sumImage, opts);
        end 
       
        function fframes = getFilteredFrames(obj, sigma, method)
            obj.checkPhase(2); 
            frames=obj.getFrames();
            P = obj.checkParamsFindMaxima();
            if nargin<3
                method = P.method;
            end
            if nargin==1
                sigma=P.filterSigmas(:,1);
            elseif isscalar(sigma)
                sigma=P.filterSigmas(:,sigma);
            else
                sigma=single(sigma(:));
            end
            if isempty(obj.boxxer)
                obj.boxxer = Boxxer3D(obj.frameSize, P.filterSigma);
            end
            switch method
                case 'DoG'
                    fframes = obj.boxxer.filterDoG(frames, sigma, 1.1);
                case 'LoG'
                    fframes = obj.boxxer.filterDoG(frames, sigma);
                otherwise
                    error('HSRPT:ParamValue','Unknown filtering method %s', method);
            end
        end

        function hstm = makeTrackMovie(obj, trackIds)
            obj.checkPhase(7); 
            if nargin==1 || isempty(trackIds)
                trackIds = 1:obj.nTracks;
            end
            Ts = obj.getTracks(trackIds);
            min_track_length = 10;
            trackIds(min_track_length > cellfun(@(t) size(t,1),Ts(trackIds))) = []; %remove short tracks
            hstm = HSTrackMovie(obj, trackIds);
            hstm.setTrackColorMethod('Wavelength')
            hstm.viewSequence();
        end


    end % Public Methods
    
    %% Dependent properties
    methods
        function val = get.frameSize(obj)
            % Frame size [L Y X] in pixels
            if isempty(obj.ROI)
                val = [];
            else
                val = [obj.ROI(6)-obj.ROI(5)+1, obj.ROI(4)-obj.ROI(3)+1, obj.ROI(2)-obj.ROI(1)+1];
            end
        end
        
        function N = get.nFrames(obj)
            if isempty(obj.ROI)
                N = 0;
            else
                N = obj.ROI(8)-obj.ROI(7)+1;
            end
        end
        
        function origin = get.ROIOrigin(obj)
            % shift [X,Y,L] in pixels that should be added to local ROI pixel
            % coords to translate into global ROI coords
            if isempty(obj.ROI)
                origin = [];
            else
                origin = obj.ROI([1,3,5])-1;
            end
        end
        
        function pROI = get.ROIPhysical(obj)
            % The ROI in physical coordinates of um, nm and s instead of pixels
            % and frames [X Y L T]
            if isempty(obj.ROI)
                pROI = [];
            else
                pS = obj.data.pixelSize;
                dT = obj.data.frameT;
                pROI = [(obj.ROI(1)-1)*pS, obj.ROI(2)*pS,...
                        (obj.ROI(3)-1)*pS, obj.ROI(4)*pS,...
                        obj.data.lambdaPixelsToWavelength([obj.ROI(5),obj.ROI(6)+1]),...
                        (obj.ROI(7)-1)*dT, obj.ROI(8)*dT];
            end
        end

        function porigin = get.ROIPhysicalOrigin(obj)
            % shift [X,Y] in um that should be added to local physical
            % coords to translate into global physical coords
            if isempty(obj.ROI)
                porigin = [];
            else
                porigin = (obj.ROI([1,3])-1).*obj.data.pixelSize;
            end
        end
    
        function cm = get.ROIColorMap(obj)
            % Get an appropriate colormap for this ROI
            cm=obj.data.colorMap(obj.ROI(5):obj.ROI(6),:);
        end
    end %public methods
    
    methods (Static=true)
        function tableT = convertTracksTable(Ts)
            % This is static method to preform the reformatting of the HSRPT tracks format into a single table 
            % where all tracks have been concatenated and distinguished by an extra trackID column.  
            % This is slower an more inefficient than the default HSRPT Tracks format, but can be convenient as Matlab
            % tables have the ability to store the names, units, and descritptions of all columns.
            % [in]
            %   Ts - A cellarray of track matricies in normal HSRPT Tracks format
            % [out] 
            %  tableT - The tracks converted into a table giving each tracks information using an inital column for trackID
            tableT = array2table(cellmatfun(@(i) [repmat(i,size(Ts{i},1),1), Ts{i}]' ,1:length(Ts))');
            tableT.Properties.DimensionNames = HSRPT.TrackTableDimensionNames;
            tableT.Properties.Description = HSRPT.TrackTableTitle;
            %Add a Track index column to the table output.
            tableT.Properties.VariableNames = ['trackID',HSRPT.TrackColumnNames];
            tableT.Properties.VariableDescriptions = ['Track Index', HSRPT.TrackColumnDescriptions];
            tableT.Properties.VariableUnits = [{''}, HSRPT.TrackColumnUnits];
        end
        
        function structT = convertTracksStruct(Ts)
            % This is static method to preform the reformatting of the HSRPT tracks format into a structure array,
            % where each track is stored as a structure with field corresponding to the columns of the normal HSRPT
            % track format.
            % [in]
            %   Ts - A cellarray of track matricies in normal HSRPT Tracks format
            % [out] 
            %  tableT - The tracks converted into a structre array with one structure per track
            cols = cellmap(@(i) cellmap(@(t) t(:,i),Ts), 1:HSRPT.NTrackColumns);
            args = [HSRPT.TrackColumnNames; cols];
            structT = struct(args{:});
        end
        
        function trackHs = plotTracks3D(Ts, Cs, ROIphysical, im, opts)
            % This operates in physical units of micron(um) and seconds.  Draws tracks using the surface
            % command which is OpenGL accelerated.
            %
            % [IN]
            % Ts - Tracks in HSRPT Track format - cell array of matricies of localizations
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

            %Determine color range
            minCs = min(cellfun(@min,Cs));
            maxCs = max(minCs+1,max(cellfun(@max,Cs))); %assure span is at least 1
            spanCs = maxCs-minCs;

            %Plot Image
            im = im - min(im(:)); %shift to 0
            im = im ./ max(im(:)); %Scale to [0,1]
            im = (im.*spanCs.*(image_cRange/track_cRange))+maxCs*1.0001; %Shift image up to the top of the colorspace so it gets mapped to gray
            minT = min(cellfun(@(T) min(T(:,1)),Ts)); %Min time
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
            for i = length(Ts):-1:1
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
            axis('tight');

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
       end

        
         function checkTracks(tracks)
            %check the tracks are OK
            if ~iscell(tracks)
                error('HSRPT:checkTracks','Tracks must be a cell-array track format');
            elseif ~all(cellfun(@(t) size(t,2),tracks) == HSRPT.NTrackColumns)
                error('HSRPT:checkTracks','All tracks must have %i colums',HSRPT.NTrackColumns);
            end
        end
    end % Public static methods



    methods (Access=protected)
%         function L=getmakeLocalizations(obj, E)
%             pS = obj.data.pixelSize;
%             lambda = obj.data.lambdaPixelsToWavelength(E(3,:));
%             Slambda = obj.data.lambdaPixelsToWavelength(E(3,:)-E(8,:))-lambda;
%             SE_lambda = obj.data.lambdaPixelsToWavelength(E(3,:)-E(11,:))-lambda;
%             SE_Slambda = obj.data.lambdaPixelsToWavelength(E(3,:)-E(16,:))-lambda;
%             Fs = E(17,:)-1+obj.ROI(7);
%             ROIx=obj.ROI(1)-1;
%             ROIy=obj.ROI(2)-1;
%             L = [(E(1,:)+ROIx)*pS; (E(2,:)+ROIy)*pS; lambda; E(4:5,:); E(6:7,:)*pS; Slambda;...
%                  E(9:10,:)*pS; SE_lambda; E(12:13,:); E(14:15,:)*pS; SE_Slambda; Fs]';
%         end

        function loadData(obj, data, roi_in, roi_name)
            [obj.frames_, obj.ROI] = data.getFrames(roi_in);
            if isempty(obj.frames_)
                error('RPT:loadSPData','unable to load frames');
            end
            obj.frames_loaded = true;
            obj.data = data;
            obj.workingDir=obj.data.getFilePath('HSRPT');
            if isempty(obj.workingDir)
                obj.workingDir=fullfile(obj.data.workingDir,'HSRPT');
            end
            obj.Paths.data=relativepath(obj.workingDir, obj.data.saveFilePath);
            if isscalar(roi_in) && roi_in>0 && roi_in<=length(obj.data.ROI)  %Index into predefined ROI             
                obj.ROIname = obj.data.ROIname{roi_in};
                obj.Paths.saveFile = sprintf('%s_%s%s',obj.data.saveFileBaseName,obj.ROIname,obj.saveFileExt);
            elseif isempty(roi_in) %Full frame ROI
                obj.ROIname = 'FullFrame';
                obj.Paths.saveFile = sprintf('%s_%s%s',obj.data.saveFileBaseName,obj.ROIname,obj.saveFileExt);               
            else %Manual ROI
                obj.ROIname = 'manualROI';
                pattern = sprintf('%s_%s%%i%s',obj.data.saveFileBaseName,obj.ROIname,obj.saveFileExt);
                obj.Paths.saveFile = Pickle.findUnusedFileName(obj.workingDir,pattern);
            end
            obj.sumImage = obj.data.getSumImageRGB(obj.ROI);
            obj.initialized=true;
            obj.setPhase(2); %Initialized
        end
        
        function filterFrames(obj)
            % Actually do the filtering of the frames and set the properties
            obj.checkPhase(2); 
            P=obj.checkParamsFindMaxima();

            % Boxxer uses natural, [X Y] sized images.  We use flip to
            % convert to [Y X] sized images.
            if isempty(obj.boxxer)
                obj.boxxer=Boxxer3D(   );
            end
            switch P.method
                case 'LoG'
                    obj.filtered_frames_=obj.boxxer.filterLoG(obj.getFrames());
                case 'DoG'
                    obj.filtered_frames_=obj.boxxer.filterDoG(obj.getFrames());
                otherwise
                    error('HSRPT:filterFrames','Unknown filter method "%s"',P.method);
            end
            if ~isempty(obj.filtered_frames_)
                obj.filtered_frames_loaded=true;
            end
        end


        function im = computeMaximaImage(obj, maxima, maximaVals)
            mean_im = zeros(obj.frameSize);
            max_im = zeros(obj.frameSize);
            fs = obj.getFilteredFrames();
            for mi=1:length(maximaVals)
                Li = maxima(1, mi);
                yi = maxima(2, mi);
                xi = maxima(3, mi);
                si = maxima(4, mi);
                mean_im(Li,yi,xi) = mean_im(Li, yi,xi)+maximaVals(mi);
                max_im(Li,yi,xi) = max(max_im(Li, yi,xi),maximaVals(mi));
            end
            max_im = max_im./max(max_im(:));
            mean_im = mean_im/obj.nFrames;
            mean_im = mean_im./max(mean_im(:));
            im = 0.5*(mean_im+max_im);
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
