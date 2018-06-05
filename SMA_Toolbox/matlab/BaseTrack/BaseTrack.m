classdef BaseTrack < handle
    %BaseTrack Basic Coordinate Tracking class
    %   Basic Tracking Class that does exactly what SPT RC1 did on the
    %   connect fits portion
    %   Peter Relich April 2015, UNM
    %   Mark J. Olah (mjo@cs.unm.edu) 2015, UNM
    properties (Constant=true, Hidden=true)
        MAX_BIRTHS = 15000; %More than this many births and we start to have a memory problem now.  Must use the sparseC++ version for these cases
    end
    properties
        %Emitter data:
        % These 5 properties are input at construction and specify the emitter's position.
        % In general emitters are represented by rows, and dimensions or different features are
        % represented by columns.  If there are no additional features for use, features and SE_featrues
        % can be omitted and left empty and they will not be used.   On the other hand, when using
        % features, the featureTol member variable should also be set. 
        ROI; % [xmin,xmax,ymin,ymax,tmin,tmax] x & y should be in pixels.  t's should be in frames.
        frameIdx; %The index of each emitter size:[Nemitters,1]
        position; %The estimated (2D) position of the emitter.   size:[Nemitters,Ndims]
        SE_position; %The estimated standard error on the position estimation in each dimension (i.e., sqrt(var(estimation)) );
        feature; % The feature values for each localization (e.g., lambda)  size:[Nemitters,Nfeatures].
        SE_feature; % The estimated measurement Standard Error for feature values for each localization (e.g., sqrt(var(est_lambda)) )  size:[Nemitters,Nfeatures].

        % Track Parameters.  This set of properties can be freely set by the user to control behavior.
        D= 1; % Diffusion Coefficient px^2/frame
        featureTol; %Effective feature deviation tolerated as a gaussian sigma parameters size:[1,Nfeatures], Units should match those of the feature
        kon= 0.1; % Blink On Rate [1/frame]
        koff= 0.1; % Blink Off Rate [1/frame]
        maxSpeed = 2.5; %Maximum speed in px/frame
        maxGapCloseFrames = 25;  %Maximum frame span over which to consider connections
        minGapCloseTrackLength = 1; %Minimum length of tracks to allow into the gap-closing phase.  Increasing will emilinate any short tracks BEFORE gap closing.
        minFinalTrackLength = 1; %Minimum length of tracks to return after gap closing.  It effectively must be >=minGapCloseTrackLength.)


        %Frame information
        frameStart; % Start frame index
        frameEnd; % End frame index

        % Counts
        Nemitters; %Number of emitters
        Ndims; %Number of spatial dimensions (at present Ndims=2 is only valid value.  Easily extendable to 3 at some future point.)
        Nfeatures; %number of features (e.g., lambda is a feature, also potentialy intensity or apparent sigma, etc.)
        Ntracks=0; %Number of produced tracks
        Nframes; %Total number of frames from first frame to last frame
        
        %Computed parameters
        ROIarea; % area of the tracking Region of Interest
        rho; % Particle Density [#particles/px]
                
        % Intermediate data
        frameRef; % cell array size [1,Nframes].  For each frame holds a vector of row (emitter) indexs in that frame. 
        trackStartFrameIdx; % frame index for a track start size:[Ntracks,1]
        trackEndFrameIdx; % frame index for a track end size:[Ntracks,1]
        trackStartLocIdx; % localization(emitter) index for a track start size:[Ntracks,1]
        trackEndLocIdx; % localization(emitter) index for a track end size:[Ntracks,1]

        links; % particle indexing property (matrix: Ntracks x frames)
        TrackLinks; % final matrix of track links
        GapCost; %Final gap closing costs
    end
    
    methods
        %% Constructor
        % The inputs will have to be changed, more information will be
        % required...
        function obj = BaseTrack(ROI, frameIdx, position, SE_position, feature, SE_feature)
            obj.ROI = ROI;
            obj.Nemitters = size(position,1); %emitters are rows
            obj.Ndims = size(position,2); %dimensions are cols
            assert(obj.Ndims==2); %Don't handle 3D yet...
            assert(all(size(position)==size(SE_position)));
            obj.position = position;
            obj.SE_position = SE_position;
            if nargin==4
                obj.Nfeatures=0;
            else
                obj.Nfeatures = size(feature,2);
                obj.feature = feature;
                obj.SE_feature = SE_feature;
                assert(all(size(obj.feature)==size(obj.SE_feature)));
                assert(size(obj.feature,1)==obj.Nemitters);
                obj.featureTol = ones(1,obj.Nfeatures); %defualt featureTol, should be changed later.
            end
            
            obj.frameIdx = frameIdx(:);
            assert(size(obj.frameIdx,1) == obj.Nemitters);

            obj.ROIarea = (obj.ROI(4)-obj.ROI(3)+1) * (obj.ROI(2)-obj.ROI(1)+1); %area in px^2
            obj.frameStart = obj.ROI(5);
            obj.frameEnd = obj.ROI(6);
            obj.Nframes = obj.frameEnd-obj.frameStart+1;
            % associate the coordinates by order by frame
            obj.computeFrameRef();
            % calculate avg. particle density
            obj.rho = obj.computeParicleDensity();
        end
         
        %% highest level of tracking code
        function doLAP(obj) 
            % intermediate variables of importance
            % perform rescaling for variables of interet
            n = zeros(obj.frameEnd-obj.frameStart+1,1);
            m = nnz(obj.frameIdx == obj.frameStart); % first set of coords
            Ind{obj.frameEnd-obj.frameStart} = [];
            NTracks = m;
            % loop through all the frames and run the frame 2 frame analysis
%             tic;
            for tt = obj.frameStart:obj.frameEnd-1
%                 if mod(tt-obj.frameStart,500)==0
%                     fprintf('Frame2Frame: %i\n',tt);
%                 end
                % make a relative time variable
                trel = tt-obj.frameStart+1;
                n(trel) = m;
                m = nnz(obj.frameIdx == tt+1); % get coords for next frame
                costFrame = obj.Frame2Frame(trel,n,m);
                % solve the LAP on the matrix if the matrix is not empty
                % otherwise skip, everything is empty for those frames
                if ~isempty(costFrame)
                [links12, links21, cost] = lap(costFrame,0,0,0);
                % store associations in the properties
                Ind{trel} = [links12, links21, cost];
                NTracks = NTracks + sum(links21(1:m)>n(trel));
                end
            end
            n(trel+1) = m; % populate last frame info
%             fprintf('Frame2Frame Time [Nframes=%i]: %.3fs\n', obj.Nframes, toc);
            if NTracks>obj.MAX_BIRTHS
                error('BaseTrack:ProblemSize','Got %i births which is sadly too many to fit in memory.  Try the C++ tracker',NTracks);
            end
            % build the link matrix, required for gap closing
%             fprintf('Link Matrix Creation\n');
%             tic;
            obj.linkMatrix(Ind,NTracks,n) 
%             fprintf('LinkMatrix Time: %.3fs\n',toc);

            % generate gap closing matrix 
%             fprintf('Making GapClosing Matrix\n');
%             tic;
            costGap = obj.GapClose(NTracks);
%             fprintf('GapClose Time[NTracks:%i]: %.3fs\n',NTracks,toc);
            % solve the gap closing matrix
%             fprintf('Solving GapClosing Matrix\n');
%             tic;
            [links12, links21, cost] = lap(costGap,0,0,0);
%             fprintf('SolveGapClose Time[NTracks:%i]: %.3fs\n',NTracks,toc);
            % export track data
            obj.GapCost = cost;
            obj.getOutput(NTracks, links12, links21);
        end
        
        %% lower level methods follow in the lines below
        
        %% Associate matrix coordinates with each frame for LAP calls
        function computeFrameRef(obj)
            T = obj.frameEnd-obj.frameStart+1; % number of frames to span
            obj.frameRef = cellmap(@(i) find(obj.frameIdx == i+obj.frameStart-1),1:T);
        end
        
        %% frame 2 frame code
        function costFrame = Frame2Frame(obj,tt,n,m)
            % populate birth and death and connection matrices
            % added logic to deal with blank frames
            caseval = logical(n(tt))*2 + logical(m);           
            switch caseval
                case 0 % there are no points in any frame, blank cost matrix
                    costFrame = {};
                case 1 % there are only points in the second frame, birth only matrix
                    costFrame = diag(ones(1,m));
                case 2 % there are only points in the first frame, death only matrix
                    costFrame = diag(ones(1,n(tt)));
                case 3 % there are points in both frames to connect
                    bM = obj.frameBirth(m);
                    dM = obj.frameDeath(n(tt));
                    cM = obj.frameConnect(tt,n(tt),m);
                    % build junk LR matrix
                    jM = (cM>0)';
                    jM = jM*eps; % minimize impact of costs with a small number
                    % build output cost matrix
                    costFrame = [cM dM; bM jM];
            end
        end
        
        function bM = frameBirth(obj,msz)
            val = -log(obj.rho) - log(obj.kon);
            bM = diag(val*ones(msz,1));
        end
        
        function dM = frameDeath(obj,nsz)
            val = -log(obj.koff);
            dM = diag(val*ones(nsz,1));
        end
        
        function distM = frameConnect(obj,tt,nsz,msz)
            % get coordinate positions
            nind = obj.frameRef{tt};
            mind = obj.frameRef{tt+1};
            % hard-coding coordinates now, will go dynamic once the tests
            % are satisfactory.
            Xn = obj.position(nind,1);
            Yn = obj.position(nind,2);
            Xm = obj.position(mind,1);
            Ym = obj.position(mind,2);
            varXn = obj.SE_position(nind,1);
            varYn = obj.SE_position(nind,2);
            varXm = obj.SE_position(mind,1);
            varYm = obj.SE_position(mind,2);
            % pre-initialize separable cost matrix components
            distM = zeros(nsz,msz);
            Xmat = zeros(nsz,msz);
            Ymat = zeros(nsz,msz);
            log2pi = log(2*pi);

            if obj.Nfeatures >= 1 %e.g., for hyperspectral data Lambda is a feature
                for k=1:obj.Nfeatures
                    Fn = obj.feature(nind,k);
                    Fm = obj.feature(mind,k);
                    varFn = obj.SE_feature(nind,k);
                    varFm = obj.SE_feature(mind,k);
                    Fmat = zeros(nsz,msz);
                    for jj = 1:msz
                        dF = Fn-Fm(jj);
                        varF = obj.featureTol(k)^2 + varFn + varFm(jj);
                        Fmat(:,jj) = 0.5*( dF.^2./varF + log2pi + log(varF)); 
                    end
                    distM = distM + Fmat;
                end
            end
            for jj = 1:msz
                dx = Xn-Xm(jj);
                varX = 2*obj.D + varXn + varXm(jj);
                Xmat(:,jj) = 0.5*( dx.^2./varX + log2pi + log(varX)); 
                dy = Yn-Ym(jj);
                varY = 2*obj.D + varYn + varYm(jj);
                Ymat(:,jj) = 0.5*( dy.^2./varY + log2pi + log(varY)); 
            end
            distM = distM + Xmat + Ymat - log(1-obj.koff);
        end
        
        %% build the link matrix
        function linkMatrix(obj,Ind,births,n)
            T = length(Ind)+1; % length of frames spanned
            % initialize links matrix and start and end indices
            obj.links = zeros(births,T);
            obj.links(1:n(1),1) = 1:n(1);
            obj.trackStartFrameIdx = zeros(births,1);
            obj.trackEndFrameIdx = zeros(births,1);
            % initialize frame start at 1 and correct with the real frame
            % start time at the end of the function
            obj.trackStartFrameIdx(1:n(1)) = 1;            
            mt = n(1); % number of assigned tracks, grows with birth
            % build links matrix, and trackStartFrameIdx,trackEndFrameIdx vectors
            for tt = 1:T-1                
               % skip the loop if there are no links
               if isempty(Ind{tt}) || length(Ind{tt}) == 1
                   continue;
               end               
               % define intermediate variables
               links12 = Ind{tt}(:,1);
               links21 = Ind{tt}(:,2);               
               % death counter or localization linkage
               obj.deathorlink(n,links12,tt);
               % birth counter, linkage start
               mt = obj.birthorlink(n,links21,mt,tt);
            end
            % get the starts and ends relative to the frameStart
            obj.trackEndFrameIdx(obj.trackEndFrameIdx==0) = T;
            %Use the pre-computed track endpoints to speed up the gap-closing
            %localization (row) index for the start of each track into the obj.position matrix, et. al.
            obj.trackStartLocIdx = arrayfun(@(i) obj.frameRef{obj.trackStartFrameIdx(i)}( obj.links(i,obj.trackStartFrameIdx(i)) ), ...
                                        1:births); 
            %localization (row) index for the start of each track into the obj.position matrix, et. al.
            obj.trackEndLocIdx = arrayfun(@(i) obj.frameRef{obj.trackEndFrameIdx(i)}( obj.links(i,obj.trackEndFrameIdx(i)) ),...
                                        1:births); 
        end
        
        %% death counter or link
        function deathorlink(obj,n,links12,tt)
            kkind = find(obj.links(:,tt))';
            % determine if there is a death or if a particle is linked
            for kk = kkind
                % death count
                if links12(obj.links(kk,tt)) > n(tt+1)
                    obj.links(kk,tt+1) = 0; % death
                    obj.trackEndFrameIdx(kk) = tt;
%                     obj.trackEndLocIdx(kk) = tt;
                    
                    % link to next id
                else
                    obj.links(kk,tt+1) = links12(obj.links(kk,tt));
                end
            end
        end
        
        %% birth counter, start a link point
        function mt = birthorlink(obj,n,links21,mt,tt)
            for bb = 1:n(tt+1)
                % congratulations, a birth!
                if links21(bb) > n(tt)
                    mt=mt+1; % increase birth count
                    obj.links(mt,tt+1) = bb;
                    obj.trackStartFrameIdx(mt) = tt+1;
                end
            end            
        end
        
        %% gap closing code
        function costGap = GapClose(obj,NTracks)
            % build birth matrix
            bM = obj.gapBirth;
            % build death matrix
            dM = obj.gapDeath;
            % build conenction matrix
            cM = obj.gapConnect(NTracks);
            % build junk LR matrix
            jM = (cM>0)';
            jM = jM*eps; % minimize impact of costs with a small number
            % build output cost matrix
            % Note to self: I need to make the gap closing matrix a sparse
            % off of 3 vectors, value (v), row (r) and column (c).  For now
            % I'm using a full matrix to get things to work, but the cross
            % over will be essential!
            costGap = [cM dM; bM jM];          
        end
        
        % birth matrix for gap closing
        function bM = gapBirth(obj)
            stTmp = obj.trackStartFrameIdx-1; % put track starts at base 0
            bM = diag(-log(obj.rho)-log(obj.kon)-(stTmp)*log(1-obj.kon));
        end
        % death matrix for gap closing
        function dM = gapDeath(obj)
            % determine frame size
            enTmp = obj.Nframes - obj.trackEndFrameIdx; % figure out gap lengths from track ends
            dM = diag(eps-(enTmp)*log(1-obj.kon));
        end

        function cM = gapConnect(obj,NTracks)
            cM = zeros(NTracks); % initialize connection matrix
            % This is a double for loop for now
            log2pi = log(2*pi);
            logkon = -log(obj.kon);
            logkoff= -log(1-obj.kon);
            gapCloseMaxDistance = 20*obj.D;
            for mm=1:NTracks
                %joining end of track nn to start of track mm
%                 if mod(mm,2000)==0 || mm==NTracks
%                     fprintf('[GapClosing] Handling Track: %i/%i\n',mm,NTracks);
%                 end
                stid=obj.trackStartLocIdx(mm);

                %vectorized temporal cuttoffs
                dT = obj.trackStartFrameIdx(mm)-obj.trackEndFrameIdx;
                temporally_feasible = find(dT>0 & dT<obj.maxGapCloseFrames);
                Ntemporally_feasible = numel(temporally_feasible);
                if Ntemporally_feasible == 0
                    continue;
                end
                
                %vectorized spatial cuttoffs
                sqdisp = (repmat(obj.position(stid,:),Ntemporally_feasible,1) - obj.position(obj.trackEndLocIdx(temporally_feasible),:)).^2;
                spatially_feasible = all(sqdisp<=repmat(min(gapCloseMaxDistance, 2*obj.D*dT(temporally_feasible)),1,2), 2);
                sqdisp = sqdisp(spatially_feasible);
                feasible = temporally_feasible(spatially_feasible);
                if isempty(feasible)
                    continue;
                end
                for ii = 1:numel(feasible)
                    nn = feasible(ii);
                    enid=obj.trackEndLocIdx(nn); % get the actual matrix index
                    
                    %calculate the squared spatial distance between ends and starts
                    posVar = 2*obj.D*dT(nn) + obj.SE_position(stid,:).^2 + obj.SE_position(enid,:).^2;
                    %compute costMat according to diffusion
                    cost = 0.5*(sum(sqdisp(ii,:)./posVar) + obj.Ndims*(log2pi + sum(log(posVar))));
                    if obj.Nfeatures > 1
                        for k = 1:obj.Nfeatures
                            Fdisp = obj.feature(enid,k) - obj.feature(stid,k);
                            varF = obj.featureTol(k) + obj.SE_feature(stid,k)^2 +obj.SE_feature(enid,k)^2;                        
                            cost = cost + 0.5*(Fdisp^2/varF + log2pi + log(varF));
                        end
                    end
                    cM(nn,mm) = cost + logkon+ logkoff; %Hacked this line to change to make longer connections more costly
                end
            end           
        end

        %% output formatting code
        function getOutput(obj,Nsegments,links12,links21)
            %Nsegments - number of pre-gap closing segments
            % determine number of tracks
            obj.Ntracks = sum(links21(1:Nsegments)>Nsegments);
            % get final track info
            obj.getTrackLinks(Nsegments,links12,obj.Ntracks)
        end
        
        function getTrackLinks(obj,Nsegments,links12,Ntracks)
            % takes the gap closing links and builds a final association
            % matrix for the original input matrix (TrackLinks)
            
            obj.TrackLinks=zeros(Ntracks,obj.Nframes); % build the TrackLinks matrix           
            st_grouping = obj.groupShortTracks(Nsegments,links12,Ntracks);
            % build the full track link matrix prior to frame ref
            % associaiton
            c_links = zeros(Ntracks,obj.Nframes);
            for nn = 1:Ntracks
               for mm = st_grouping{nn}
                    c_links(nn,:) = c_links(nn,:)+obj.links(mm,:);                
               end
            end
            % assign coordinate locations of tracks to the final link matrix
            for tt = 1:obj.Nframes
                tempc_links = c_links(:,tt);
                link_logical = tempc_links>0; % so we index track links properly!
                tempc_links = tempc_links(tempc_links>0); % so we don't index 0's!                
                obj.TrackLinks(link_logical,tt) = obj.frameRef{tt}(tempc_links);                
            end           
        end
        
        function st_grouping = groupShortTracks(~,Nsegments,links12,Ntracks)
            % group the assosciations of short tracks into separate cell
            % arrays for compressing the intermediate link matrix into NtracksxT size
            st_grouping{Ntracks} = [];
            % group short tracks together
            count = 1; % start a counter of the grouping matrix
            remvec = (1:Nsegments)';
            id = remvec(1); % start with an id here
            for nn = 1:Nsegments
                if id > Nsegments
                    count = count+1;
                    id = remvec(1);
                end
                st_grouping{count} = [st_grouping{count} id];
                remvec(remvec == id) = []; % remove chosen id from vector
                id = links12(id);
            end
        end
        
        function Tracks=makeTracksArray(obj,L)
            Tracks{obj.Ntracks} = []; % initialize cell array
            % loop through cell array indices
            for nn = 1:obj.Ntracks
                % output format of cell array
                % 'x', 'y', 'lambda', 'I', 'std_x', 'std_y', 'std_lambda', 'std_I', 'frame'
                %find all links for track nn
                TrackLinksnn = obj.TrackLinks(nn,:);
                TrackLinksnn = TrackLinksnn(TrackLinksnn>0); % remove 0 indices (intermittent gaps)
                Tracks{nn} = L(TrackLinksnn,:);                
            end
            if obj.minFinalTrackLength>1
                Tracks(cellfun(@(t) size(t,1)<obj.minFinalTrackLength, Tracks))=[];
            end
        end
        
        %% calculate the average density of the coordinates
        function rho=computeParicleDensity(obj)
            count = histc(obj.frameIdx, obj.frameStart:obj.frameEnd);
            % remove frames without localizations -- they won't be tracked
            rho = mean(count(count>0))/obj.ROIarea; % estimate density counts/pixel
        end
    end
end

