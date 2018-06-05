% File: SRRender2D.m
%
% Mark J. Olah (mjo@cs.unm.edu)
% 03/2015
%
% A class for super-resolution rendering of gaussians and histrograms.  
% We provide methods to render histograms and gaussians as a single 2D image or as a movie sequence of 
% 2D images.
%
% This class is a C++/Matlab hybrid using the MexIface class interface.  All computation is done in C++ in
% parallel using openMP.

classdef SRRender2D < IfaceMixin
    %
    % * Units: The spatial units used are arbitrary and can be pixels, microns, nm, or whatever.  Just
    %          make sure the same "world" units are used consitently thoughout.
    % * Points Format:  The format of the input points is a matrix where each row is a localization and
    %                   the columns are [I, x, y, sigma_x, sigma_y, frameIdx].  frameIdx is optional and
    %                   is only needed for generating "movie" outputs.  The intensity units should be all
    %                   ones if you don't want to distinguish the points by intensity, but can also
    %                   represent any other feature which will be mapped to effectrive inensity of
    %                   output.
    % * Image Format: The resulting image always has size [sizeY sizeX] where the increasing rows
    % correspond to increasing Y-coordinates, and increasing columns correspond to increasing column
    % coordinates.  This is the default coordiante system for DipImaga and matlab's image toolbox.  Of
    % note, the coordiantes of the edge of the image will exactly match the ROI, except because they
    % must be integer coordinates the resulting coordinates may need to be rounded up to the next pixel
    % boundary.
    %  (see: http://www.mathworks.com/help/images/ref/imref2d-class.html)
            
    properties
        sigmaAccuracy = 5; %The accuracy with which gaussians will be rendered in multiples of sigma
    end
    properties (SetAccess=protected)
        ROI; % [xmin xmax ymin ymax] size in units matching the emitters that will be passed in
        datatype; %This should be 'single' or 'double'
    end

    properties (Access = protected)
        datacaster; %Handle to method to convert data to this type, @single or @double matchin obj.datatype
    end
    
    methods
        function obj=SRRender2D(roi, type)
            % Construct a new SRRender object or a blank object.
            %
            % Usage:
            %   obj=SRRender(size, type);
            %
            %   [in] ROI=[xmin, xmax, ymin, ymax] - 4-element vector giving effective region of interest that
            %                                       will act as the default ROI for the images generated.
            %   [in] (optional) type = 'single' or 'double'.  Single will use less memory.  
            %                   Important: all inputs to methods should be of this datatype.
            %                   [default='single']
            %
            if nargin==1
                type = 'single';
            end
            switch type
                case 'double'
                    iface=@SRRender2DDouble_Iface;
                    datatype='double';
                    datacaster=@double;
                case 'single'
                    iface=@SRRender2DSingle_Iface;
                    datatype='single';
                    datacaster=@single;
                otherwise
                    error('SRRender2D:constructor','Unsupported type: %s', class(a));
            end
            obj = obj@IfaceMixin(iface);
            obj.datatype=datatype;
            obj.datacaster=datacaster;
            obj.ROI=obj.checkROI(roi);
        end
        
        function [im, imCoords] = renderHist(obj, points, imSizePx, roi)
            % Generates a single histogram image of the points weigheted by intensity.  The points columns
            % for sigma_x and sigma_y are irrelevent for this method.
            %
            % [in] points - matrix of points in standard format
            % [in] imSizePx - (scalar) The image size of resulting image along the longest coordinate of roi.
            % [in] roi - [optional] The roi to image for with format [xmin, xmax, ymin, ymax].
            %            [Default: use the class's obj.ROI property]
            % [out] im - A hisrogram image with maximum dimension given by imSizePx.
            % [out] imCoords - An imref2d object describing the mapping from pixels to real world coordinates in
            %                  the resulting image. 
            if nargin<4
                roi = obj.ROI;
            else
                roi = obj.checkROI(roi);
            end
            points = obj.checkPoints(points);
            [im, imCoords] = obj.makeImage(imSizePx, roi);
            effectiveROI = obj.datacaster([imCoords.XWorldLimits, imCoords.YWorldLimits]);
            obj.callstatic(obj.ifaceHandle,'renderHist', points, effectiveROI, im);
        end

        function [im, imCoords] = renderGauss(obj, points, imSizePx, roi)
            % Generates a single frame of a gaussian blob representation of the given points.
            %
            % Note: to control width of rendered gaussians in sigma multiples set obj.sigmaAccuracy
            %
            % [in] points - matrix of points in standard format
            % [in] imSizePx - (scalar) The image size of resulting image along the longest coordinate of roi.
            % [in] roi - [optional] The roi to image for with format [xmin, xmax, ymin, ymax].
            %            [Default: use the class's obj.ROI property]]
            % [out] im - A hisrogram image with maximum dimension given by imSizePx.
            % [out] imCoords - An imref2d object describing the mapping from pixels to real world coordinates in
            %                  the resulting image. 
            if nargin<4
                roi = obj.ROI;
            else
                roi = obj.checkROI(roi);
            end
            points = obj.checkPoints(points);
            [im, imCoords] = obj.makeImage(imSizePx, roi);
            effectiveROI = obj.datacaster([imCoords.XWorldLimits, imCoords.YWorldLimits]);
            obj.callstatic(obj.ifaceHandle,'renderGauss', points, effectiveROI, obj.sigmaAccuracy, im);
        end

        function [im, imCoords] = renderHistMovie(obj, points, imSizePx, roi)
            % Generates a sequence of histogram images of the points weigheted by intensity.  The points columns
            % for sigma_x and sigma_y are irrelevent for this method, but must still be included.  The
            % column format is [I x y sigma_x sigma_y frameIdx].  frameIdx starts at 1 and must be an
            % integer.
            %
            % [in] points - matrix of points in standard format with last column (col 6) giving frame
            %               index.  Frame indexes should start at 1 and be integer valued.  They
            %               correspond directly to the resulting image sequence returned.
            % [in] imSizePx - (scalar) The image size of resulting image along the longest coordinate of roi.
            % [in] roi - [optional] The roi to image for with format [xmin, xmax, ymin, ymax].
            %            [Default: use the class's obj.ROI property]
            % [out] im - A hisrogram image with maximum dimension given by imSizePx.
            % [out] imCoords - An imref2d object describing the mapping from pixels to real world coordinates in
            %                  the resulting image. 
            if nargin<4
                roi = obj.ROI;
            else
                roi = obj.checkROI(roi);
            end
            points = obj.checkPointsMovie(points);
            nFrames = max(points(:,6));
            [im, imCoords] = obj.makeImage(imSizePx, roi, nFrames);
            effectiveROI = obj.datacaster([imCoords.XWorldLimits, imCoords.YWorldLimits]);
            points(:,6) = points(:,6) -1; %Convert to 0-based indexing for C++
            obj.callstatic(obj.ifaceHandle,'renderHistMovie', points, effectiveROI, im);
        end

        function [im, imCoords] = renderGaussMovie(obj, points, imSizePx, roi)
            % Generates a sequence of gaussian blob images of the points weigheted by intensity.  The
            % column format is [I x y sigma_x sigma_y frameIdx].  frameIdx starts at 1 and must be an
            % integer.
            %
            % [in] points - matrix of points in standard format with last column (col 6) giving frame
            %               index.  Frame indexes should start at 1 and be integer valued.  They
            %               correspond directly to the resulting image sequence returned.
            % [in] imSizePx - (scalar) The image size of resulting image along the longest coordinate of roi.
            % [in] roi - [optional] The roi to image for with format [xmin, xmax, ymin, ymax].
            %            [Default: use the class's obj.ROI property]
            % [out] im - A hisrogram image with maximum dimension given by imSizePx.
            % [out] imCoords - An imref2d object describing the mapping from pixels to real world coordinates in
            %                  the resulting image.
            if nargin<4
                roi = obj.ROI;
            else
                roi = obj.checkROI(roi);
            end
            points = obj.checkPointsMovie(points);
            nFrames = max(points(:,6));
            [im, imCoords] = obj.makeImage(imSizePx, roi, nFrames);
            effectiveROI = obj.datacaster([imCoords.XWorldLimits, imCoords.YWorldLimits]);
            points(:,6) = points(:,6) -1; %Convert to 0-based indexing for C++
            obj.callstatic(obj.ifaceHandle,'renderGaussMovie', points, effectiveROI, im);
        end
    end %public methods

    methods (Access=protected)
        function points = checkPoints(obj,points)
            if any( points(:,1)<=0 )
                error('SRRender2D:checkPoints','Bad intensity values');
            end
            points = obj.datacaster(points);
        end

        function points = checkPointsMovie(obj,points)
            points = obj.datacaster(points);
            if size(points,2)<6
                error('SRRender2D:checkPointsMovie','Not enough columns.  Got size:%s',num2str(size(points)));
            end
            fs=points(:,6);
            if any(fs<=0) || any(~isfinite(fs)) || any(round(fs)~=fs)
                error('SRRender2D:checkPointsMovie','Bad frame indexes in column 6');
            end
        end

        function [im, imCoords] = makeImage(obj, imSizePx, roi, nFrames)
            % Make a blank image or video sequence of the correct size and type
            if nargin==3
                nFrames = 1;
            end
            xdist = roi(2)-roi(1);
            ydist = roi(4)-roi(3);
            if xdist>=ydist
                xSize = imSizePx;
                ySize = ceil(imSizePx*ydist/xdist);
                xWorldLimits = [roi(1), roi(2)];
                yWorldLimits = [roi(3), roi(3)+ySize*(xdist/xSize)];      
            else
                xSize = ceil(imSizePx*xdist/ydist);
                ySize = imSizePx;
                xWorldLimits = [roi(1), roi(1)+xSize*(ydist/ySize)];
                yWorldLimits = [roi(3), roi(4)];
            end
            im = zeros(ySize,xSize,nFrames,obj.datatype);
            imCoords = imref2d(size(im),xWorldLimits,yWorldLimits);
        end
        
        function roi=checkROI(obj,roi)
            %check the formatting or roi is correct
            roi = obj.datacaster( roi(:)');
            if any(~isfinite(roi)) || length(roi)~=4 || roi(1)>=roi(2) || roi(3)>=roi(4) 
                error('SRRender2D:checkROI','Invalid roi: %s',num2str(roi));
            end
        end

    end %protected methods

    methods (Static=true)
        function points = simulateData(roi, nPoints, meanSigma)
            intensities = rand(nPoints,1)+1;
            pos = SRRender2D.makeTestGrid(roi,nPoints);
            sigmas = gamrnd(4,meanSigma/4,nPoints,2);
            points = single([intensities, pos, sigmas]);
        end

        function points = simulateMovieData(roi, nPoints, nFrames, meanSigma)
            points=zeros(nPoints*nFrames, 6, 'single');
            for n = 1:nFrames
                intensities = rand(nPoints,1)+1;
                pos = SRRender2D.makeTestGrid(roi,nPoints);
                sigmas = gamrnd(4,meanSigma/4,nPoints,2);
                frameIdxs = n*ones(nPoints,1);
                points((n-1)*nPoints+(1:nPoints),:) = [intensities, pos, sigmas, frameIdxs];
            end
        end

        function demo()
            roi = [243, 536, 10, 158]; %[xmin xmax ymin ymax]%
            nFrames = 50;
            nPoints = 10000;
            meanSigma = 0.3;
            imSizePx = 12048; %number of picels along longest edge in resulting images
            movieSizePx = 4096; %number of picels along longest edge in resulting images
    
            %Make a new SRRender object
            srr = SRRender2D(roi);     
            %Simualate points wich include a frameIdx column suitable for movie plots
            tic;
            points = SRRender2D.simulateMovieData(roi, nPoints, nFrames, meanSigma);
            tPoints = size(points,1);
            fprintf('Simulate N=%i points. Time:%.5fs\n',tPoints, toc);
            %Single histogram image
            tic;
            [histim, histImRef] = srr.renderHist(points,imSizePx);
            elapsed = toc;
            fprintf('Render Histogram Image [N=%i] Time:%.5fs Emitters/s:%.4g\n', tPoints, elapsed, tPoints/elapsed);
            
            %Single gaussian image
            tic;
            [gaussim, gaussImRef] = srr.renderGauss(points,imSizePx);
            elapsed = toc;
            fprintf('Render Gaussian Image [N=%i] Time:%.5fs Emitters/s:%.4g\n', tPoints, elapsed, tPoints/elapsed);

            %Histogram movie
            tic;
            hist_movie = srr.renderHistMovie(points,movieSizePx);
            elapsed = toc;
            fprintf('Render Histogram Movie [N=%i] Time:%.5fs Emitters/s:%.4g\n', tPoints, elapsed, tPoints/elapsed);

            %Gaussian movie
            tic;
            gauss_movie = srr.renderGaussMovie(points,movieSizePx);
            elapsed = toc;
            fprintf('Render Gaussian Movie [N=%i] Time:%.5fs Emitters/s:%.4g\n', tPoints, elapsed, tPoints/elapsed);
            srr.viewDipImage(histim);
            srr.viewDipImage(gaussim);
            srr.viewDipImage(hist_movie);
            srr.viewDipImage(gauss_movie);
            histImRef
        end

        function renderSRTest(srtest, imSizePx)
            % This will render data retrieved from an SRTest object.
            %
            % [in] srtest - An SRTest object
            % [in] scale - integer giving maginfication scale factor.  1=actual pixel size.  10=10x resolution on
            %              both X and Y scales
            % [in] intensityMode - 
            % [out] im - The super-res image
            roi = [0, srtest.Dim1Size, 0, srtest.Dim2Size];           
            srr = SRRender2D(roi,'single');
            min_frame=min(srtest.Results.FrameNum);
            movieSize = 2048;
            frameCompressFactor = 10; %How many frames of emitters to group together
            frameIdx = ceil((srtest.Results.FrameNum-min_frame+1)/frameCompressFactor);
            nPoints = length(srtest.Results.X);

            %Change intensities
            %intensities = ones(nPoints, 1); %uniform intensity
            intensities =srtest.Results.Photons; %scale intensity by photon count
            
            points = [intensities srtest.Results.X srtest.Results.Y srtest.Results.X_STD srtest.Results.Y_STD,  frameIdx];
            tPoints = size(points,1);
            %Single histogram image
            tic;
            [histim, histImRef] = srr.renderHist(points,imSizePx);
            elapsed = toc;
            fprintf('Render Histogram Image [N=%i] Time:%.5fs Emitters/s:%.4g\n', tPoints, elapsed, tPoints/elapsed);
            
            %Single gaussian image
            tic;
            [gaussim, gaussImRef] = srr.renderGauss(points,imSizePx);
            elapsed = toc;
            fprintf('Render Gaussian Image [N=%i] Time:%.5fs Emitters/s:%.4g\n', tPoints, elapsed, tPoints/elapsed);

            %Zoomed single gaussian image
            tic;
            roi = [3/8 5/8 3/8 5/8]*256; % center of image area
            [zoomGaussim, zoomGaussImRef] = srr.renderGauss(points,imSizePx,roi);
            elapsed = toc;
            fprintf('Render Zoomed Gaussian Image [N=%i] Time:%.5fs Emitters/s:%.4g\n', tPoints, elapsed, tPoints/elapsed);

            
            %Histogram movie
            tic;
            hist_movie = srr.renderHistMovie(points,movieSize);
            elapsed = toc;
            fprintf('Render Histogram Movie [N=%i] Time:%.5fs Emitters/s:%.4g\n', tPoints, elapsed, tPoints/elapsed);

            %Gaussian movie
            tic;
            gauss_movie = srr.renderGaussMovie(points,movieSize);
            elapsed = toc;
            fprintf('Render Gaussian Movie [N=%i] Time:%.5fs Emitters/s:%.4g\n', tPoints, elapsed, tPoints/elapsed);
            srr.viewDipImage(histim);
            srr.viewDipImage(gaussim);
            srr.viewDipImage(zoomGaussim);
            srr.viewDipImage(hist_movie);
            srr.viewDipImage(gauss_movie);
        end
        
        function fig = viewDipImage(image, fig)
            if nargin==1
                fig=figure();
            end
            if ~isa(image,'dip_image')
                image = dip_image(image);
            end
            dipshow(fig,image);
            diptruesize(100*900./max(size(image)))
            dipmapping(fig,'colormap',hot);
        end
    end %public static methods

    methods (Access=protected, Static=true)

        function ps = makeTestGrid(roi, nPoints)
            xs = rand(nPoints,1)*(roi(2)-roi(1))+roi(1);
            ys = rand(nPoints,1)*(roi(4)-roi(3))+roi(3);
            nGrids=5;
            noise = 0.025;
            gridSpacing = min( roi([4,2])-roi([3,1]))/nGrids;
            psum = xs+ys;
            delta_sum = (round(psum/gridSpacing)+noise*randn(nPoints,1))*gridSpacing - psum;
            pdiff = xs-ys;
            delta_diff = (round(pdiff/gridSpacing)+noise*randn(nPoints,1))*gridSpacing - pdiff;
            shiftSumMask = abs(delta_sum)<abs(delta_diff);
            shiftDiffMask = ~shiftSumMask;
            xs(shiftSumMask)=xs(shiftSumMask)+delta_sum(shiftSumMask)/2;
            ys(shiftSumMask)=ys(shiftSumMask)+delta_sum(shiftSumMask)/2;
            xs(shiftDiffMask)=xs(shiftDiffMask)+delta_diff(shiftDiffMask)/2;
            ys(shiftDiffMask)=ys(shiftDiffMask)-delta_diff(shiftDiffMask)/2;
            xs = max(xs,roi(1));
            xs = min(xs,roi(2));
            ys = max(ys,roi(3));
            ys = min(ys,roi(4));
            ps=[xs ys];
        end
    end %protected static methods
end %classdef
