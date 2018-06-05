% Boxxer.m - The base class for 2D and 3D scale space filtering and box-finiding
% 05/2015
% Mark J. Olah (mjo@cs.unm.edu)
%
% This uses the MexIface C++/Matlab class to implement the scale space filtering in parallel with OpenMP

classdef Boxxer < IfaceMixin
    properties (Constant=true)
        ValidMaximaNeighborhoodSizes=[3, 5];
    end

    properties (SetAccess=protected)
        %dimension of each frame
        dim;
        % shape (dimensions) of image unit:pixels.  This is the size of each image.
        imsize;  % equivalent to size(im) [nrows, ncols] or [nrows, ncols, nslices]
        % sigma size of gaussian blobs to filter for order matches order of dimensions in image
        sigma; %size: [ndim x nscales] dims correspond to rows, cols, slices.
        nScales; %number of scales
        datatype='single';
        datacaster=@single;
    end

    properties (Access = protected)
        initialized=false; %True once the object is correctly initialized
    end
 
    methods
        function obj = Boxxer( iface )
            obj=obj@IfaceMixin(iface);
        end

        function success=initialize(obj, imsize, sigma)
            % [in] imsize - equivalent to size(image).  E.g. [nrows, ncols] for 2D  or [nrows, ncols, nslices] for
            %               3D (HS) data.
            obj.dim=length(imsize);
            if isscalar(sigma)
                sigma=[sigma;sigma]; %Scalar sigma is assumed to be isotropic
            end
            if size(sigma,1)~=obj.dim
                if size(sigma,1)==1 && size(sigma,2)==obj.dim
                    sigma=sigma(:);
                else
                    error('Boxxer:initialize','sigma rows expected: %i, got:%i',obj.dim,size(sigma,1));
                end
            end
            if ~all(imsize>0);
                error('Boxxer:initialize','Got non-positive image size');
            end
            if ~all(sigma>0);
                error('Boxxer:initialize','Got non-positive sigma');
            end
            obj.imsize = int32(imsize(:));
            obj.sigma = obj.datacaster(sigma);
            obj.nScales = size(obj.sigma,2);
               
            obj.initialized=obj.openIface(obj.imsize, obj.sigma);
            success=obj.initialized;
        end

        function setDoGSigmaRatio(obj, sigma_ratio)
            if ~isscalar(sigma_ratio) || sigma_ratio<=1
                error('Boxxer:ParamValue','sigma_ratio should be >1');
            end
            obj.call('setDoGSigmaRatio',sigma_ratio);
        end

        function fimage=filterScaledLoG(obj, image)
            % fimage=obj.filterLoG(image)
            % Filter using a Laplacian of Gaussian filter to detect blobs of size (scale)
            % given by sigma.
            %  [in] image: a single stack of imsize shaped frames, last dimension is time
            %  [out] fimage:a single hyperstack of imsizeY x imsizeX x nScales x nFrames shaped filtered frames
            obj.checkImage(image);
            fimage=obj.call('filterScaledLoG', image);
        end

        function fimage=filterScaledDoG(obj, image)
            % fimage=obj.filterDoG(image, sigma_ratio)
            % Filter using a Difference of Gaussian approximation to LoG filter
            %  [in] image: a single stack of imsize shaped frames, last dimension is time
            %  [out] fimage:a single hyperstack of imsizeY x imsizeX x nScales x nFrames shaped filtered frames
            obj.checkImage(image);
            fimage=obj.call('filterScaledDoG', image);
        end

        function [maxima, max_vals] = scaleSpaceLoGMaxima(obj, image, neighborhoodSize, scaleNeighborhoodSize)
            % Filter using a seperable LoG implementation at multiple scales
            %  [in] image: a single stack of imsize shaped frames, last dimension is time
            %  [in] neighborhoodSize: The size of the neighborhood for local maxima finding (use: 3 or 5)
            %  [in] scaleNeighborhoodSize: The size of the neighborhood for maxima finding  over scales(use: 3 or 5)
            %  [out] maxima: 4xN matrix of maxima rows are [xpos, ypos, scale, frame].  Scales and frames are
            %                1-based indexes
            %  [out] max_vals: 1xN vector of maxima values at each local maxima found.
            obj.checkImage(image);
            if nargin<4
                scaleNeighborhoodSize=3;
            end
            if nargin<3
                neighborhoodSize=5;
            end
            
            [maxima, max_vals] = obj.call('scaleSpaceLoGMaxima', image, int32(neighborhoodSize),...
                                                                int32(scaleNeighborhoodSize));
            maxima=maxima+1; %correct for 0-based C++ coords to 1-based Matlab coords;
        end

        function [maxima, max_vals] = scaleSpaceDoGMaxima(obj, image, neighborhoodSize, scaleNeighborhoodSize)
            % Filter using a seperable DoG implementation at multiple scales
            %  [in] image: a single stack of imsize shaped frames, last dimension is time
            %  [in] neighborhoodSize: The size of the neighborhood for local maxima finding (use: 3 or 5)
            %  [in] scaleNeighborhoodSize: The size of the neighborhood for maxima finding  over scales(use: 3 or 5)
            %  [out] maxima: 4xN matrix of maxima rows are [xpos, ypos, scale, frame].  Scales and frames are
            %                1-based indexes
            %  [out] max_vals: 1xN vector of maxima values at each local maxima found.
            obj.checkImage(image);
            if nargin<4
                scaleNeighborhoodSize=3;
            end
            if nargin<3
                neighborhoodSize=5;
            end
            
            [maxima, max_vals] = obj.call('scaleSpaceDoGMaxima', image, int32(neighborhoodSize),...
                                                                int32(scaleNeighborhoodSize));
            maxima=maxima+1; %correct for 0-based C++ coords to 1-based Matlab coords;
        end

        function fimage=filterLoG(obj, image, sigma)
            % fimage=obj.filterLoG(image)
            % Filter using a Laplacian of Gaussian filter to detect blobs of size (scale)
            % given by sigma.
            %  [in] image: a single stack of imsize shaped frames, last dimension is time
            %  [in] sigma: [optional] [nDim x 1] vector of sigmas to use
            %  [out] fimage: a single stack of imsize shaped filtered frames
            obj.checkImage(image);
            if nargin==2
                sigma=obj.sigma(:,1);
            else
                sigma=single(sigma(:));
                if length(sigma)~=obj.dim
                    error('Boxxer:BadParameter','Sigma value must be size: %ix1',obj.dim);
                end
            end
            fimage=obj.callstatic(obj.ifaceHandle,'filterLoG', image, sigma);
        end
        
        function fimage=filterDoG(obj, image, sigma, sigmaRatio)
            % fimage=obj.filterDoG(image, sigma_ratio)
            % Filter using a Difference of Gaussian approximation to LoG filter
            % to look for blobs of size (scale) given by sigma, and using ratio sigma_I/sigma_E=sigma_ratio
            %  [in] image: a single stack of imsize shaped frames, last dimension is time
            %  [in] sigma: [optional] [nDim x 1] vector of sigmas to use
            %  [in] sigmaRatio: a scalar giving the ratio of sigmas in the DoG method:
            %             sigma_I/sigma_E (I=inhibitory, E=excititory)
            %             [Default=obj.DefaultDoGSigmaRatio].
            %  [out] fimage:a single stack of imsize shaped filtered frames
            obj.checkImage(image);
            if nargin<=3
                sigmaRatio=1.1;
            end
            if nargin==2
                sigma=obj.sigma(:,1);
            else
                sigma=single(sigma(:));
                if length(sigma)~=obj.dim
                    error('Boxxer:BadParameter','Sigma value must be size: %ix1',obj.dim);
                end
            end
            if ~isscalar(sigmaRatio)
                error('Boxxer:filterDoG','Unexpectedly got vector for sigmaRatio.');
            end
            if sigmaRatio<=1
                error('Boxxer:filterDoG','Sigma ratio must be >1.');
            end
            fimage=obj.callstatic(obj.ifaceHandle,'filterDoG', image, sigma, single(sigmaRatio));
        end

        function fimage=filterGauss(obj, image, sigma)
            % fimage=obj.filterGauss(image)
            % Filter using a Gaussian to smooth image
            %  [in] image: a single stack of imsize shaped frames, last dimension is time
            %  [in] sigma: [optional] [nDim x 1] vector of sigmas to use
            %  [out] fimage: a single stack of imsize shaped filtered frames
            obj.checkImage(image);
            if nargin==2
                sigma=obj.sigma(:,1);
            else
                sigma=single(sigma(:));
                if length(sigma)~=obj.dim
                    error('Boxxer:BadParameter','Sigma value must be size: %ix1',obj.dim);
                end
            end
            fimage=obj.callstatic(obj.ifaceHandle,'filterGauss', image, sigma);
        end

        function [maxima, max_vals]=enumerateImageMaxima(obj, image, neighborhoodSize )
            % maxima=obj.enumerateImageMaxima(image, neighborhoodSize)
            % Use the Non-maximum suppression algorithm to return all pixel coordinates that are local maxima of
            % their neighborhood.  The neighborhood size is the linear dimension of the neighborhood, so in 3D
            % a neighborhoodSize of 3 implies that the maxima is bigger than any of the 26 pixels that surround it.
            %
            % [in] image: a single stack of imsize shaped frames, last dimension is time
            % [in] neighborhoodSize: an odd integer. Will be converted to uint32.  Acceptable values are in ValidMaximaNeighborhoodSizes.  (default=3)
            % [out] maxima: a (dim+1)xN matrix of maxima where rows are X, Y, ..., T and columns are different maxima detected.
            % [out] max_vals; a Nx1 Vector giving the value at each maxima.
            if nargin==2
                neighborhoodSize=obj.ValidMaximaNeighborhoodSizes(1);
            end
            obj.checkImage(image);
            neighborhoodSize=int32(neighborhoodSize);
            if mod(neighborhoodSize,2)~=1 || neighborhoodSize<3
                error('Boxxer:ParmaValue','Got invalid neighborhoodSize: %i',neighborhoodSize);
            end
            [maxima, max_vals]=obj.callstatic(obj.ifaceHandle,'enumerateImageMaxima', image, neighborhoodSize);
            maxima=maxima+1; %correct for 0-based C++ coords to 1-based Matlab coords;
        end
      
        function box_coords = generateBoxCoords(obj, maxima, Nframes, optimalBoxSize)
            %
            % [out] box_coords -
            if nargin==4
                box_coords = BoxCoords(maxima, obj.imsize, obj.sigma, Nframes, optimalBoxSize);
            else
                box_coords = BoxCoords(maxima, obj.imsize, obj.sigma, Nframes);
            end
        end
        

    end %public methods

    methods (Access=protected)
        function checkImage(obj, im)
            if ~isa(im,'single')
                error('Boxxer:filterLoG','Expected type "%s" got "%s".', obj.datatype, class(im));
            end
            input_size=size(im);
            if ~all(input_size(1:obj.dim)==obj.imsize')
                error('Boxxer:filterLoG','Got incorrect sized image %s',mat2str(input_size));
            end
        end
    end

    methods (Access=public, Abstract=true)
        checkMaxima(obj, image, maxima, max_vals);
    end

    methods (Static=true)

        function boxsize = defaultBoxSize(sigma)
            boxsize = ceil(6*sigma+1);
        end

        function [fmaxima, fmax_vals, filter] = filterMaxima(maxima, max_vals, thresh_val)
            % [in] maxima - size [3,N] for non-scale space or [4,N] for scale space.  Postions of maxima
            %                          rows: [x_px, y_px, frame_idx, scale_idx]
            % [in] max_vals - size [N] list of maximum values a maxima positions
            % [in] thresh_val - [optional] A threshold value above which the maxima will be selected.
            %                              if ommitted then use triThres to estimate threshold values.
            % [out] fmaxima - size [3,K] or [4,K] The filtered list of maxima same number of rows as maxima.
            % [out] fmax_vals - size [K] The filtered list of maxima values corrsponding to the maxima
            % [out] filter - size [N].  1 if maxima was selected.  0 otherwise
            if nargin < 3
                [smax, sidx] = sort(max_vals,'descend');
                threshIdx = triThres(smax, false);
                fmaxima = maxima(:,sidx(1:threshIdx));
                fmax_vals = max_vals(sidx(1:threshIdx));
                filter = zeros(1,numel(max_vals));
                filter(sidx(1:threshIdx)) = 1;
            else
                filter = max_vals>=thresh_val;
                fmaxima = maxima(:,filter);
                fmax_vals = max_vals(filter);
            end
        end
        
         function [fmaxima, fmax_vals] = filterMaximaAuto(maxima, max_vals, minNum)
             % Return at least minNum if possible
             % [in]
             %  maxima - [nD x N] vector of indexs to maxima
             %  max_vals - [1 xN] vector of max_values associated with maxima
             %  minNum - Integer minum to return
             % [out]
             %  fmaxima - Filtered maxima
             %  fmax_vals - Filtered maxima_vals
            [smax, sidx] = sort(max_vals,'descend');
            threshIdx = triThres(smax, false);
            if nargin==3 && minNum>0 && threshIdx<minNum
                threshIdx=min(minNum,numel(max_vals));
            end
            fmaxima = maxima(:,sidx(1:threshIdx));
            fmax_vals = max_vals(sidx(1:threshIdx));
        end


        function overlayIm = plotScaleMaximaRPT(im, maxima, selected)
            % This works in the RPT format where the order of maxima is [x y frame scale]
            % [IN]
            %  im - A grayscale image
            %  maxima - 4xN coordinate list of maxima
            %  maximaVals - length N vector of filter values at maxima
            %  selected - [optional] length N vector of booleans 1 if maxima is selected 0 otherwise
            % [OUT]
            %  overlayIm - dip_image RGB image
            imScales = size(im,size(maxima,1)); % number of scales present in image
            overlayIm = colorspace(dip_image(cosmicNorm(im)*255),'RGB');
            maximaScales = maxima(4,:);
            for n=1:imScales
                index_cols=num2cell([maxima(1:3,:); n*ones(1,size(maxima,2))],2); %cell array of rows of maxima
                atScale = (maximaScales == n)';
                inds = sub2ind(size(im),index_cols{:})-1; %minus 1 for dipimage coords
                normVals=255; %Added this to highlight maxima better
                if nargin<4
                    overlayIm{1}(inds) = normVals;
                    overlayIm{2}(inds) = 0;
                    overlayIm{3}(inds) = normVals .* ~atScale;
                else
                    overlayIm{1}(inds) = normVals .* ~selected;
                    overlayIm{2}(inds) = normVals .* selected;
                    overlayIm{3}(inds) = normVals .* (selected & ~atScale);
                end
            end
        end

        function overlayIm = plotScaleMaxima(im, maxima, selected)
            % This works in the Boxxer format where the order of maxima is [x y scale frame]
            % [IN]
            %  im - A grayscale image size 4-D [x y scale frames]
            %  maxima - 4xN coordinate list of maxima
            %  maximaVals - length N vector of filter values at maxima
            %  selected - [optional] length N vector of booleans 1 if maxima is selected 0 otherwise
            % [OUT]
            %  overlayIm - dip_image RGB image
            imScales = size(im,3); % number of scales present in image
            imFrames = size(im,4);
            if imFrames==1
                maxima=maxima(1:3,:);
            end
            overlayIm = colorspace(dip_image(cosmicNorm(im)*255),'RGB');
            index_cols=num2cell(maxima,2); %cell array of rows of maxima
            inds = sub2ind(size(im),index_cols{:})-1; %minus 1 for dipimage coords
            normVals=255; %Added this to highlight maxima better
            if nargin<3
                overlayIm{1}(inds) = normVals;
                overlayIm{2}(inds) = 0;
                overlayIm{3}(inds) = 0;
            else
                overlayIm{1}(inds) = normVals .* ~selected;
                overlayIm{2}(inds) = normVals .* selected;
                overlayIm{3}(inds) = 0;
            end
        end

        
        function overlayIm = plotMaxima(im, maxima, selected)
            % [IN]
            %  im - A grayscale image
            %  maxima - 3xN coordinate list of maxima
            %  maximaVals - length N vector of filter values at maxima
            %  selected - [optional] length N vector of booleans 1 if maxima is selected 0 otherwise
            % [OUT]
            %  overlayIm - dip_image RGB image
            overlayIm = colorspace(dip_image(im.*(255/max(im(:)))),'RGB');
            index_cols = num2cell(maxima,2); %cell array of rows of maxima
            inds = sub2ind(size(im),index_cols{:})-1; %minus 1 for dipimage coords
            normVals=255; %Added this to highlight maxima better
            if nargin<3
                overlayIm{1}(inds) = normVals;
                overlayIm{2}(inds) = 0;
                overlayIm{3}(inds) = 0;
            else
                overlayIm{1}(inds) = normVals .* ~selected;
                overlayIm{2}(inds) = normVals .* selected;
                overlayIm{3}(inds) = 0;
            end
        end



        function overlayim = plotTrueBoxCoordsDIP(image, box_coords, true_coords)
            overlayim=Booxer.plotBoxCoordsDIP(image, box_coords);
            true_coords_mask=zeros(size(image),'uint16');
            true_coords_mask(sub2ind(size(image),true_coords(:,2),true_coords(:,1),true_coords(:,3)))=1;
            overlayim{3}=overlayim{3}+255*true_coords_mask; %draw true coords in blue
        end
        
        function overlayim = plotBoxCoordsDIP(im, box_coords, maxima, maximaVals, filter)
            maxima_index_cols=num2cell(maxima(1:3,:),2); %cell array of rows of maxima
            maxima_idx=sub2ind(size(im),maxima_index_cols{:})-1;
            if nargin>=4
                normVals = cosmicNorm(maximaVals)*255;
            else
                normVals=255;
            end    
            mask=128*box_coords.overlayMask();
            overlayim = colorspace(dip_image(cosmicNorm(im)*255),'RGB');
            overlayim{1}(maxima_idx)=0;
            overlayim{2}(maxima_idx)=0;
            overlayim{3}(maxima_idx)=normVals;
            overlayim{1}=overlayim{1}+mask;
            if nargin==5
                fmask=128*box_coords.overlayMask(filter);
                overlayim{2}=overlayim{2}+fmask;
            end
        end

        function [im, cm]=plotBoxCoordsOverlay(image, box_coords)
            index_cols=num2cell([box_coords.boxCenter; box_coords.boxFrameIdx],2); %cell array of rows of maxima
            inds=sub2ind(size(image),index_cols{:});            
            mask=box_coords.overlayMask();
            im=cosmicNorm(image);
            cm=[gray(64); gray(64)+repmat([1 0 0],64,1);gray(64)+repmat([0.5 0 0],64,1);gray(64)+repmat([1 1 0],64,1)  ];
            cm(:)=min(cm(:),1);
            im(mask==1)=mod(im(mask==1),1)+1;
            im(mask==0.5)=mod(im(mask==0.5),1)+2;
            im(inds)=mod(im(inds),1)+3;
        end
        
        function im=plotEmitterMovie(box_coords, box_images)
            im=zeros([box_coords.imsize, box_coords.Nframes]);
            for szi=1:box_coords.NsizeCategories            
                ims = cosmicNorm(box_images{szi});
                
                idxs = box_coords.sizeIndexes{szi};
                for ni=1:length(idxs)
                    n = idxs(ni);
                    o = box_coords.boxOrigin(:,n);
                    s = box_coords.boxSize(:,n);
                    f = box_coords.boxFrameIdx(n);
                    e = o+s-1;
                    if box_coords.Ndim==2
                        im(o(1):e(1),o(2):e(2),f) = max(im(o(1):e(1),o(2):e(2),f), ims(:,:,ni)); 
                    else
                        im(o(1):e(1),o(2):e(2),o(3):e(3),f) = im(o(1):e(1),o(2):e(2),o(3):e(3),f) + ims(:,:,:,ni); 
                    end
                end
            end
        end

        function overlayim=plotRGBBoxCoordsDIP(image, box_coords, cm)
            %This is for hyperspectral data.  Try to combine this with the
            %2D plotting
            rgbimage=makeRGB(image, cm);
            maxima_mask=zeros(size(image),'uint16');
            index_cols=num2cell([box_coords.boxCenter; box_coords.boxFrameIdx],2); %cell array of rows of maxima
            inds=sub2ind(size(image),index_cols{:});
            maxima_mask(inds)=1;
            maxima_mask=squeeze(sum(maxima_mask,1));
            mask=box_coords.overlayMask();
            mask=squeeze(sum(mask,1));
            overlayim=rgbimage;
            overlayim{1}=overlayim{1}+(mask./max(mask(:)));
            overlayim{1}=overlayim{1}+(maxima_mask./max(maxima_mask(:)));
            overlayim{2}=overlayim{2}+(maxima_mask./max(maxima_mask(:)));
            overlayim{3}=overlayim{3}+(maxima_mask./max(maxima_mask(:)));
        end
    end
    
end %classdef
