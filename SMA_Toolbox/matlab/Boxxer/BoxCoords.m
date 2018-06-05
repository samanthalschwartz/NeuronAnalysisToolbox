%
% Mark J. Olah (mjo@cs.unm.edu)
% 08/06/14
%
% This class is small class for holding and manipulating box coordinates  for a sequences
% of 2-D or 3-D images in the many formats necessary for visualization and plotting.  
% This removes the need to
% manipulate the box coordinates in opaque formats, and standardizes the storage
% and manipulation of coordinates

classdef BoxCoords < matlab.mixin.Copyable
    properties (SetAccess=protected)
        %Properties determined by the dataset type
        Ndim; % 2 or 3
        Nscales; % number of different scales
        Nframes; %number of frames
        imsize; %size of an individual frame as returned by size() function. 2D: [Y, X] = [#rows, #cols]
        scaleSigmas; %size:[Ndim,Nscales] where there are N different scales
        optimalBoxSize; %size:[Ndim,Nscales] giving optimal box size for each scale
        optimalBoxSigmaWidth=3; % if we need to determine optimal box size we use 2*optimalBoxSigmaWidth*psfSigma+1

        %properties that describe the boxes
        boxMaxima; % size:[Ndim x Npoints] locations of local maxima (points of interest) [row1=center row idx, row2=center col_idx]
        boxOrigin; % size:[Ndim x Npoints] -- upper left-hand corner of the box (the corner with smallest coordinates in all dims) [row1=origin row idx, row2=origin col_idx]
        boxSize; % size:[Ndim x Npoints] -- size [x y] for each box
        boxFrameIdx; %vector size:Npoints -- giving frame index for each box
        boxScaleIdx; %vector size:Npoints -- giving scale index for each box
        boxSizeCategoryIdx; %vector size:Npoints -- giving index into sizeCategories column corresponding to this size

        %Vital statistics about the box distributions
        Npoints; %number of points of interest
        NsizeCategories; %integer.  Number of different size categories.
        sizeCategories; % 2xN matrix giving the potential box sizes. rows=[X Y] or rows=[X Y L]
        
        %index lists to allow reverse-lookup of box indexes.
        frameIndexes; % cell array of vectors length:Nframes -- list of box indexes for each frame
        scaleIndexes; % cell array of vectors length:Nscales -- list of box indexes for each scale
        sizeIndexes; % cell array of vectors length:NsizeCategories -- list of box indexes for each size category 
    end

    properties (Access=protected, Transient=true)
        overlayMaskImage;
    end

    methods
        function obj=BoxCoords(maxima, imsize, scaleSigmas, Nframes, optimalBoxSize)
            obj.imsize = double(imsize(:)');
            obj.Ndim = numel(obj.imsize);
            if isscalar(scaleSigmas)
                scaleSigmas=[scaleSigmas;scaleSigmas];
            elseif size(scaleSigmas,1)==1 && size(scaleSigmas,2)==obj.Ndim
                scaleSigmas=scaleSigmas(:);
            elseif size(scaleSigmas,1)~=obj.Ndim;
                error('BoxCoords:BoxCoords','Bas scaleSigmas: %s', num2str(scaleSigmas));
            end
            obj.scaleSigmas = double(scaleSigmas);
            obj.Nscales = size(obj.scaleSigmas,2);
            if nargin<4
                obj.Nframes=max(maxima(end,:));
            else
                obj.Nframes=Nframes;
            end
            if nargin<5
                obj.optimalBoxSize = ceil(obj.scaleSigmas.*(obj.optimalBoxSigmaWidth*2)+1);
            else
                obj.optimalBoxSize = ceil(optimalBoxSize);
                if isvector(obj.optimalBoxSize) % Ensure if there is just one boxsize that it is a column vector
                    obj.optimalBoxSize = obj.optimalBoxSize(:);
                end
            end
            %Make sure the boxes are smaller than the image dimensions
            for d=1:obj.Ndim
                obj.optimalBoxSize(d,:) = min(obj.optimalBoxSize(d,:),double(obj.imsize(d)));
            end

            obj.createBoxes_Basic(maxima);
       end

        function [roi, findex]=makeROIsingle(obj, frames, idx)
            % Make a singleton ROI image for a given box index.
            sizeIdx = obj.boxSizeCategoryIdx(idx);
            roi = zeros(obj.sizeCategories(:,sizeIdx)','double');                
            o = obj.boxOrigin(:,idx);
            s = obj.boxSize(:,idx);
            e = o+s; %extent
            findex = obj.boxFrameIdx(idx);
            if obj.Ndim==2
                roi(:) = frames(o(1):e(1)-1,o(2):e(2)-1,findex); 
            else
                roi(:) = frames(o(1):e(1)-1,o(2):e(2)-1,o(3):e(3)-1,findex); 
            end
            roi(:) = max(0,roi(:)); %remove negative values
        end


        function [roi, findex]=makeROI(obj, frames)
            roi = cell(1,obj.NsizeCategories);
            findex = cell(1,obj.NsizeCategories);
            for i = 1:obj.NsizeCategories
                idxs = obj.sizeIndexes{i};
                N = length(idxs);
                ims = zeros([obj.sizeCategories(:,i)', N]);                
                for ni=1:N
                    n=idxs(ni);
                    o = obj.boxOrigin(:,n);
                    s = obj.boxSize(:,n);
                    e = o+s; %extent
                    f = obj.boxFrameIdx(n);
                    if obj.Ndim==2
                        ims(:,:,ni)=frames(o(1):e(1)-1,o(2):e(2)-1,f); 
                    elseif obj.Ndim==3
                        ims(:,:,:,ni)=frames(o(1):e(1)-1,o(2):e(2)-1,o(3):e(3)-1,f); 
                    end
                end
                ims(:) = max(0,ims(:)); %remove negative values
                roi{i} = ims;
                findex{i} = obj.boxFrameIdx(idxs);
            end
            if obj.NsizeCategories==1
                roi = roi{1};
                findex = findex{1};
            end
        end

        function imcube = makeROIcube(obj, frames)
            %make the rois in a cube as big as the biggest box size in each dimension where boxes
            %smaller than box size are inserted in the upper-left of the frame and the rest of the unused
            %pixels are 0
            maxsize = max(obj.sizeCategories,[],2);
            imcube = zeros([maxsize' obj.Npoints]);
            for n =1:obj.Npoints;
                o = obj.boxOrigin(:,n);
                s = obj.boxSize(:,n);
                e = o+s; %extent
                f = obj.boxFrameIdx(n);
                if obj.Ndim==2
                    imcube(1:s(1),1:s(2),n)=frames(o(1):e(1)-1,o(2):e(2)-1,f); 
                else
                    imcube(1:s(1),1:s(2),1:s(3),n)=frames(o(1):e(1)-1,o(2):e(2)-1,o(3):e(3)-1,f); 
                end
            end
            imcube(:) = max(0,imcube(:)); %remove negative values
        end
    
        function absEmitterTheta = covertToEmitterAbsoluteLocations(obj, estEmitterTheta)
            % [in]
            %   estEmitterTheta - Nparams x Npoints - Theta where first two rows are x and y poisition within the box
            % [out]
            %   absEmitterTheta - Nparams x Npoints - Theta where the x & y positions have been converted
            %                                       into aboslute coordinates
            shift = obj.boxOrigin([2,1],:)-1; %shift switches x/y since boxxer deals in row/col and locs are in x/y.  Offset is 1 pixel so the upper left of image is (0,0)
            absEmitterTheta = estEmitterTheta;
            absEmitterTheta(1:2,:) = estEmitterTheta(1:2,:) + shift;
        end

        function relEmitterPos = covertToEmitterBoxRelativeLocations(obj, absEmitterPos)
            % [in]
            %   absEmitterPos - Nparams x Npoints - Theta where the x & y positions have been converted
            %                                       into aboslute coordinates
            % [out]
            %   relEmitterPos - Nparams x Npoints - Theta where first two rows are x and y poisition within the box
            shift = obj.boxOrigin([2,1],:)-1; %shift switches x/y since boxxer deals in row/col and locs are in x/y.  Offset is 1 pixel so the upper left of image is (0,0)
            relEmitterPos = absEmitterPos - shift;
        end

       function relEmitterPos = covertToEmitterBoxCenterRelativeLocations(obj, absEmitterPos)
            % [in]
            %   absEmitterPos - Nparams x Npoints - Theta where the x & y positions have been converted
            %                                       into aboslute coordinates
            % [out]
            %   relEmitterPos - Nparams x Npoints - Theta where first two rows are x and y poisition within the box
            shift = (obj.boxOrigin([2,1],:)-1 + obj.boxSize([2,1],:)./2)'; %shift switches x/y since boxxer deals in row/col and locs are in x/y.  Offset is 1 pixel so the upper left of image is (0,0)
            relEmitterPos = absEmitterPos - shift;
        end


        function mask=overlayMask(obj,filter)
            %Number of frames to make mask
            sz=[obj.imsize, obj.Nframes];
            if nargin==1 
                if ~isempty(obj.overlayMaskImage) && ndims(obj.overlayMaskImage)==numel(sz) && all(size(obj.overlayMaskImage)==sz)
                    mask=obj.overlayMaskImage; %return cached copy
                    return
                else
                    filter=true(1,obj.Npoints);
                end
            end
            mask=zeros(sz,'single');
            for n=1:obj.Npoints
                if filter(n)
                    o = obj.boxOrigin(:,n);
                    s = obj.boxSize(:,n);
                    e = o+s-1; %extent
                    f = obj.boxFrameIdx(n);
                    if obj.Ndim==2
                        mask(o(1):e(1),o(2),f)=max(mask(o(1):e(1),o(2),f),1); 
                        mask(o(1):e(1),e(2),f)=max(mask(o(1):e(1),e(2),f),1); 
                        mask(o(1),o(2):e(2),f)=max(mask(o(1),o(2):e(2),f),1); 
                        mask(e(1),o(2):e(2),f)=max(mask(e(1),o(2):e(2),f),1); 
                        mask(o(1)+1:e(1)-1,o(2)+1:e(2)-1,f)= max(mask(o(1)+1:e(1)-1,o(2)+1:e(2)-1,f), 0.2);
                    elseif obj.Ndim==3
                        mask(o(1):e(1),o(2):e(2),o(3):e(3),f)=1; 
                        mask(o(1)+1:e(1)-1,o(2)+1:e(2)-1,o(3)+1:e(3)-1,f)=0.2;
                    end
                end
            end
            if nargin==1
                 obj.overlayMaskImage=mask; %only remember for the unfiltered case
            end
        end

        function overlayim = plotBoxes(obj, frames, filter, filterColorLabels)
            maxima_idx=sub2ind(size(frames),obj.boxMaxima(1,:), obj.boxMaxima(2,:), obj.boxFrameIdx)-1;
            mask=128*obj.overlayMask();
            overlayim=colorspace(dip_image(cosmicNorm(frames)*255),'RGB');
            overlayim{1}(maxima_idx) = 0;
            overlayim{2}(maxima_idx) = 0;
            overlayim{3}(maxima_idx) = overlayim{3}(maxima_idx)+128;
            if nargin>=3
                nColors = max(filter)+1;
                colors = [[0,1,0]; prism( ceil((nColors-1)/2)); cool( floor((nColors-1)/2))];
                for c=1:nColors
                    colorFilter = filter==c-1;
                    if ~any(colorFilter); continue; end
                    fmask=128*obj.overlayMask(colorFilter);
                    overlayim{1} = overlayim{1}+fmask*colors(c,1);
                    overlayim{2} = overlayim{2}+fmask*colors(c,2);
                    overlayim{3} = overlayim{3}+fmask*colors(c,3);
                end
                if nargin==4
                    pos = [100,100,300,500];
                    figure('Position',pos);
                    colormap(colors);
                    H=colorbar();
                    axis('off');
                    H.Position = [0.725 0.025 0.25 0.95];
                    H.Direction = 'reverse';
                    H.YAxisLocation='left';
                    H.Label.String='Box Filter Color Index';
                    H.Label.FontSize = 14;
                    H.Ticks=linspace(1/(2*nColors),1-1/(2*nColors),nColors);
                    H.TickLabels = filterColorLabels;
                    H.TickDirection = 'out';
                end
            else
                overlayim{1}=overlayim{1}+mask;
            end
        end

         function organizeBoxes(obj)
            %refresh the index lists for frames scales and sizes
            obj.frameIndexes = cellmap(@(f) find(obj.boxFrameIdx==f), 1:obj.Nframes);
            obj.scaleIndexes = cellmap(@(s) find(obj.boxScaleIdx==s), 1:obj.Nscales);
            obj.sizeIndexes  = cellmap(@(n) find(obj.boxSizeCategoryIdx==n), 1:obj.NsizeCategories);
        end
    end%public methods

    methods (Static = true)
        function newcoords = filter(boxcoords, keep_idx)
            %Create a new BoxCoords object from an old BoxCoords object and a boolean filter specifiying
            % boxes to keep.
            newcoords=copy(boxcoords);
            newcoords.boxMaxima = newcoords.boxMaxima(:,keep_idx);
            newcoords.Npoints = size(newcoords.boxMaxima,2);
            newcoords.boxScaleIdx = newcoords.boxScaleIdx(:,keep_idx);
            newcoords.boxFrameIdx = newcoords.boxFrameIdx(:,keep_idx);
            newcoords.boxSizeCategoryIdx = newcoords.boxSizeCategoryIdx(:,keep_idx);
            newcoords.boxSize = newcoords.boxSize(:,keep_idx);
            newcoords.boxOrigin = newcoords.boxOrigin(:,keep_idx);            
            newcoords.organizeBoxes();
        end
    end
    methods (Access=private)
        function createBoxes_Basic(obj, maxima)
            obj.Npoints = size(maxima,2);
            obj.boxMaxima = double(maxima(1:obj.Ndim,:));
            if size(maxima,1)==obj.Ndim+1
                obj.boxScaleIdx = ones(1,obj.Npoints);
            else
                obj.boxScaleIdx = maxima(obj.Ndim+2,:);
            end
            obj.boxFrameIdx = maxima(obj.Ndim+1,:);
            obj.sizeCategories = obj.optimalBoxSize;
            obj.NsizeCategories = size(obj.sizeCategories,2);
            obj.boxSize = obj.optimalBoxSize(:,obj.boxScaleIdx);
            obj.boxSizeCategoryIdx = obj.boxScaleIdx;
            obj.boxOrigin = obj.boxMaxima - fix(obj.boxSize/2);
            %shift any boxes which are over the edge
            for d=1:obj.Ndim
                obj.boxOrigin(d,:) = min(max(obj.boxOrigin(d,:), 1), obj.imsize(d)-obj.boxSize(d,:)+1);
            end
            obj.organizeBoxes();
        end

        %To be completed
%         function cleanBoxes(obj)
%             maxOverlap = 1;
%             for f=1:obj.Nframes
%                 boxIds = obj.frameIndicies{f};
%                 Nboxes = numel(boxIds);
%                 overlap = zeros(Nboxes,Nboxes);
%                 for n = 1:Nboxes
%                     on = obj.boxOrigin(:,n);
%                     sn = obj.boxSize(:,n);
%                     
% 
%                     for m = 1:Nboxes
%                     e = o+s; %extent
%                         overlap = 
%                     end
%                 end
%             end
%         end

    end
end %classdef
