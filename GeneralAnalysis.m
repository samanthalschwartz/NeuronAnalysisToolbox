classdef GeneralAnalysis < handle
properties
    test = [];
end

methods (Static)
    function img_out = cropimage(img_in,cutoff)
        % creates an output image that is a subset of the input image
        % inputs:
        %   img_in - dipimage or matrix image (converts to type dipimage)
        %   cutoff - structure array with fields
        %    cutoff.xrange - can be single number represented number of outer edge pixels to remove or the vector of pixels to include in dimension 1
        %    cutoff.yrange - can be single number represented number of outer edge pixels to remove or the vector of pixels to include in dimension 2
        %    cutoff.trange - can be single number represented number of outer edge pixels to remove or the vector of pixels to include in dimension 3. Use
        %               dipimage numbering so first frame is index 0!
        %   range values must match the dimensionality of the image
        %       example1: to crop 15 pixels from the edges of an image but
        %       keep all frames
        %           xrange = 15:size(img_in,1)-15;
        %           cutoff.yrange = 15:size(img_in,1)-15;
        %           cutoff.trange = 0:size(img_in,3)-15;
        %       example2: to crop 15 pixels from the edges along the x-axis
        %       and only include the first half of the frames
        %           xrange = 15
        %           yrange = [];
        %           trange = 0:(floor(size(img_in,3)/2)-1))
        % output:
        %   img_out - cropped dipimage object. use single(dipimage) to convert to matlab matrix.
        if ~isa(img_in,'dip_image')
            try
                img_in = dip_image(img_in);
            catch
                warning('input must be an image matrix');
                    return;
            end
        end
        if isempty(cutoff)
            img_out = img_in;
            return;
        end
        if isfield(cutoff,'xrange') && ~isempty(cutoff.xrange)
            if numel(cutoff.xrange)==1
                img_out = img_in(cutoff.xrange:(size(img_in,1)-cutoff.xrange-1),:,:);
            else
            img_out = img_in(cutoff.xrange,:,:);
            end
        end
        if isfield(cutoff,'yrange') && ~isempty(cutoff.yrange)
            if numel(cutoff.yrange)==1
                img_out = img_out(:,cutoff.yrange:(size(img_out,2)-cutoff.yrange-1),:);
            else
            img_out = img_out(:,cutoff.yrange,:);
            end
        end
        if isfield(cutoff,'trange') && ~isempty(cutoff.trange)
            if numel(cutoff.trange)==1
                img_out = img_in(:,:,cutoff.trange:size(img_in,3)-cutoff.trange-1,:,:);
            else
            img_out = img_in(:,:,cutoff.trange);
            end
        end
    end
    function img = old_loadtiff(filepath) %to delete
        % this function loads a tiff file into matlab and generates a dipimage
        % must have the bioformats function bfopen: download at https://docs.openmicroscopy.org/bio-formats/5.7.0/developers/matlab-dev.html
        imgbefore = bfopen(filepath);
        img = dip_image(zeros([size(imgbefore{1,1}{1}),size(imgbefore{1,1},1)]));
        for ii = 1:size(img,3)
            img(:,:,ii-1) = imgbefore{1,1}{ii};
        end
    end 
    function ov = displaytiff(image)
        if ndims(image) > 3
            switch size(image,4)
                case 1
                    ov = joinchannels('rgb',image(:,:,:,1));
                case 2
                    ov = joinchannels('rgb',image(:,:,:,1),image(:,:,:,2));
                case 3
                    ov = joinchannels('rgb',image(:,:,:,1)+image(:,:,:,3),image(:,:,:,2)+image(:,:,:,3),image(:,:,:,3));
            end
        else
            ov = joinchannels('rgb',image);
        end
    end
    function im_array = loadtiff_3ch(filepath)
        % requires loadtiff function from % Copyright (c) 2012, YoonOh Tak
        oimg = loadtiff(filepath);
        frames3 = size(oimg,3);
        ch1 = oimg(:,:,1:3:frames3);
        ch2 = oimg(:,:,2:3:frames3);
        ch3 = oimg(:,:,3:3:frames3);
        im_array = cat(4,ch1,ch2,ch3);
    end
    function im_array = loadtiff_2ch(filepath)
        % requires loadtiff function from % Copyright (c) 2012, YoonOh Tak
        oimg = loadtiff(filepath);
        frames2 = size(oimg,3);
        ch1 = oimg(:,:,1:2:frames2);
        ch2 = oimg(:,:,2:2:frames2);
        im_array = cat(4,ch1,ch2);
    end
    function ch = loadtiff_1ch(filepath)
        % requires loadtiff function from % Copyright (c) 2012, YoonOh Tak
        ch = dip_image(loadtiff(filepath));
    end
    function im_array = splitANDsavetiff_3ch(filepath)
       im_array = loadtiff_3ch(filepath);
       for i =1:3
          image = dip_image(im_array(:,:,:,i));
          [fpath,name,~] = fileparts(filepath);
          save(fullfile(fpath,[name '_ch' num2str(i)]),'image');
       end
    end
    function img_gauss = imgGauss(img_in,gsig)
        img_gauss = gaussf(img_in,gsig);        
    end
    function img_lapl = imgLaplace(img_in,lsig,gsig)
        % gaussian filters image and then calculates laplacian
        % inputs: 
            % img_in - dipimage or matrix image (converts to type dipimage)
            % lsig - kernal for laplacian. must be
            %           same dimension as img_in. example: [1 1 0] is  
            %           transfrom in x and y but not time.
            % gsigma - optional input to set the gaussian kernal. must be
            %           same dimension as img_in. example: [1 1 0] is gaussian 
            %           smoothing in x and y but not time.
        % outputs:
            % img_lapl - filtered dipimage. to convert to matlab array use
            %           single(img_out).
        if nargin<3
            gsig = ones(1,numel(size(img_in)));
            if nargin<2
                lsig = ones(1,numel(size(img_in)));
            end
        end
        img_g = gaussf(img_in,gsig);
        img_lapl = dxx(img_g,lsig)+dyy(img_g,lsig);
    end
    function img_laplcutoff = imgLaplaceCutoff(img_in,lsig,gsig)
        if nargin<3
            gsig = ones(1,numel(size(img_in)));
            if nargin<2
                lsig = ones(1,numel(size(img_in)));
            end
        end
        img_lapl = GeneralAnalysis.imgLaplace(img_in,lsig,gsig);
        img_laplcutoff = -img_lapl;
        img_laplcutoff(img_laplcutoff<0) = 0;
    end
    function img_dcc = imgDcc(img_in,gsig)
        if nargin<2
            gsig = 1;
        end
        img_dcc = dcc(img_in,gsig);
    end
    function img_dcccutoff = imgDccCutoff(img_in,gsig) %good for edge detection!
        if nargin<2
            gsig = 1;
        end
        img_dcc = GeneralAnalysis.imgDcc(img_in,gsig);
        img_dcccutoff = -img_dcc;
        img_dcccutoff(img_dcccutoff<0) = 0;
    end
        
    function mask = imgThreshold(img_in)
        threshval = multithresh(single(img_in),2);
        mask = img_in>=threshval(1);
    end
    function [mask,threshval] = imgThreshold_fixedUserInput(img_in)
        uiwait(msgbox('Select a representative background region','Title','modal'));
        h = dipshow(img_in,'log');
        diptruesize(h,125);
        [a,b] = dipcrop(h);
        threshval = max(a);
        mask = threshold(img_in,'fixed',threshval);
        close(h);
    end
    function newmask = cleanUpMask_manual(underimgin,mask_in)
        %        lb = label(mask_in);
%         ov = overlay(underimgin,mask_in);
%         h = dipshow(ov,'log');
%         dipmapping(h,'global')
        %       h = dipshow(lb,'labels');
%         while(ishandle(h))
%             [a b] = dipcrop(h);
%             mask_in(b(1,1):b(1,1)+b(2,1),b(1,2):b(1,2)+b(2,2),:) = 0;
%             close(h);
%             ov = overlay(underimgin,mask_in);
%             h = dipshow(ov,'log');
%             dipmapping(h,'global');
%         end
        lb = label(mask_in);
        ov = underimgin;
        ov(lb~=0) = 0;
        g = dipfig('ov');
        dipshow(ov,'log');
        diptruesize(g,200);
        while(ishandle(g))
            try
                v = dipgetcoords(g,1);
            catch
                break;
            end
            val = single(lb(v(1),v(2),v(3)));
            lb(lb == val) = 0; 
            ov = underimgin;
            ov(lb~=0) = 0
        end
      dipfig -unlink
      newmask = logical(lb);  
    end
    function perim = maskperim(mask_in)
         perim = dt(mask_in);
         perim = (perim==1);     
    end
    function [labeledim] = labelmask_byframe(mask_in,conn,minSize,maxSize)
        % this function labels a binary image/mask and assigns the label
        % value as the frame number
        % inputs:
        %   mask_in - 3D dipimage
        %   conn - connectivity for dipimage label function (type help label
        %               to get more info)
        % outputs:
        %   labeledim - label dipimage object
        if nargin<4
            maxSize = 0; %default no cutoff
            if nargin<3
                minSize = 0; %default no min
                if nargin<2
                    conn = 1;
                end
            end
        end
        labeledim = dip_image(zeros(size(mask_in)));
        if numel(size(mask_in))<3
            disp('This method is not useful for 2D images, because it will just return the mask. Try with a 3D image series.')
            return;
        end
        for ii = 1:size(mask_in,3)
            mask = dip_image(mask_in);
            temp = GeneralAnalysis.labelmask(mask(:,:,ii-1),conn,minSize,maxSize);
            temp(temp>0) = ii;
            labeledim(:,:,ii-1) = temp;
        end
    end
    function [labeledim] = labelmask(mask_in,conn,minSize,maxSize)
        % labels a binary image/mask and assigns the label
        % value as the frame number
        % inputs:
        %   mask_in - 3D dipimage
        %   conn - connectivity for dipimage label function (type help label
        %               to get more info)
        % outputs:
        %   labeledim - label dipimage object
        if nargin<4
            maxSize = 0; %default no cutoff
            if nargin<3
                minSize = 0; %default no min
                if nargin<2
                    conn = 1;
                end
            end
        end
        labeledim = label(mask_in,conn,minSize,maxSize);
    end
    function [labeledim] = labelmask_unique(mask_in,conn,minSize,maxSize)
        % this function labels a binary image/mask without connecting
        % labeled objects between fraame
        % inputs:
        %   mask_in - 3D dipimage
        %   conn - connectivity for dipimage label function (type help label
        %               to get more info)
        % outputs:
        %   labeledim - label dipimage object
        if nargin<4
            maxSize = 0; %default no cutoff
            if nargin<3
                minSize = 0; %default no min
                if nargin<2
                    conn = 1;
                end
            end
        end
        if numel(size(mask_in))<3 %image is just 2D
            [labeledim] = GeneralAnalysis.labelmask(mask_in,conn,minSize,maxSize);
            return;
        end
        labeledim = mask_in*0;
        for ii = 0:(size(mask_in,3)-1)
            labeledim(:,:,ii) = GeneralAnalysis.labelmask(mask_in(:,:,ii),conn,minSize,maxSize);
        end
    end
    function mask_out = bwmorph_timeseries(mask_in,fun_string,n_repeats)
        % applies Matlab bwmorph function to a each frame of a time series
        % type help bwmorph for all the great options!
        % examples include:
        %   'fill','bridge','close','branchpoints','endpoints','skel','thicken'
        if nargin<3
            n_repeats = Inf;
        end
        mask_out = dip_image(zeros(size(mask_in,2),size(mask_in,1),size(mask_in,3)));
        for ii = 1:size(mask_in,3)
            bwmframe = bwmorph(single(mask_in(:,:,ii-1)),fun_string,n_repeats);
            mask_out(:,:,ii-1) = bwmframe;
        end
    end
    function bgim_out = regionfill_timeseries(image_in,mask_in)
       % image_in: image with holes (0 values) where values are
       % interpolated from
       % mask_in: mask over which the values need to be interpolated 
       maxframe = size(image_in,3);
       
       if size(mask_in,3)==1
           mask_in = repmat(mask_in,1,1,maxframe);
       end
       im = zeros(size(single(image_in)));
       for tt = 1:maxframe
          I = single(image_in(:,:,tt-1));
          w = single(mask_in(:,:,tt-1));
          im(:,:,tt) = regionfill(I,w);
       end
        bgim_out = dip_image(im);
    end
    function wshed = watershed_timeseries(image_in,conn)
        wshed = dip_image(zeros(size(image_in,2),size(image_in,1),size(image_in,3)));
        for ii = 1:size(image_in,3)
            wshedframe = watershed(single(image_in(:,:,ii-1)),conn);
            wshed(:,:,ii-1) = wshedframe;
        end  
    end
    
    function mask_thick = thicken(mask_in,numpix)
        sumproj_out = GeneralAnalysis.sumproj_masktimeseries(mask_in);
        sumproj_out_thick = GeneralAnalysis.bwmorph_timeseries(sumproj_out,'thicken',numpix);
        mask_thick = GeneralAnalysis.bwmorph_timeseries(sumproj_out_thick,'bridge');
    end
    function sumproj_out = sumproj_masktimeseries(mask_in)
       sm = sum(mask_in,[],3);
       sm(sm>0) = 1;
       sumproj_out = repmat(sm,[1 1 size(mask_in,3)]);
    end
    function lbl_out = findLabelsInMask(lbl_in,mask)
        % Excludes labels in a labeled image that are exclusively out of the bounds of an input mask.
        % or Includes labels in a labeled image if any part of the label is within the bounds of an input mask.
        % inputs:
        %   lbl_in - labeled 3D dipimage (integer label values are pixel value for all labeled objects in the image) 
        %   mask_in - 3D dipimage
        % outputs:
        %   lbl_out - label dipimage object. label numbering is preserved
        if numel(size(lbl_in))<3
            tsize = 1;
        else
            tsize= size(lbl_in,3);
        end
        lbl_out = lbl_in*0;
        for tt = 0:(tsize-1)
            lblframe = lbl_in(:,:,tt);
            maskframe = mask(:,:,tt);
            test = lblframe*maskframe;
            lbid = unique(single(test));
            ids2remove = find(~ismember(1:max(lblframe),lbid));
            lbl_outframe = lblframe;
            for ii = ids2remove
                lbl_outframe(lblframe == ii)=0;
            end
            lbl_out(:,:,tt) = lbl_outframe;
        end
    end
    function distMat = geodesic_seedDistfromMask(sink_mask,seed_mask,geom_mask,plotflag,plotsavedir,saveflag)
        if nargin<6
            saveflag = 0;
        end
        if nargin<5
            plotsavedir = pwd;
        elseif nargin<4
            plotflag = 0;
            plotsavedir = [];
        end
        assert(isequal(size(sink_mask),size(seed_mask)) & isequal(size(sink_mask),size(geom_mask)));
        if numel(size(sink_mask))<3
            tsize = 1;
        else
            tsize= size(sink_mask,3);
        end
        wb = waitbar(0,'Analyzing Distances (this may take a while...)');
        if plotflag
            pathfig = figure();
            pathax = gca;
        end
        % make a guess at distance matrix by using # of objects at a
        % random frame * # of frames
        temp = seed_mask(:,:,end);
        temp_labeled = label(temp,1);
        numrows = max(temp_labeled)*tsize;
        distMat = nan(numrows,2); %this is matrix of all distances matched to frame
        dm_id = 1;
        %          distMap = zeros(size(single(sink_mask))); %this is a movie of all the distance images to check how good it did. only creates if plotFlag = 1
        for tt = 0:tsize-1
            geoframe = bclosing(logical(squeeze(geom_mask(:,:,tt))));
            sinkframe = sink_mask(:,:,tt);
            seedframe = seed_mask(:,:,tt);
            % calc dist map for geom_mask in frame tt
            sinkDist = bwdistgeodesic(logical(geoframe),logical(sinkframe),'quasi-euclidean');
            %             dipshow(sinkDist,'labels')
            % now label seeds from seed_mask===__-----------------
            seedlbl = label(seedframe,1);
            if ~any(seedlbl)
                continue;
            end
            labeledmat = zeros(max(seedlbl),2);
            labeledmat(:,1) = tt;
            if plotflag
                P = false(size(logical(seedlbl)));
                P = imoverlay(P, ~logical(geoframe), [1 1 1]);
                P = imoverlay(P, logical(sinkframe), [0 0 1]);
            end
            for ll = 1:max(seedlbl)
                seedlbl_ll = squeeze((seedlbl == ll));
                seedmsr = measure(seedlbl_ll,0*seedlbl_ll,'Center');
                seeds = floor(seedmsr.Center);
                rows = seeds(1,:)+1;
                cols = seeds(2,:)+1;
                seedDist = bwdistgeodesic(logical(geoframe),rows,cols, 'quasi-euclidean');
%                 seedDist = bwdistgeodesic(logical(geoframe),logical(seedlbl_ll), 'quasi-euclidean');
                D = sinkDist+seedDist;
                % actual distance value should be minimum. save this value
                distval = min(D(:));
                D(isnan(D)) = inf;
                D = round(D * 32) / 32;  % to correct for roundoff errors in bwdistgeodesic, see comment on next line
                %                 mindistmask = D <= (distval + sqrt(eps));   % interestingly this doesn't work, so issue in bwdistgeodesic is really bad (terrible implementation - should just keep track of sqrt(2) indices to fix problem)
                mindistmask = D==min(D(:));
                if distval> prod(size(geoframe)) %sanity check for distance size - prod*size because dipimage doesn't support numel function
                    display(['No connectivity for label # ' num2str(ll) ' in frame ' num2str(tt)]);
                    continue;
                end
                %                 This is code from https://blogs.mathworks.com/steve/2011/12/13/exploring-shortest-paths-part-5/
                %                 However, more straightforward to do as above (don't need
                %                 imregional min because ALL correct paths are minimum.
                %                 paths = imregionalmin(D);
                %                 paths_thinned_many = bwmorph(paths, 'thin', inf); -- thin
                %                 does not do what 'Steve' thinks here. See how to do in
                %                 plotting function
%                 track = typicalShortestPath(sinkDist,[rows,cols],min(D(:)));
                closedmask = bwmorph(mindistmask,'fill',inf);
                paths_thinned_many = bwmorph(closedmask, 'thin', inf);
%                 typicalpaths = dip_image(false(size(sinkDist)));
%                 for ii = 1:size(track,1)
%                    typicalpaths(track(ii,1),track(ii,2)) = true; 
%                 end
                
                
                if plotflag
                    P = imoverlay(P, paths_thinned_many, [.5 .5 .5]);
%                     P = imoverlay(P, logical(typicalpaths), [0 1 0]);
                    P = imoverlay(P, logical(seedlbl_ll), [1 0 0]);
                end
                %                 dist = size(find(paths_thinned_many),1);
                labeledmat(ll,2) = distval;
            end
            if (dm_id + size(labeledmat,1) -1) > size(distMat,1)  % need to make more space
                addons = nan(numrows,2);
                distMat = [distMat;addons];
            end
            distMat(dm_id:(dm_id+size(labeledmat,1)-1),:) = labeledmat;
            dm_id = dm_id+size(labeledmat,1);
            if plotflag
                imshow(P,'InitialMagnification', 'fit','Parent',pathax);
                savefig(pathfig,fullfile(plotsavedir,['MinPath_frame#' num2str(tt) '.fig']));
                saveas(pathfig,fullfile(plotsavedir,['MinPath_frame#' num2str(tt) '.tif']));
            end
            waitbar(tt/(size(seed_mask,3)-1),wb);
        end
        close(wb);
        close(pathfig);
        if saveflag
            distMovie = readtimeseries(fullfile(plotsavedir,['MinPath_frame#'  '*.tif']),'',[],1,0);
            save(fullfile(plotsavedir,'MinPaths'),'distMovie','distMat','sink_mask','seed_mask','geom_mask');
        end
    end
    
     function [img_out,sv_arr] = timedriftCorrect(img_in)
        img_out = 0*img_in;
        imref = squeeze(img_in(:,:,0));
        img_out(:,:,0) = imref;
        sv_arr = nan(2,size(img_in,3)-1);
        for ii = 1:(size(img_in,3)-1)
            imcurr= squeeze(img_in(:,:,ii));
            sv1 = findshift(imref,imcurr,'iter',0);
            shiftim = shift(imcurr,sv1,1);
            img_out(:,:,ii) = shiftim;
            sv_arr(:,ii) = sv1;
        end
     end
     function img_out = applydriftCorrect(img_in,sv_arr)
         img_out = 0*img_in;
         img_out(:,:,0) = img_in(:,:,0);
        for ii = 1:(size(img_in,3)-1)
            currframe = squeeze(img_in(:,:,ii));
            shiftcurrframe = shift(currframe,sv_arr(:,ii));
            img_out(:,:,ii) = shiftcurrframe;            
        end
     end
     function [h,overlayim] = overlay(grey_im,bin_im,cm,mskcol)
         % this function overloads the dipimage overlay method but with a
         % better colormapping
         if nargin<4
             mskcol = [1 0 0]; %make mask perim red
         end
         if nargin<3
             cm = bone(256);
         end
         cm(1,:) = mskcol;
         overlayim = grey_im;
         overlayim(bin_im) = 0;
         h = dipshow(overlayim,cm);
         dipmapping(h,'global');
         dipmapping(h,'lin');
         diptruesize(h,200);
     end
     
     function stitchimage = stitch2images(im1,im2)
         % this functions uses matlab's normxcorr2 function to combine
         % images at the maximum cross correlation position. the larger of
         % the two images serves as the base image and then the smaller
         % image is added on around --> the base image is used in regions of overlap 
         if numel(im1)>numel(im2)
             image = im1;
             template = im2;
         else
             template = im1;
             image = im2;
         end
         gputemplate = gpuArray(template);
         gpuimage = gpuArray(image);
         cc = normxcorr2(template,image);
         [xpeak, ypeak] = find(cc==max(cc(:)));
         xadd = size(template,1) - xpeak;
         yadd = size(template,2) - ypeak;
         stitchimage=zeros(xadd+size(image,1),yadd+size(image,2));
         stitchimage(1:xadd+xpeak,1:yadd+ypeak) = gather(gputemplate);
         stitchimage(xadd+1:end,yadd+1:end) = gather(gpuimage);
     end
     
     function stitchmovie2 = stitch2movies(mov1,mov2)
         % put this on the GPU
         assert(size(mov1,3) == size(mov2,3));
         stitchmovie = cell(1,size(mov1,3));
         lastframe = size(mov1,3);
         parfor ff = 1:lastframe
             im1 = squeeze(mov1(:,:,ff));
             im2 = squeeze(mov2(:,:,ff));
             stitchimage = GeneralAnalysis.stitch2images(im1,im2);
             stitchmovie{ff} = stitchimage;
         end
         xsizes = cellfun(@(x) size(x,1),stitchmovie);
         ysizes = cellfun(@(x) size(x,2),stitchmovie);
         stitchmovie2 = zeros(max(xsizes),max(ysizes),lastframe);
         for ff = 1:lastframe
             currframe = stitchmovie{ff};
             stitchmovie2(1:size(currframe,1),1:size(currframe,2),ff) = currframe;
         end
         % clean up the image
         fulltest = sum(stitchmovie2,3);
         test1 = sum(fulltest,1);
         ylast = find(test1>0,1,'last');
         test2 = sum(fulltest,2);
         xlast = find(test2>0,1,'last');
         stitchmovie2 = stitchmovie2(1:xlast,1:ylast,:);
     end
     
     function [h,overlayim] = viewMaskOverlayPerimStatic(image,mask,cm,mskcol)
         if nargin<4
             mskcol = [1 1 1];
         end
         if nargin<3
             cm = hot(256);
         end
         perim = dt(mask);
         bin_im = (perim==1);
         [h,overlayim] = GeneralAnalysis.overlay(image,bin_im,cm,mskcol);
     end
end
end