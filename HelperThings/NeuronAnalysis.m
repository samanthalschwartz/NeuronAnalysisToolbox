classdef NeuronAnalysis < handle
properties
    path_channel_cellfill = ''; %'F:\Sam\050118 NL1 insertion\cell4-C2.tif'
    path_channel_DHFR = ''; %'F:\Sam\050118 NL1 insertion\cell4-C3.tif';
    path_channel_TfR = ''; %'F:\Sam\050118 NL1 insertion\cell4-C1.tif';
    channel_cellfill = []; % dipimage
    channel_DHFR = []; % dipimage
    channel_TfR = []; % dipimage
    mask_cellfill = [];
    mask_DHFR = [];
    mask_TfR = [];
    cutoff = [];
    soma = [];
    distance_mat = [];
end

methods
    function loadimgfrompaths(obj)
        if ~isempty(obj.path_channel_cellfill) 
            obj.channel_cellfill = obj.loadtiff(obj.path_channel_cellfill);
        end
        if ~isempty(obj.path_channel_DHFR)
            obj.channel_DHFR = obj.loadtiff(obj.path_channel_DHFR);
        end
        if ~isempty(obj.path_channel_TfR)
            obj.channel_TfR  = obj.loadtiff(obj.path_channel_TfR);       
        end
        obj.cropchannels();
    end
    
    function cropchannels(obj)
        if ~isempty(obj.channel_cellfill)
            obj.channel_cellfill = obj.cropimage(obj.channel_cellfill,obj.cutoff);
        end
        if ~isempty(obj.channel_DHFR)
            obj.channel_DHFR = obj.cropimage(obj.channel_DHFR,obj.cutoff);
        end
        if ~isempty(obj.channel_TfR)
            obj.channel_TfR = obj.cropimage(obj.channel_TfR,obj.cutoff);
        end
    end
    
    function make_masks(obj)
        obj.make_mask_cellfill();
        obj.make_mask_DHFR();
        obj.make_mask_TfR();
    end
    function h = maskoverlay_cellfill(obj)
       cm4overlay = hot(256);
       cm4overlay(1,:) = [1 1 1];
       perim = dt(obj.mask_cellfill);
       overlayim = obj.channel_cellfill;
       overlayim(perim==1) = 0;
       h = dipshow(overlayim,cm4overlay);
       dipmapping(h,'lin');
       dipmapping(h,'global');
       diptruesize(h,100);
    end
    
    function h = maskoverlay_DHFR(obj)
       cm4overlay = grey(256);
       cm4overlay(1,:) = [1 0 0];
       overlayim = obj.channel_DHFR;
       overlayim(obj.mask_DHFR) = 0;
       h = dipshow(overlayim,cm4overlay);
       dipmapping(h,'lin');
       diptruesize(h,100);
    end
    
    function make_mask_cellfill(obj)
        if isempty(obj.channel_cellfill)
           if isempty(obj.path_channel_cellfill)
               disp(['No cellfille image to mask'...
               ' Try setting path to .tif file with obj.path_channel_cellfill']);
                return;
           else
               obj.loadimgfrompaths();
           end
        end     
        gsigma = [2 2 1];
        lsigma = [1 1 1];
        LoG = obj.LoGcall(obj.channel_cellfill,gsigma,lsigma);
        [othermasks.mask_raw, othermasks.thresh_raw] = obj.threshold_image(gaussf(obj.channel_cellfill,gsigma));
        [othermasks.mask_LoGcutoff, othermasks.thresh_LoGcutoff] = obj.threshold_image(LoG);
        mask_out = othermasks.mask_raw|othermasks.mask_LoGcutoff;
%         [mask_out,~] = obj.mask_image(obj.channel_cellfill,gsigma);
        obj.mask_cellfill = mask_out;
        disp(['Cellfill image mask set as obj.mask_cellfill']);
    end
    function make_mask_DHFR(obj)
        if isempty(obj.channel_DHFR)
           if isempty(obj.path_channel_DHFR)
               disp(['No DHFR image to mask'...
               ' Try setting path to .tif file with obj.path_channel_cellfill']);
                return;
           else
               obj.loadimgfrompaths();
           end
        end     
        gsigma = 1;
        lsigma = 2;
        LoG = obj.LoGcall(obj.channel_cellfill,gsigma,lsigma);
        [mask, ~] = obj.threshold_image(LoG);
        if ~isempty(obj.mask_cellfill)
            mask = mask*obj.mask_cellfill;
        end
        obj.mask_DHFR = mask;
        disp(['DHFR image mask set as obj.mask_DHFR']);
    end
     function make_mask_TfR(obj)
        if isempty(obj.channel_TfR)
           if isempty(obj.path_channel_TfR)
               disp(['No TfR image to mask'...
               ' Try setting path to .tif file with obj.path_channel_cellfill']);
                return;
           else
               obj.loadimgfrompaths();
           end
        end     
        gsigma = [1 1 0];
        [mask_out,~] = obj.mask_image(obj.channel_TfR,gsigma);
        obj.mask_TfR = mask_out;
        disp(['TfR image mask set as obj.mask_TfR']);
     end
    
     %---
     function h = plot_DHFR_minFrame(obj)
         [labeledim] = obj.labelmask_byframe(obj.mask_DHFR);
         test = min(labeledim,labeledim>0,3);
         test(test>(size(labeledim,3)+1)) = 0;
         %-- make colormap for plotting
         blackjet = flip(jet(255));
         blackjet(1,:) = [0 0 0]; blackjet(end,:) = [0 0 0];
         %-- now plot results
         h = dipshow(test,blackjet);
         dipmapping(h,[0 size(labeledim,3)]);
         diptruesize(h,100);
         % get colorbar tick info
         colorunit = size(labeledim,3)/255;
         numofcolbarval = 4;
         colbarplace =[0:numofcolbarval]*255/4;
         colbarval = floor(colbarplace * colorunit);
         c = colorbar;
         c.Location = 'WestOutside';
         c.Ticks = colbarplace;
         c.TickLabels = colbarval;
         c.FontSize = 16;
         c.Label.String = 'First Frame with DHFR Insertion';
         c.Label.FontSize = 16;
         h.OuterPosition = h.OuterPosition + [0 0 400 50];
     end
     
     function calculate_DHFRdistancefromSoma(obj,plotflag)
         if nargin<2
             plotflag = 0;
         end
         [FILEPATH,NAME,EXT] = fileparts(obj.path_channel_DHFR);
         
         [savefilename, savepathname] = uiputfile('*.mat','Select a directory to save these results',fullfile(FILEPATH, [NAME '_DHFRdistfromSoma_',...
             datestr(datetime('now'),'yy-mm-dd-hh-MM-ss')]));
         % first select what is considered soma is
         if isempty(obj.soma)
             lastframe = size(obj.channel_cellfill,3);
             cellfill_LoG = obj.laplaceofgauss(obj.channel_cellfill,[1 1 0]);
             h = dipshow(cellfill_LoG(:,:,floor(lastframe/2)));
             diptruesize(h,200);
             [soma, vs] = diproi(h);
             obj.soma = soma;
             close(h);
         else
             soma = obj.soma;
         end
         wb = waitbar(0,'Analyzing Distances (this may take a while...)');
         bw2 = soma;
         % make a guess at distance matrix
         temp = obj.mask_DHFR(:,:,floor(size(obj.channel_DHFR,3)/2));
         temp_labeled = label(temp,1);
         numrows = max(temp_labeled)*size(obj.channel_DHFR,3);
         dist_mat = nan(numrows,2);
         dm_id = 1;
         % now loop through each frame, make a mask, label, calculate
         % distance to soma along the cell mask
         for ii = 0:(size(obj.channel_DHFR,3)-1)
             cmframe = closing(closing(obj.mask_cellfill(:,:,ii)));
             maskframe = obj.mask_DHFR(:,:,ii);
             dhfr_labeled = label(maskframe,1);
             D2 = bwdistgeodesic(logical(cmframe), logical(bw2), 'quasi-euclidean');
             if ~any(dhfr_labeled>0)
                 continue;
             else
                 labeledmat = zeros(max(dhfr_labeled),2);
                 labeledmat(:,1) = ii;
                 if plotflag
                 pathfig = figure;
                 P = false(size(logical(dhfr_labeled)));
                 P = imoverlay(P, ~logical(cmframe), [1 1 1]);
                 P = imoverlay(P, logical(bw2), [0 0 1]);
                 end
                 for ll = 1:max(dhfr_labeled)
                     bw1 = squeeze((dhfr_labeled == ll));
                     D1 = bwdistgeodesic(logical(cmframe), logical(bw1), 'quasi-euclidean');
                     D = D1+D2;
                     D = round(D * 8) / 8;
                     D(isnan(D)) = inf;
                     paths = imregionalmin(D);
                     paths_thinned_many = bwmorph(paths, 'thin', inf);
                     if plotflag
%                      P = imoverlay(P, paths, [.5 .5 .5]);
                     P = imoverlay(P, paths_thinned_many, [0 1 0]);
                     P = imoverlay(P, logical(bw1), [1 0 0]);
                     end
                     dist = size(find(paths_thinned_many),1);
                     labeledmat(ll,2) = dist;
                 end
                 if plotflag
                 imshow(P, 'InitialMagnification', 'fit');
                 savefig(pathfig,fullfile(savepathname,[savefilename '_path_frame#' num2str(ii) '.fig']));
                 saveas(pathfig,fullfile(savepathname,[savefilename '_path_frame#' num2str(ii) '.png']));
                 close(pathfig);
                 end
                 if (dm_id + size(size(labeledmat,1))) > size(dist_mat,1)
                     % need to make more space
                     addons = nan(numrows,2);
                     dist_mat = [dist_mat;addons];
                 end
                 dist_mat(dm_id:(dm_id+size(labeledmat,1)-1),:) = labeledmat;
                 dm_id = dm_id+size(labeledmat,1);
             end
             waitbar(ii/(size(obj.channel_DHFR,3)-1),wb);
         end
         close(wb);
         save(fullfile(savepathname,savefilename),'dist_mat');
         obj.dist_mat = dist_mat;
     end
end
methods(Static)
    function img = loadtiff(filepath)
        % this function loads a tiff file into matlab and generates a dipimage
        % must have the bioformats function bfopen: download at https://docs.openmicroscopy.org/bio-formats/5.7.0/developers/matlab-dev.html
        imgbefore = bfopen(filepath);
        img = dip_image(zeros([size(imgbefore{1,1}{1}),size(imgbefore{1,1},1)]));
        for ii = 1:size(img,3)
            img(:,:,ii-1) = imgbefore{1,1}{ii};
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
            gsig = zeros(size(img_in,3),1,3);
            if nargin<2
                lsig = [1 1 0];
            end
        end
        img_g = gaussf(img_in,gsig);
        img_lapl = dxx(img_g,lsig)+dyy(img_g,lsig);
    end
    function img_laplcutoff = imgLaplaceCutoff(img_in,lsig,gsig)
        if nargin<3
            gsig = zeros(size(img_in,3),1,3);
            if nargin<2
                lsig = [1 1 0];
            end
        end
        img_lapl = imgLaplace(img_in,lsig,gsig);
        img_laplcutoff = -img_lapl;
        img_laplcutoff(img_laplcutoff<0) = 0;        
    end
    function mask = imgThreshold(img_in)
        threshval = multithresh(single(img_in),2);
        mask = img_in>=threshval(1);
    end
    function img_out = laplaceofgauss(img_in,gsigma)
        % gaussian filters image and then calculates laplacian
        % inputs: 
            % img_in - dipimage or matrix image (converts to type dipimage)
            % gsigma - optional input to set the gaussian kernal. must be
            %           same dimension as img_in. example: [1 1 0] is gaussian 
            %           smoothing in x and y but not time.
        % outputs:
            % img_out - filtered dipimage. to convert to matlab array use
            %           single(img_out).
            if ~isa(img_in,'dip_image')
                try
                    img_in = dip_image(img_in);
                catch
                    warning('input must be an image matrix');
                    return;
                end
            end
            if nargin<2
                gsigma = [1 1 0];
            end
            %             gim = gaussf(img_in,gsigma);
            gim = img_in;
            img_out = dxx(gim,gsigma) + dyy(gim,gsigma);
    end
    function img_out = laplaceofgauss_mag(img_in,gsigma)
         % gaussian filters image and then calculates laplacian
        % inputs: 
            % img_in - dipimage or matrix image (converts to type dipimage)
            % gsigma - optional input to set the gaussian kernal. must be
            %           same dimension as img_in. example: [1 1 0] is gaussian 
            %           smoothing in x and y but not time.
        % outputs:
            % img_out - filtered dipimage. to convert to matlab array use
            %           single(img_out).
            if ~isa(img_in,'dip_image')
                try
                    img_in = dip_image(img_in);
                catch
                    warning('input must be an image matrix');
                    return;
                end
            end
            if nargin<2
                gsigma = [1 1 0];
            end 
            img_out = NeuronAnalysis.LoGcall(img_in,gsigma,lsigma);
    end
    
    function LoG = LoGcall(im_in,gsigma,lsigma)
        if nargin<3
            lsigma = [1 1 0];
            if nargin<2
                gsigma = [0 0 0]
            end
        end
        gim = gaussf(im_in,gsigma);
        LoG = dxx(gim,lsigma) + dyy(gim,lsigma);
        LoG = -LoG;
        LoG(LoG<0)=0;
    end
    function [mask, threshval] = threshold_image(img_in)
        % threshold image to generate a mask - this can be improved!
        % inputs:
           % img_in - dipimage or matrix image (converts to type dipimage)
           if ~isa(img_in,'dip_image')
                try
                    img_in = dip_image(img_in);
                catch
                    warning('input must be an image matrix');
                    return;
                end
            end
        threshval = multithresh(single(img_in),2);
        mask = img_in>=threshval(1);
    end
    
    function [labeledim] = labelmask_byframe(mask_in)
        % this function labels a binary image/mask and assigns the label
        % value as the frame number
        % inputs:
        %   mask_in - 3D dipimage
        % outputs:
        %   labeledim - label dipimage object
        labeledim = dip_image(zeros(size(mask_in)));
        if numel(size(mask_in))<3
            disp('This method is not useful for 2D images, because it will just return the mask. Try with a 3D image series.')
            return;
        end
        for ii = 1:size(mask_in,3)
            temp = label(mask_in(:,:,ii-1),1);
            temp(temp>0) = ii;
            labeledim(:,:,ii-1) = temp;
        end
    end  
    
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
    
end
end