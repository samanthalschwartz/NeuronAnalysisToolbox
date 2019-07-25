classdef AshleyAnalysis < handle
    properties
        savename = [];
        path_channel_cellfill = []; %'F:\Sam\050118 NL1 insertion\cell4-C2.tif'
        path_channel_surfaceCargo = []; %'F:\Sam\050118 NL1 insertion\cell4-C3.tif';
        path_channel_TfR = []; %'F:\Sam\050118 NL1 insertion\cell4-C1.tif';
        path_3ch_datafile = [];
        ROI_trim =[];
        cellFill = [];
        surfaceCargo = [];
        TfR = [];
        inCellSurfaceCargo = [];
        distancematrix = [];
        pxsize = 1/3.5;
        M = [];
        M_AIS = [];
        M_noAIS = [];
        cleanedcargomask = [];
        distmask = [];
        distmaskPart = [];
        distmaskAIS = [];
        cargo_heatmap = [];
        % release frame is first frame post release (so surface cargo in
        % this frame counts)
        % first frame of baseline (ie: 1);
        % last frame of baseline (ie: 6, frame 6 counted as baseline - occurs at t=0)
        % frame rate in min/frame (ie: 1, 1 minute interval between frames)
        imagingparams = struct(...
            'baseline', struct('frame_start',[],'frame_end',[], 'framerate',[]),...
            'releaseframe',[],... % duration of 405 light in min (ie: 1, means 1 minute between end of baseline and first frame of post release-- this dictates the time at postreleaseframe_start
            'postrelease', struct('frame_start',[],'frame_end',[], 'framerate',[])); % first frame of baseline (ie: 7, frame 7 is first frame of post release, occurs at time = releasetime
        
        guiFig;
    end
   
   methods
       function obj = AshleyAnalysis()
           obj.cellFill = channelCellFill();
           obj.surfaceCargo = channelSurfaceCargo();
           obj.TfR = channelTfR();
       end
       function load(obj)
          % cast all dipimage objects back to dip_images or logicals
          try
          obj.cellFill.soma_mask = logical(obj.cellFill.soma_mask);
          catch
          end
          try
          obj.cellFill.image = dip_image(obj.cellFill.image);
          catch
          end
          try
          obj.cellFill.mask = logical(obj.cellFill.mask);
          catch
          end
          try
          obj.surfaceCargo.image = dip_image(obj.surfaceCargo.image);
          catch
          end
          try
          obj.surfaceCargo.mask = logical(obj.surfaceCargo.mask);
          catch
          end
          try
          obj.cargo_heatmap = dip_image(obj.cargo_heatmap); 
          catch
          end
       end
       function loadImages(obj)
           if ~isempty(obj.path_3ch_datafile)
               im_array = loadtiff_3ch(obj.path_3ch_datafile);
               obj.cellFill.setimage(im_array(:,:,:,1));
               obj.surfaceCargo.setimage(im_array(:,:,:,2));
               obj.TfR.setimage(im_array(:,:,:,3));
           else
               obj.cellFill.setfilepath(obj.path_channel_cellfill);
               obj.surfaceCargo.setfilepath(obj.path_channel_surfaceCargo);
               obj.TfR.setfilepath(obj.path_channel_TfR);
               if ~isempty(obj.path_channel_cellfill)
                   obj.cellFill.loadimage();
               end
               if ~isempty(obj.path_channel_surfaceCargo)
                   obj.surfaceCargo.loadimage();
               end
               if ~isempty(obj.path_channel_TfR)
                   obj.TfR.loadimage();
               end
           end
       end
       function setCellFill(obj,image_in)
           obj.cellFill.setimage(obj,image_in)
       end
       
       function setSurfaceCargo(obj,image_in)
           obj.surfaceCargo.setimage(obj,image_in)
       end
       function setTfR(obj,image_in)
           obj.TfR.setimage(obj,image_in)
       end
       function setFilePaths(obj)
           obj.cellFill.setfilepath(obj.path_channel_cellfill);
           obj.surfaceCargo.setfilepath(obj.path_channel_surfaceCargo);
           obj.TfR.setfilepath(obj.path_channel_TfR);
       end
       function cropimage(obj)
           try
               obj.cellFill.trim_rawimage(1);
               obj.surfaceCargo.ROI_trim = obj.cellFill.ROI_trim;
               obj.surfaceCargo.cropimagefromtrimROI();
           catch
               disp('Cropping aborted, leaving image as it was');
           end
           resetimages(obj);
       end
       function resetimages(obj)
            obj.cleanedcargomask = [];
            obj.cargo_heatmap = [];
            obj.cellFill.resetImage();
            obj.cellFill.resetCellFill();
            obj.surfaceCargo.resetImage();
       end
       function clearimages(obj)
          obj.resetimages;
          obj.cellFill.clearImage;
          obj.surfaceCargo.clearImage;
       end
       function maskCargoInsideCell(obj)
           obj.inCellSurfaceCargo = obj.surfaceCargo.mask.*obj.cellFill.mask;
       end
       function maskImages(obj)
           obj.cellFill.mask_img();
           obj.surfaceCargo.mask_img();
       end
       function showMasks(obj)
          obj.cellFill.viewMaskOverlayPerim();
          obj.surfaceCargo.viewMaskOverlayFill();
       end
       function ll = labelExample(obj)
           img4wtsd = gaussf(obj.surfaceCargo.image);
           wshed = GeneralAnalysis.watershed_timeseries(-img4wtsd,1);
           newsfmask = obj.surfaceCargo.mask;
           newsfmask(wshed==1) = 0;
         ll =  label(berosion(newsfmask(:,:,end-3))*obj.cellFill.mask(:,:,end-3)) ;
       end
       function h = calc_cargo_minFrame(obj)
           % cellperim is boolean for including cell perimeter in image
           if ~isempty(obj.cleanedcargomask)
               maskimg = dip_image(logical(obj.cleanedcargomask));
           else
               maskimg =  dip_image(logical(obj.surfaceCargo.mask));
           end
           
           [lbl_out] = GeneralAnalysis.labelmask_byframe(maskimg);
           
           cfmask = obj.cellFill.mask;
           
           labeledim = lbl_out.*cfmask;
           %          lbl_out = GeneralAnalysis.findLabelsInMask(labeledim,obj.cellFill.mask);
           test = min(labeledim,labeledim>0,3);
           test(test>(size(labeledim,3)+1)) = 0;
           if nargin==3
               dist = dt(cfmask(:,:,0));
               test(dist==1) = max(test)+1;
           end
           obj.cargo_heatmap = test;
           
       end
       function [h,im] = plotCargoHeatMap(obj,reset,imgparams)
% input params:
% -- reset: boolean 0 or 1 to recalculate with cal_cargo_minFrame
% -- imgparams: structure with
%           imgparams.maxtime - maxtime range for the colomap
%           imgparams.colmap - input colormap

           if isempty(obj.cargo_heatmap)
               obj.calc_cargo_minFrame();
           end
           if  nargin>1 && reset
               obj.calc_cargo_minFrame();
           end
                   
           % release frame is first frame post release (so surface cargo in
           % this frame counts)
           working_image = obj.cargo_heatmap - (obj.imagingparams.releaseframe - 1);
           working_image(working_image<0) = 0;
               blackjet = flip(jet(255));
               blackjet(1,:) = [0 0 0];% blackjet(end,:) = [1 1 1];
               colmap = blackjet;
               if nargin<3
                   maxtime = 160;
               else
                   if isfield(imgparams,'maxtime')
                       maxtime  = imgparams.maxtime;
                   end
                   if isfield(imgparams,'colmap')
                       colmap = imgparams.colmap;
                   end
               end
               
         
           %-- make colormap for plotting
         
           
           % get colorbar tick info
           %            colorunit = size(labeledim,3)/255;
%            colorunit = obj.cargo_heatmap.timerange(end)/255;
           colorunit = maxtime/255;
           numofcolbarval = 4;
           colbarplace =[0:numofcolbarval]*255/4;
           colbarval = [floor(colbarplace * colorunit)];
           %-- now plot results
           h = dipshow(working_image,blackjet);
           top_diprange = maxtime ./ obj.imagingparams.postrelease.framerate;
           im = working_image;
           im(im>top_diprange) = 0;
           dipmapping(h,[0 top_diprange]);
%            dipmapping(h,[0 max(working_image(:))]);
           diptruesize(h,80);
           
           
           c = colorbar;
           c.Location = 'WestOutside';
           c.Ticks = colbarplace;
           c.TickLabels = colbarval;
           c.TickLength = 0.02
           c.FontSize = 18;
           c.LineWidth = 1.5;
           c.FontName = 'Arial';
           %c.FontWeight = 'bold';
           c.Label.String = 'time of first appeareance after release (min)';
           c.Label.FontSize = 18;
           c.Label.FontName = 'Arial';
           h.OuterPosition = h.OuterPosition + [-300 -300   300   300]; 
%            set(findall(h,'-property','FontName'),'FontName','Arial')
       end
       
       function calcDensityperTime(obj,distances)
           %input:  distances in microns as an 1 x n vector of max values, values between are used
           %    example: distances = [100,200,300, inf]; interval is
           %    0<val<=100, 100<val<=200, 200<val<=300, val>300;
           % Average intensity within the distance is
           % make the distance mask
           
           if nargin<2
               distances = [5 40 200];
           end
           if isempty(obj.distmask)
               obj.makeDistanceMask();
           end
           if isempty(obj.cellFill.AIS_mask)
               regions = {'Total'};
           else
               regions = {'Total','No AIS','AIS only'};
           end
           if isempty(obj.surfaceCargo.bgsub_image)
               obj.makeBackgroundSubSurfaceCargoImage();
           end
           currM.distance = distances;
           currM.rawintensity = zeros(numel(distances),size(obj.surfaceCargo.bgsub_image,3));
           currM.areanormintensity = zeros(numel(distances),size(obj.surfaceCargo.bgsub_image,3));
           for rr = 1:numel(regions)
               switch regions{rr}
                   case 'Total'
                       surfaceCargoMask = obj.cleanedcargomask;
                   case 'No AIS'
                       surfaceCargoMask = obj.cleanedcargomask.*~obj.cellFill.AIS_mask;
                   case 'AIS only'
                       surfaceCargoMask = obj.cleanedcargomask.*obj.cellFill.AIS_mask;
               end
               for ii = 1:numel(distances)
                   if ii == 1
                       currmask = obj.distmask<=distances(ii);
                   else
                       currmask = obj.distmask>distances(ii-1) & obj.distmask<=distances(ii);
                   end
                   themask = currmask.*surfaceCargoMask;
                   signal = sum(obj.surfaceCargo.bgsub_image,themask,[1 2]);
                   signal = single(squeeze(signal));
                   mask_area = sum(themask(:)); % area in number of pixels
                   % sliding window
                   windowsize = 4;
                   currM.rawintensity(ii,:) = movmean(signal,windowsize);
                   currM.areanormintensity(ii,:) = currM.rawintensity(ii,:)/mask_area;
                   currM.mask_area(ii) = mask_area;
               end
               switch regions{rr}
                   case 'Total'
                       obj.M = currM;
                   case 'No AIS'
                       obj.M_noAIS = currM;
                   case 'AIS only'
                       obj.M_AIS = currM;
               end
           end
       end
       function plotDensityperTime(obj,distances)
           %input:  distances in microns as an 1 x n vector of max values, values between are used
           %    example: distances = [100,200,300, inf]; interval is
           %    0<val<=100, 100<val<=200, 200<val<=300, val>300;
           % Average intensity within the distance is
           % make the distance mask
           
           if nargin<2
               distances = [5 40 200];
           end
           if isempty(obj.distmask)
               obj.makeDistanceMask();
           end
           if isempty(obj.cellFill.AIS_mask)
               regions = {'Total'};
           else
               regions = {'Total','No AIS','AIS only'};
           end
           
           for rr = 1:numel(regions)
               switch regions{rr}
                   case 'Total'
                       surfaceCargoMask = obj.cleanedcargomask;
                   case 'No AIS'
                       surfaceCargoMask = obj.cleanedcargomask-obj.cellFill.AIS_mask;
                   case 'AIS only'
                       surfaceCargoMask = obj.cleanedcargomask.*obj.cellFill.AIS_mask;
               end
               interestmsk = dip_image(surfaceCargoMask.*obj.cellFill.mask);
               sfim = obj.surfaceCargo.image*interestmsk;
               scsums = sort(single(squeeze(sum(sfim,[],[1 2]))));
               currM.maxintensity = mean(scsums(end-2:end));
               scmask = sum(dip_image(surfaceCargoMask),[],3);
               % now go through for each time point and calculate densities.
               bgmaskthin = isnan(obj.distmask) & ~scmask;
               bgmask = berosion(bgmaskthin,5);
               bgmask = repmat(bgmask,1,1,size(sfim,3));
               fullbackgroundimage = GeneralAnalysis.regionfill_timeseries(obj.surfaceCargo.image*bgmask,~bgmask);
               backgroundimage = fullbackgroundimage*interestmsk;
               currM.distance = distances;
               currM.rawintensity = zeros(numel(distances),size(sfim,3));
               currM.areanormintensity = zeros(numel(distances),size(sfim,3));
               for ii = 1:numel(distances)
                   if ii == 1
                       currmask = obj.distmask<=distances(ii);
                   else
                       currmask = obj.distmask>distances(ii-1) & obj.distmask<=distances(ii);
                   end
                   mskarea = sum(currmask(:));
                   newsfmask = repmat(currmask,1, 1, size(sfim,3));
                   sumcargoinmask = sum(sfim,newsfmask,[1 2]);
                   sumbginmask = sum(backgroundimage,newsfmask,[1 2]);
                   thisplot = squeeze(sumcargoinmask) - squeeze(sumbginmask);
                   % sliding window
                   origvals = single(squeeze(thisplot));
                   windowsize = 4;
                   currM.rawintensity(ii,:) = movmean(origvals,windowsize);
                   currM.areanormintensity(ii,:) = currM.rawintensity(ii,:)/mskarea;
               end
               
               switch regions{rr}
                   case regions{1}
                       obj.M = currM;
                   case regions{2}
                       obj.M_noAIS = currM;
                   case regions{3}
                       obj.M_AIS = currM;
               end
           end
       end
       function makeBackgroundSubSurfaceCargoImage(obj)
           %--- calculate background of SurfaceCargo Image
           max_cleanedcargomask = max(dip_image(obj.cleanedcargomask),[],3); %max project down to 2d imgae
           bgmaskthin = isnan(obj.distmask) & ~max_cleanedcargomask; %get region without cleaned cargo mask and cellfill
           bgmask = berosion(bgmaskthin,5); %errode mask
           bgmask = repmat(bgmask,1,1,size(obj.cleanedcargomask,3)); %remake 3d
           
           fullbackgroundimage = GeneralAnalysis.regionfill_timeseries(obj.surfaceCargo.image*bgmask,~bgmask);% use regionfill to generate bg image
           obj.surfaceCargo.bgsub_image = obj.surfaceCargo.image-fullbackgroundimage;
           obj.surfaceCargo.bgsub_image(obj.surfaceCargo.bgsub_image<0) = 0;
       end
       function save(obj,savename)
           if nargin < 2
               disp('Need to add a savepath');
               return;
           end
          % cast all dipimage objects back to singles
          obj.cellFill.soma_mask = single(obj.cellFill.soma_mask);
          obj.cellFill.image = single(obj.cellFill.image);
          obj.cellFill.mask = single(obj.cellFill.mask);
          obj.surfaceCargo.image = single(obj.surfaceCargo.image);
          obj.surfaceCargo.mask = single(obj.surfaceCargo.mask);
%           obj.cargo_heatmap = single(obj.cargo_heatmap);
           obj.guiFig =[];
          aa = obj;
          
          save(savename,'aa');
       end
       function makeDistanceMask(obj)
           if isempty(obj.cellFill.mask)
               errordlg('Need to make a Cell Region Mask');
               return
           end
           if isempty(obj.cellFill.soma_mask)
               obj.cellFill.selectSoma;
           end
           sinkframe = squeeze(obj.cellFill.soma_mask(:,:,1));
           % make the distance mask for the full cellfill mask area (AIS included)
           sums = bdilation(logical(obj.cellFill.mask),1);
           geoframe = sum(sums,[],3);
           obj.distmask = dip_image(bwdistgeodesic(logical(geoframe),logical(sinkframe),'quasi-euclidean'));
           clear sums geoframe;
       end

       function [h,lagim] = plot_cargo_minFrameMovie(obj,savename, framelag)
           %            cellperim is boolean for including cell perimeter in image
           if nargin<3
               framelag = 4;
           end
           if nargin<2
               savename = fullfile(pwd,'movie');
           end
           [lbl_out] = GeneralAnalysis.labelmask_byframe(obj.cleanedcargomask);
           labeledim = lbl_out.*obj.cellFill.mask;
           %          lbl_out = GeneralAnalysis.findLabelsInMask(labeledim,obj.cellFill.mask);
           test = min(labeledim,labeledim>0,3);
           test(test>(size(labeledim,3)+1)) = 0;
           
           blackjet = flip(jet(max(test)));
           blackjet(1,:) = [0 0 0]; blackjet(end,:) = [1 1 1];
           videoout = labeledim*0;
           % make video per frame
           for ff = 0:(max(test)-1)
               videoout(:,:,ff) = (test == (ff+1));
           end
           % now make arbitrary sliding windw
           lagim = 0*videoout;
           for ff = 1:size(videoout,3)
               if ff<=framelag
                   lagim(:,:,ff-1) = sum(videoout(:,:,0:ff-1),[],3);
               else
                   lagim(:,:,ff-1) = sum(videoout(:,:,ff-framelag-1:ff-1),[],3);
               end
           end
           
           lagim(lagim~=0) = 30;
           tcm = jet(30);
           tcm(30,:) = [0 1 1];
           tcm(1,:) = [0 0 0];
           h = dipshow(lagim,tcm);
           diptruesize(h,100);
           
           %save the movie
           options.framerate = 18; %default
           writeDipImageMovie(h,savename,options)
       end
       
       function calculateSurfaceCargoDistances(obj,plotflag,savedir)
           if isprop(obj,'cleanedcargomask') && ~isempty(obj.cleanedcargomask) 
               seed_mask = obj.cleanedcargomask;
           else
               if isempty(obj.inCellSurfaceCargo)
                   obj.maskCargoInsideCell;
               end
               seed_mask = logical(obj.inCellSurfaceCargo);
           end
           if isempty(obj.cellFill.mask_thick)
               obj.cellFill.make_thickMask();
           end
           sink_mask = logical(obj.cellFill.fullsoma_mask);
           geom_mask = logical(obj.cellFill.mask_thick>0);
           if nargin==3
               distMat = GeneralAnalysis.geodesic_seedDistfromMask(sink_mask,seed_mask,geom_mask,plotflag,savedir);
           elseif nargin==2
               distMat = GeneralAnalysis.geodesic_seedDistfromMask(sink_mask,seed_mask,geom_mask,plotflag);
           else
               distMat = GeneralAnalysis.geodesic_seedDistfromMask(sink_mask,seed_mask,geom_mask);
           end
           obj.distancematrix = distMat;
       end
       
       function cleanSurfaceCargoMask(obj)
           img4wtsd = gaussf(obj.surfaceCargo.image);
           wshed = GeneralAnalysis.watershed_timeseries(-img4wtsd,1);
           newsfmask = obj.surfaceCargo.mask;
           newsfmask(wshed==1) = 0;
           cargomask = newsfmask*obj.cellFill.mask;
           obj.cleanedcargomask = cargomask;
       end
       function cleanSurfaceCargoMask_Manual(obj,reset)
           % uses current cleaned surface mask to start with if it exists
           % any input in the function call acts as a reset
%            [h,overlayim] = GeneralAnalysis.viewMaskOverlayPerimStatic(obj.surfaceCargo.image,obj.cellFill.mask);
%            close(h);
%            overlayim(overlayim==0) = max(overlayim)*2;
            % - for resetting image use 1
           overlayim = obj.surfaceCargo.image;
           if isprop(obj,'cleanedcargomask') && ~isempty(obj.cleanedcargomask) && nargin<2 || ~reset
               cargomask = GeneralAnalysis.cleanUpMask_manual_square(overlayim,obj.cleanedcargomask.*obj.cellFill.mask,100);
           else
               cargomask = GeneralAnalysis.cleanUpMask_manual_square(overlayim,single(obj.surfaceCargo.mask).*obj.cellFill.mask,100);
           end
           obj.cleanedcargomask = cargomask;
       end
       function cleanSurfaceCargoMaskbyFrame_Manual(obj,reset)
           overlayim = obj.surfaceCargo.image;
           if isprop(obj,'cleanedcargomask') && ~isempty(obj.cleanedcargomask) && nargin<2 || ~reset
               cargomask = GeneralAnalysis.cleanUpMask_byframe_square(overlayim,obj.cleanedcargomask.*obj.cellFill.mask,100);
           else
               cargomask = GeneralAnalysis.cleanUpMask_byframe_square(overlayim,single(obj.surfaceCargo.mask).*obj.cellFill.mask,100);
           end
           obj.cleanedcargomask = dip_image(logical(cargomask));
       end
       function cleanCellFillMask_Manual(obj)
           overlayim = obj.cellFill.image;
           mask = GeneralAnalysis.cleanUpMask_manual_square(overlayim,obj.cellFill.mask);
           obj.cellFill.mask = mask;
       end
       function [h,overlayim] = viewCleanedSurfaceCargoMask(obj,cm,mskcol)
           if nargin<3
               mskcol = [1 1 1];
           end
           if nargin<2
               cm = hot(256);
           end
           
           if isempty(obj.cleanedcargomask)
               obj.cleanSurfaceCargoMask;
           end
           [h,overlayim] = channelBase.viewMaskOverlayPerimStatic(obj.cleanedcargomask,obj.surfaceCargo.image);
       end
       
       
       function viewsurfaceCargowithCellFillOutline(obj,cm,mskcol)
           if nargin<3
               mskcol = [1 1 1];
           end
           if nargin<2
               cm = hot(256);
           end
           
          
           perim = dt(obj.cellFill.mask);
           bin_im = (perim==1);
           bin_im;
           [h,overlayim] = GeneralAnalysis.overlay(obj.surfaceCargo.image,bin_im,cm,mskcol)
%            [h,overlayim] = channelBase.viewMaskOverlayPerimStatic(cargomask,obj.surfaceCargo.image);
           
       end
   end
   
end