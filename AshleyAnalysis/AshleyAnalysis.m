classdef AshleyAnalysis < handle
   properties
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
    cleanedcargomask = [];
    distmask = [];
    distmaskFull = [];
    distmaskAIS = [];
    cargo_heatmap = [];
    % first frame of baseline (ie: 1);
    % last frame of baseline (ie: 6, frame 6 counted as baseline - occurs at t=0)
     % frame rate in min/frame (ie: 1, 1 minute interval between frames)
    imagingparams = struct(...
        'baseline', struct('frame_start',[],'frame_end',[], 'framerate',[]),... 
        'releaseframe',[],... % duration of 405 light in min (ie: 1, means 1 minute between end of baseline and first frame of post release-- this dictates the time at postreleaseframe_start 
        'postrelease', struct('frame_start',[],'frame_end',[], 'framerate',[])); % first frame of baseline (ie: 7, frame 7 is first frame of post release, occurs at time = releasetime
   end
   
   methods
       function obj = AshleyAnalysis()
           obj.cellFill = channelCellFill();
           obj.surfaceCargo = channelSurfaceCargo();
           obj.TfR = channelTfR();
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
       function recrop(obj)
           
           
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
       function h = calc_cargo_minFrame(obj,maxtime,perimbool)
           % cellperim is boolean for including cell perimeter in image
           if ~isempty(obj.cleanedcargomask)
               maskimg = dip_image(logical(obj.cleanedcargomask));
           else
               maskimg =  dip_image(logical(obj.surfaceCargo.mask));
           end
           
           if isfield(obj.imagingparams,'baselineframe_start') && ~isempty(obj.imagingparams.baselineframe_start)
               startfrm = obj.imagingparams.postrelease(1).frame_start;
               postrelease = obj.imagingparams.postrelease;
               endfrm = postrelease(end).frame_end;
               if nargin>1
                   endfrm = maxtime;
               elseif ischar(endfrm)
                   if strcmp(endfrm,'end')
                       endfrm = size(obj.cellFill.image,3);
                       postrelease(end).frame_end = endfrm;
                   end
               end
               if isa(maskimg,'dip_image')
                   [lbl_out] = GeneralAnalysis.labelmask_byframe(maskimg(:,:,startfrm-1:endfrm-1));
               else
                   [lbl_out] = GeneralAnalysis.labelmask_byframe(maskimg(:,:,startfrm:endfrm));
               end
               
               if isa(obj.cellFill.mask,'dip_image')
                   cfmask = obj.cellFill.mask(:,:,startfrm-1:endfrm-1);
               else
                   cfmask = obj.cellFill.mask(:,:,startfrm:endfrm);
               end
               firstfrmtime = obj.imagingparams.releaseframe;
               
               currfirst = firstfrmtime;
               timerange = [];
               for ff = 1:numel(postrelease)
               duration = (postrelease(ff).frame_end - postrelease(ff).frame_start)*postrelease(ff).framerate;
               newend = duration+currfirst;
               timerange = [timerange, currfirst:postrelease(ff).framerate:newend];
               currfirst = newend;
               end
               
           else %no imaging params set
               [lbl_out] = GeneralAnalysis.labelmask_byframe(maskimg);
               if nargin>1
                   endfrm = maxtime;
               else
                   endfrm = size(obj.cellFill.image,3);
               end
               cfmask = obj.cellFill.mask;
               timerange  = 1:endfrm;
           end
           labeledim = lbl_out.*cfmask;
           %          lbl_out = GeneralAnalysis.findLabelsInMask(labeledim,obj.cellFill.mask);
           test = min(labeledim,labeledim>0,3);
           test(test>(size(labeledim,3)+1)) = 0;
           if nargin==3
               dist = dt(cfmask(:,:,0));
               test(dist==1) = max(test)+1;
           end
           obj.cargo_heatmap.image = test;
           obj.cargo_heatmap.timerange = timerange;
           
       end
       function h = plotCargoHeatMap(obj,reset)
           if isempty(obj.cargo_heatmap) || nargin>1
               obj.calc_cargo_minFrame();
           end
           fixedtimerange = 160;
           %-- make colormap for plotting
           blackjet = flip(jet(255));
           blackjet(1,:) = [0 0 0];% blackjet(end,:) = [1 1 1];
           
           
           % get colorbar tick info
           %            colorunit = size(labeledim,3)/255;
%            colorunit = obj.cargo_heatmap.timerange(end)/255;
           colorunit = fixedtimerange/255;
           numofcolbarval = 4;
           colbarplace =[0:numofcolbarval]*255/4;
           colbarval = [floor(colbarplace * colorunit)];
           
           %-- now plot results
           h = dipshow(obj.cargo_heatmap.image,blackjet);
           dipmapping(h,[obj.imagingparams.releaseframe size(obj.cargo_heatmap.timerange,2)]);
           diptruesize(h,120);
           
           
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
           h.OuterPosition = h.OuterPosition + [0 0 400 50]; 
       end
       
       function M = plotDensityperTime(obj,distances)
           %input:  distances in microns as an 1 x n vector of max values, values between are used
           %    example: distances = [100,200,300, inf]; interval is
           %    0<val<=100, 100<val<=200, 200<val<=300, val>300;
           % Average intensity within the distance is
           % make the distance mask
           
           if nargin<2
               distances = [18/obj.pxsize 100/obj.pxsize 200/obj.pxsize];
           end
           clear M;
           if isempty(obj.distmask)
               obj.makeDistanceMask();
           end
           temp = obj.distmask; temp(temp==Inf)=0;
           maxdist = max(temp);
           interestmsk = obj.cleanedcargomask.*obj.cellFill.mask;
           sfim = obj.surfaceCargo.image*interestmsk;
           scsums = sort(single(squeeze(sum(sfim,[],[1 2]))));
           M.maxintensity = mean(scsums(end-2:end));
           scmask = sum(dip_image(obj.cleanedcargomask),[],3);
           % now go through for each time point and calculate densities.
           bgmaskthin = isnan(obj.distmask) & ~scmask;
           bgmask = berosion(bgmaskthin,5);
           bgmask = repmat(bgmask,1,1,size(sfim,3));
           fullbackgroundimage = GeneralAnalysis.regionfill_timeseries(obj.surfaceCargo.image*bgmask,~bgmask);
           backgroundimage = fullbackgroundimage*interestmsk;
           M.distance = distances;
           M.rawintensity = zeros(numel(distances),size(sfim,3));
           M.areanormintensity = zeros(numel(distances),size(sfim,3));
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
               M.rawintensity(ii,:) = movmean(origvals,windowsize);
               M.areanormintensity(ii,:) = M.rawintensity(ii,:)/mskarea;
           end
           obj.M = M;
       end
       
       function M = plotAreaperTime(obj,distances)
           %input:  distances in microns as an 1 x n vector of max values, values between are used 
           %    example: distances = [100,200,300, inf]; interval is
           %    0<val<=100, 100<val<=200, 200<val<=300, val>300;
           % Average intensity within the distance is 
           % make the distance mask
           
           if nargin<2
           distances = [50/obj.pxsize 100/obj.pxsize 200/obj.pxsize];
           end
           clear M;
          
           distmask = obj.makeDistanceMask();
           interestmsk = obj.surfaceCargo.mask*obj.cellFill.mask;
           sfim = obj.surfaceCargo.image*interestmsk;
           scsums = sort(single(squeeze(sum(sfim,[],[1 2]))));
           M.maxintensity = mean(scsums(end-2:end));
           scmask = sum(dip_image(obj.surfaceCargo.mask),[],3);
           % now go through for each time point and calculate densities.
           bgmaskthin = isnan(distmask) & ~scmask;
           bgmask = berosion(bgmaskthin,5);
           bgmask = repmat(bgmask,1,1,size(sfim,3));
           fullbackgroundimage = GeneralAnalysis.regionfill_timeseries(obj.surfaceCargo.image*bgmask,~bgmask);
           backgroundimage = fullbackgroundimage*interestmsk;
           M.distance = distances;
           M.rawintensity = zeros(numel(distances),size(sfim,3));
           M.areanormintensity = zeros(numel(distances),size(sfim,3));
           for ii = 1:numel(distances)
               
               if ii == 1
                   currmask = distmask<=distances(ii);
               else
                   currmask = distmask>distances(ii-1) & distmask<=distances(ii);
               end
               mskarea = sum(currmask);
               newsfmask = repmat(currmask,1, 1, size(sfim,3));
               sumcargoinmask = sum(sfim,newsfmask,[1 2]);
               sumbginmask = sum(backgroundimage,newsfmask,[1 2]);
               thisplot = squeeze(sumcargoinmask) - squeeze(sumbginmask);
               % sliding window
               origvals = single(squeeze(thisplot));
               windowsize = 4;
               M.rawintensity(ii,:) = movmean(origvals,windowsize);
               M.areanormintensity(ii,:) = M.rawintensity(ii,:)/mskarea;
           end
           
       end
       
       function makeDistanceMasks(obj)
           % make the distance mask for the full cellfill mask area (AIS included)
           sumsFull = bdilation(obj.cellFill.mask,1);
           geoframeFull = sum(sumsFull,[],3);
           obj.distmaskFull = dip_image(bwdistgeodesic(logical(geoframeFull),logical(sinkframe),'quasi-euclidean'));
           clear sumsFull geoframeFull;
           if ~isempty(obj.cellFill.AIS_mask)
               % make the distance mask (no AIS included if it is there)
               sums = bdilation(obj.cellFill.mask-obj.cellFill.AIS_mask,1);
               geoframe = sum(sums,[],3);
               sinkframe = squeeze(obj.cellFill.soma_mask(:,:,1));
               obj.distmask = dip_image(bwdistgeodesic(logical(geoframe),logical(sinkframe),'quasi-euclidean'));
               clear sums geoframe;
               % make the distance mask for the AIS only
               sumsAIS = bdilation(obj.cellFill.AIS_mask,1);
               geoframeAIS = sum(sumsAIS,[],3);
               obj.distmaskAIS = dip_image(bwdistgeodesic(logical(geoframeAIS),logical(sinkframe),'quasi-euclidean'));
               clear sumsAIS geoframeAIS;
           end
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
           if isprop(obj,'cleanedcargomask') && ~isempty(obj.cleanedcargomask) && nargin<2
               cargomask = GeneralAnalysis.cleanUpMask_manual_square(overlayim,obj.cleanedcargomask.*obj.cellFill.mask,100);
           else
               cargomask = GeneralAnalysis.cleanUpMask_manual_square(overlayim,single(obj.surfaceCargo.mask).*obj.cellFill.mask,100);
           end
           obj.cleanedcargomask = cargomask;
       end
       function cleanSurfaceCargoMaskbyFrame_Manual(obj,reset)
           overlayim = obj.surfaceCargo.image;
           if isprop(obj,'cleanedcargomask') && ~isempty(obj.cleanedcargomask) && nargin<2
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