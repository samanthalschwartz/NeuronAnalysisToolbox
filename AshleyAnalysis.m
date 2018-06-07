classdef AshleyAnalysis < handle
   properties
    path_channel_cellfill = []; %'F:\Sam\050118 NL1 insertion\cell4-C2.tif'
    path_channel_DHFR = []; %'F:\Sam\050118 NL1 insertion\cell4-C3.tif';
    path_channel_TfR = []; %'F:\Sam\050118 NL1 insertion\cell4-C1.tif';
    path_3ch_datafile = [];
    ROI_trim =[];
    cellFill = [];
    surfaceCargo = [];
    TfR = [];
    inCellSurfaceCargo = [];
    distancematrix = [];
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
               obj.cellFill.loadimage();
               obj.surfaceCargo.loadimage();
               obj.TfR.loadimage();
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
           obj.surfaceCargo.setfilepath(obj.path_channel_DHFR);
           obj.TfR.setfilepath(obj.path_channel_TfR);
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
           newsfmask(:,:,end-3);
         ll =  label(berosion(newsfmask(:,:,end-3))*obj.cellFill.mask(:,:,end-3)) 
       end
       function h = plot_cargo_minFrame(obj)
         [lbl_out] = GeneralAnalysis.labelmask_byframe(obj.surfaceCargo.mask);
         labeledim = lbl_out.*obj.cellFill.mask;
%          lbl_out = GeneralAnalysis.findLabelsInMask(labeledim,obj.cellFill.mask);
         test = min(labeledim,labeledim>0,3);
         test(test>(size(labeledim,3)+1)) = 0;
         dist = dt(obj.cellFill.mask(:,:,0));
         test(dist==1) = max(test)+1;
         %-- make colormap for plotting
         blackjet = flip(jet(255));
         blackjet(1,:) = [0 0 0]; blackjet(end,:) = [1 1 1];
         
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
         c.Label.String = 'First Frame with Cargo Insertion';
         c.Label.FontSize = 16;
         h.OuterPosition = h.OuterPosition + [0 0 400 50];
       end
       function viewCleanedSurfaceCargoMask(obj,cm,mskcol)
          if nargin<3
               mskcol = [1 1 1];
           end
           if nargin<2
               cm = hot(256);
           end
           
           % watershed
           img4wtsd = gaussf(obj.surfaceCargo.image);
           wshed = GeneralAnalysis.watershed_timeseries(-img4wtsd,1);
           newsfmask = obj.surfaceCargo.mask;
           newsfmask(wshed==1) = 0;
%            perim = dt(obj.cellFill.mask);
%            bin_im = (perim==1);
%            bin_im = bin_im*obj.cellFill.mask;
           cargomask = newsfmask*obj.cellFill.mask;
%            cargomask(bin_im) = 1
           [h,overlayim] = channelBase.viewMaskOverlayPerimStatic(cargomask,obj.surfaceCargo.image);
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