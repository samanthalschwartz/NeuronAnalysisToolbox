classdef channelCellFill < channelBase
    properties
        mask_thick
        lsig = [1 1 0];
        gsig = [1 1 0];
        soma_mask = [];
        soma_vertices = [];
        AIS_mask = [];
        AIS_vertices = [];
        fullsoma_mask = [];
        fullsoma_vertices = [];
    end
    
    methods
       function obj = channelCellFill()
           obj = obj@channelBase();
       end
       function mask_img(obj)
           %         img_m = medif(obj.image,3);
           %         img_g = GeneralAnalysis.imgGauss(img_m,obj.gsig);
           %         img_laplcutoff = GeneralAnalysis.imgLaplaceCutoff(img_m,obj.lsig,obj.gsig);
           img_g = GeneralAnalysis.imgGauss(obj.image,obj.gsig);
           gmask = GeneralAnalysis.imgThreshold(img_g);
           img_laplcutoff = GeneralAnalysis.imgLaplaceCutoff(obj.image,obj.lsig,obj.gsig);
           [lmask,threshval] = GeneralAnalysis.imgThreshold_fixedUserInput(img_laplcutoff);
           mask_out = gmask|lmask;
           obj.mask = mask_out;
           obj.backgroundvalue = threshval;
       end       function mask_img_other(obj)
           img_1 = medif(obj.image,3);
           img_2 = gaussf(img_1,obj.gsig);
           imglcutoff = img_2;
           img_3 = GeneralAnalysis.imgLaplaceCutoff(imglcutoff,obj.gsig,obj.lsig);
           img4mask = img_3.^1.2; 
           maskpre1 = GeneralAnalysis.imgThreshold(img4mask);
           mask = GeneralAnalysis.bwmorph_timeseries(maskpre1,'thicken',1);
           obj.mask = logical(mask);
       end
             
       function mask_img_lsnr(obj)
           img_m = medif(obj.image,4);
           [out,~,~] = backgroundoffset(img_m,0,15,20,15);
           img_g = GeneralAnalysis.imgGauss(out,[2 2 2]);
           gmask =GeneralAnalysis.imgThreshold(img_g);
           img4lap = GeneralAnalysis.imgGauss(out,[1 1 4]);
           img_laplcutoff = GeneralAnalysis.imgLaplaceCutoff(img4lap);
           lmask =GeneralAnalysis.imgThreshold(img_laplcutoff);
           mask_out = gmask|lmask;
           obj.mask = mask_out;
       end
       function make_thickMask(obj)
           % make sum projection of masked image, then thicken and bridge
           sumproj_out = GeneralAnalysis.sumproj_masktimeseries(obj.mask);
           sumproj_out_thick = GeneralAnalysis.bwmorph_timeseries(sumproj_out,'thicken',2);
           obj.mask_thick = GeneralAnalysis.bwmorph_timeseries(sumproj_out_thick,'bridge');
       end
       function selectSoma(obj)
           uiwait(msgbox('Click >2 points to select a cell soma ROI','Draw Soma ROI','modal'));
           h = dipshow(gaussf(obj.image(:,:,floor(size(obj.image,3)/2)),[1 1 0]),'log');
           diptruesize(h,150);
           [roi, v] = diproi(h);      
           obj.soma_mask = repmat(roi,[1 1 size(obj.image,3)]);
           obj.soma_vertices = v;
           close(h);        
       end
       function selectAIS(obj)
           uiwait(msgbox('Click >2 points to select a cell cell AIS','Draw AIS ROI','modal'));
           h = dipshow(gaussf(obj.image(:,:,floor(size(obj.image,3)/2)),[1 1 0]),'log');
           diptruesize(h,150);
           [roi, v] = diproi(h);      
           obj.AIS_mask = repmat(roi,[1 1 size(obj.image,3)]);
           obj.AIS_vertices = v;
           close(h);        
       end
       function selectFullSoma(obj)
           uiwait(msgbox('Click >2 points to select a cell soma ROI','Draw Soma ROI','modal'));
           h = dipshow(gaussf(obj.image(:,:,floor(size(obj.image,3)/2)),[1 1 0]),'log');
           diptruesize(h,150);
           [roi, v] = diproi(h);      
           obj.fullsoma_mask = repmat(roi,[1 1 size(obj.image,3)]);
           obj.fullsoma_vertices = v;
           close(h);        
       end
       function [h,overlayim] = viewSoma(obj)
           mskcol = [1 0 0];
           cm = bone(256);
           [h,overlayim] = GeneralAnalysis.overlay(obj.image,obj.soma_mask,cm,mskcol);
       end
       
   end
    
end