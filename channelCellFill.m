classdef channelCellFill < channelBase
    properties
        mask_thick
        lsig = [1 1 0];
        gsig = [1 1 0];
    end
    
    methods
       function obj = channelCellFill()
           obj = obj@channelBase();
       end
       function mask_img(obj)
        img_m = medif(obj.image,3);
        img_g = GeneralAnalysis.imgGauss(img_m,obj.gsig);
        gmask =GeneralAnalysis.imgThreshold(img_g);
        img_laplcutoff = GeneralAnalysis.imgLaplaceCutoff(img_m,obj.lsig,obj.gsig);
        [lmask,threshval] = GeneralAnalysis.imgThreshold_fixedUserInput(img_laplcutoff);
        mask_out = gmask|lmask;
        obj.mask = mask_out;
        obj.backgroundvalue = threshval;
       end
       function mask_img_other(obj)
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
       
   end
    
end