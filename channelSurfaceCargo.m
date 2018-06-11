classdef channelSurfaceCargo < channelBase
   properties
       setpt_highsensitivity = 1;
       lsig = [1 1 0];
       gsig = [1 1 0];
   end
   
   methods
       function obj = channelSurfaceCargo()
           obj = obj@channelBase();
       end
       function mask_img(obj)
          if obj.setpt_highsensitivity
              obj.mask_img_highsens(); 
          else
              obj.mask_img_lowsens();
          end
              
       end
       function mask_img_lowsens(obj,lsig,gsig)
           if nargin<3
               gsig = [1 1 3];
           end
           if nargin<2
               lsig = [1 1 0];
           end
           img_m = medif(obj.image,4);
           img_laplcutoff = GeneralAnalysis.imgLaplaceCutoff(obj.image,lsig,gsig);
           normim = img_laplcutoff./max(img_laplcutoff);
           obj.mask = GeneralAnalysis.imgThreshold(normim*img_m);
       end
       function mask_img_highsens(obj,lsig,gsig)
           if nargin<3
               obj.gsig = [1 1 1];
           end
           if nargin<2
               obj.lsig = [1 1 0];
           end
           img_m = medif(obj.image,3);
           img_mg = gaussf(img_m,obj.gsig);
           img_laplcutoff = GeneralAnalysis.imgLaplaceCutoff(img_mg,obj.lsig,obj.gsig);
           %            glim = gaussf(img_laplcutoff);
           %            obj.mask = GeneralAnalysis.imgThreshold(img_laplcutoff);
           [maskpre,threshval] = GeneralAnalysis.imgThreshold_fixedUserInput(img_laplcutoff);
           %            maskpre = GeneralAnalysis.imgThreshold(glim.^1.2);
                     
           mask = GeneralAnalysis.bwmorph_timeseries(maskpre,'thicken',1);
           ll = slice_op('watershed',-img_mg,2);
           maskws = mask;
           maskws(ll) = 0;
           
           obj.mask = logical(maskws);
%            obj.mask = logical(mask);
           obj.backgroundvalue = threshval;
       end
   end
end