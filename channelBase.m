classdef channelBase < handle
   properties
      image = [];
      mask = [];
      ROI_trim = 0; % this follows dipcrop output format. column 1 is upper left corner. column 2 is size.
      sigG = 1;
      sigL = 1;
      backgroundvalue = 0;
       bgsub_image = [];
   end
   properties (GetAccess = 'public', SetAccess = 'private')
       filepath = [];
       check_loadstate = false;
       rawimage = [];
   end
   
   methods
       function setfilepath(obj,filepath)
          obj.filepath = filepath;
          obj.check_loadstate = false;
       end
       function setimage(obj,image_in)
           if nargin<2               
               obj.loadimage();
           else
               obj.rawimage = dip_image(image_in);
               obj.check_loadstate = true;
               obj.loadimage();
           end
       end
       function loadimage(obj)
           if ~(obj.check_loadstate == true)
           obj.rawimage  = GeneralAnalysis.loadtiff_1ch(obj.filepath);
           obj.check_loadstate = true;
           end 
           if obj.ROI_trim
           obj.trim_rawimage();
           else
               obj.image = obj.rawimage;
           end 
       end
       function trim_rawimage(obj,a)
           if isempty(obj.ROI_trim) || nargin>1
               h = dipshow(obj.rawimage);
               try
                   dipmapping(h,'log');
               catch
                   dipmapping(h,'percentile');
               end
               diptruesize(h,100);
               answer = questdlg('Select the region to use as your image', ...
                   'Options', ...
                   'OK','Oops!','OK');
               % Handle response
               switch answer
                   case 'OK'
                       [B,C] = dipcrop(h);
                       obj.image = B;
                       obj.ROI_trim = [C,[0;size(obj.rawimage,3)-1]];
                       close(h);
                   case 'Oops!'
                       close(h);
                       return;
               end
           else
               cropimagefromtrimROI(obj)
           end
       end
       function cropimagefromtrimROI(obj)
          obj.image = obj.rawimage(obj.ROI_trim(1,1):obj.ROI_trim(1,1)+obj.ROI_trim(2,1),obj.ROI_trim(1,2):obj.ROI_trim(1,2)+obj.ROI_trim(2,2),obj.ROI_trim(1,3):obj.ROI_trim(1,3)+obj.ROI_trim(2,3));   
       end
       function resetImage(obj)
          obj.image = obj.rawimage; 
          obj.ROI_trim = [];
          obj.mask = [];
          obj.bgsub_image = [];
       end
       function clearImage(obj)
          obj.resetImage;
          obj.rawimage = [];
       end
       function clipmasks(obj)
           clip = 20;
           goodmask = obj.mask.*0';
           goodmask(clip:(end-clip),clip:(end-clip),:) = 1;
           obj.mask = obj.mask.*goodmask;
       end
       
       function removeBoundaryMaskArtifact(obj)
           a = sum(obj.mask,[],3);
           g = dipfig('a');
           while(ishandle(g))
               try
                   [B,C] = dipcrop(g);
               catch
                   break;
               end
               gcfinfo = get(g,'UserData');
               a(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2)) = 0;
           end
           dipfig -unlink
%            newmask = logical(lb);
       end
    
     
       
       
       
       function [h,overlayim] = viewMaskOverlayPerim(obj,cm,mskcol)
           if nargin<3
               mskcol = [1 1 1];
           end
           if nargin<2
               cm = hot(256);
           end
           perim = dt(obj.mask);
           bin_im = (perim==1);
           [h,overlayim] = GeneralAnalysis.overlay(obj.image,bin_im,cm,mskcol);
       end
       function [h,overlayim] = viewMaskOverlayPerimstatic(image,mask,cm,mskcol)
           if nargin<3
               mskcol = [1 1 1];
           end
           if nargin<2
               cm = hot(256);
           end
           perim = dt(mask);
           bin_im = (perim==1);
           [h,overlayim] = GeneralAnalysis.overlay(image,bin_im,cm,mskcol);
       end
       function [h,overlayim] = viewMaskOverlayFill(obj,cm,mskcol)
           if nargin<3
               mskcol = [1 0 0];
           end
           if nargin<2
               cm = bone(256);
           end
           cm(1,:) = mskcol;
           [h,overlayim] = GeneralAnalysis.overlay(obj.image,obj.mask,cm,mskcol);
       end
   end
   methods (Static)
       
   end
end