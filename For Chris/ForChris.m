classdef ForChris < handle
properties
    test = [];
end

methods (Static)
    function [h,overlayarr] = viewMaskOverlay(grayim,mask)
         if ~isa(grayim,'dip_image')
            try
                grayim = dip_image(grayim);
            catch
                warning('input must be an image matrix');
                    return;
            end
        end
         assert(ndims(grayim) == ndims(mask));
         grayim_minusmask = grayim.*~mask;
         mskfrm = max(grayim)*10*mask + grayim_minusmask;
         switch ndims(grayim)
             case 2
                 rch = cat(3,mskfrm,grayim);
                 gch = cat(3,grayim_minusmask,grayim);
             case 3
                 rch = cat(4,mskfrm,grayim);
                 gch = cat(4,grayim_minusmask,grayim);
         end
         bch = gch;
         overlayarr = joinchannels('rgb',rch,gch,bch); 
         h = dipshow(overlayarr,'log');
    end
      function [mask,threshval] = imgThreshold_fixedUserInput(img_in,image4selection)
        if ~isa(img_in,'dip_image')
            img_in = dip_image(img_in);
        end
        uiwait(msgbox('Select a representative background region','Title','modal'));
        if nargin<2
            image4selection = img_in;
        end
        h = dipshow(image4selection,'log');
        diptruesize(h,125);
        [~,C] = dipcrop(h);
        if ndims(img_in)==3
        reg = img_in(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2),:);
        elseif ismatrix(img_in)
            reg = img_in(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2));
        end
        threshval = max(reg);
        mask = threshold(img_in,'fixed',threshval);
        close(h);
      end
      function [] = LibTiff(Vol,inputname)
          %LibTiff wrote by Michael @ mic.muenter@uni-luebeck.de
          % Input:uint8 Volume (x,y,z)
          % INPUT: filepath
          
          % INPUT example: Vol = 255.*ones(300,100,200); LibTiff(Vol);
          
          Vol = uint32(Vol);
          
          t = Tiff([inputname,'.tiff'],'w'); % Filename by variable name
          tagstruct.ImageLength = size(Vol,1); % image height
          tagstruct.ImageWidth = size(Vol,2); % image width
          tagstruct.Photometric = Tiff.Photometric.MinIsBlack; % https://de.mathworks.com/help/matlab/ref/tiff.html
          tagstruct.BitsPerSample = 32;
          tagstruct.SamplesPerPixel = 1;
          tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
          tagstruct.Software = 'MATLAB';
          
          tic
          setTag(t,tagstruct)
          write(t,squeeze(Vol(:,:,1)));
          for i=2:size(Vol,3) % Write image data to the file
              writeDirectory(t);
              setTag(t,tagstruct)
              write(t,squeeze(Vol(:,:,i))); % Append
          end
          close(t);
          toc
          
          disp('Tiff File saved')
      end
end
end