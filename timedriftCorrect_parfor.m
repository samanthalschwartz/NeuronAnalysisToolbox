function [img_out,sv_arr] = timedriftCorrect_parfor(img_in,shiftmeth)
         if nargin < 2
             shiftmeth = 'iter';
         end
         if ~isa(img_in,'dip_image')
            try
                img_in = dip_image(img_in);
            catch
                warning('input must be an image matrix');
                    return;
            end
         end
        imref = squeeze(img_in(:,:,0));
        sv_arr = nan(2,size(img_in,3)-1); % array of shift vectors length is one less than size(im_in,3)
        tsize = size(img_in,3);
        tic
        img_in = single(img_in); %make not a dip_image for parfor
        shiftim = zeros(size(img_in)); %create empty array for shift images
%         shiftim(:,:,1) = img_in(:,:,1);
%         shiftim = {};
        parfor (ii = 2:tsize,4)
            sv1 = findshift(imref,img_in(:,:,ii),shiftmeth,0);
            sv_arr(:,ii-1) = sv1;
            currim = shift(img_in(:,:,ii),sv1,1);
            shiftim(:,:,ii) = single(currim);
        end
        toc
        tic
        img_out = 0.*shiftim;
        for ii = 1:(tsize-1)
        sv1 = sv_arr(:,ii);
        currim = shiftim(:,:,ii+1);
            if size(sv1,1) == 2
               if sv1(1)>0
                   currim(:,1:ceil(sv1(1))) = 0;
               else
                   currim(:,end+ceil(sv1(1)):end) = 0;
               end
               if sv1(2)>0
                   currim(1:ceil(sv1(2)),:) = 0;
               else
                   currim(end+ceil(sv1(2)):end,:) = 0;
               end
            end
            img_out(:,:,ii+1) = currim;
        end
        img_out(:,:,1) = single(img_in(:,:,1));
        toc
%         img_out = test;
end