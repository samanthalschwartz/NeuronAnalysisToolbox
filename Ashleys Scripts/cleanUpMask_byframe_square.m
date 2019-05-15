function newmask = cleanUpMask_byframe_square(underimgin,mask_in,imviewsz)
         if nargin<3
             imviewsz = 150;
         end
        lb = slice_op('label',(logical(mask_in));
        ov = underimgin;
        ov(lb~=0) = 0;
        g = dipfig('ov');
        try
        dipshow(ov,'log');
        catch
            dipshow(ov,'percentile');
        end
        diptruesize(g,imviewsz);
        clmp = bone(255);
        clmp(1,:) = [1 0 0];
        while(ishandle(g))
            try
                [B,C] = dipcrop(g);
            catch
                break;
            end
            gcfinfo = get(g,'UserData');
            if ndims(B)==3
                
                currtime = gcfinfo.curslice;
                 img2remove = lb(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2),currtime);
                 lbs2remove = unique(single(img2remove));
                 
            elseif ismatrix(B)
                 img2remove = lb(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2));
                 lbs2remove = unique(single(img2remove));
            end
            for ii = lbs2remove(lbs2remove~=0)'
            lb(lb == ii) = 0; 
            end
            ov = underimgin;
            ov(lb~=0) = 0;
            diptruesize(gcf,imviewsz);
            try
            dipmapping('log')
            catch
                dipmapping('percentile')
            end
            dipmapping('colormap',clmp);
        end
      dipfig -unlink
      newmask = logical(lb);  
    end