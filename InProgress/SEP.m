classdef SEP < handle
properties
    sep
    cellfill
    msk_cellfill
    max_msk_cellfill
    sep_mask
    sep_mask_max
    sep_mask_fixed
    sep_mask_fixed_label
    sep_sums
    ord_trace
    hm
    oldROIs
    shafttrace
    spinetrace
end

methods (Static)
    function loadimages(filepath)
        image = loadtiff('F:\Hiester et al., 2017 Frontiers Submission\Matt Becker Data (For Review)\SEPGlua1_mch\20150609_SEPGlua1_mch_4_3.tif');% uiopen('E:\Matt Becker Data (For Review)\SEPGlua1_mch\20150609_SEPGlua1_mch_4_2.tif',1);
        obj.sep = image(:,:,:,3);
        obj.cellfill = image(:,:,:,1);
        obj.sepim = dip_image(sep);
    end
    
    function mask_cellfill()
       obj.msk_cellfill = threshold(gaussf(obj.cellfill),'otsu');
        obj.max_msk_cellfill = max(msk_cellfill,[],3);
    end
    
    function mask_sep()
        image_out = opening(obj.sep,5,'elliptic');
        wth = sep - image_out;
        dggsep = GeneralAnalysis.imgDggCutoff(wth);
        varsep = varif(dggsep,5,'elliptic');
        sep_mask = threshold(varsep.^.8,'otsu');
    end
    
    function mask_sep_fixed
        max_sep = max(obj.sep,[],3);
        max_spines = max(obj.mask_sep,[],3);
        wshed = watershed(-gaussf(max_sep),2);
        max_spines(wshed) = 0;
        lb = label(max_spines,1,30);
        obj.sep_mask_fixed = repmat(lb>0,1,1,size(sep,3));
        obj.sep_mask_fixed_label = repmat(lb,1,1,size(sep,3));
    end
    
    function calculate_sepintensities()
        spinelabel = obj.sep_mask_fixed_label;
        sums = zeros(max(spinelabel),size(sep,3));
        wb = waitbar(0);
        sepim = dip_image(sep);
        measure  
        
        
        for ll = 1:max(spinelabel)
            sums(ll,:) = single(sum(sepim,spinelabel == ll,[1 2]));
            waitbar(ll/max(spinelabel),wb);
        end
        close(wb)
        obj.sep_sums = sums;
       
        
    end
    
    
    
    
    
    function img_out = cropimage(img_in,cutoff)
    end
    function viewspines
    end
   
end
end