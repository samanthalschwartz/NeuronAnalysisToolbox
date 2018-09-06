classdef SEP < handle
properties
    % file info
    filepath
    filename
    % images
    sep
    cellfill
    % image masks
    msk_cellfill
    max_msk_cellfill
    sep_mask
    sep_mask_fixed
    % labeled mask
    sep_mask_fixed_label
    % labeled measured info
    sep_sums
    sep_sizes
    % full image measurements
    shafttrace_intensity_varsize
    shafttrace_size_varsize
    shafttrace_intensity_fixedsize
    shafttrace_size_fixedsize
    spinetrace_intensity_varsize
    spinetrace_size_varsize
    spinetrace_intensity_fixedsize
    spinetrace_size_fixedsize
    % heatmap info
    ord_trace
    hm
    % for viewing
    oldROIs
end

methods (Static)
    function loadimages(filepath)
        image = loadtiff(fullfile(filepath,filename));% uiopen('E:\Matt Becker Data (For Review)\SEPGlua1_mch\20150609_SEPGlua1_mch_4_2.tif',1);
        obj.sep = image(:,:,:,3);
        obj.cellfill = image(:,:,:,1);
    end
    
    function make_mask_cellfill(obj)
       obj.msk_cellfill = threshold(gaussf(obj.cellfill),'otsu');
        obj.max_msk_cellfill = max(msk_cellfill,[],3);
    end
    
    function make_mask_sep(obj)
        image_out = opening(obj.sep,5,'elliptic');
        wth = sep - image_out;
        dggsep = GeneralAnalysis.imgDggCutoff(wth);
        varsep = varif(dggsep,5,'elliptic');
        obj.sep_mask = threshold(varsep.^.8,'otsu');
    end
    
    function make_mask_sep_fixed(obj)
        max_sep = max(obj.sep,[],3);
        max_spines = max(obj.mask_sep,[],3);
        wshed = watershed(-gaussf(max_sep),2);
        max_spines(wshed) = 0;
        lb = label(max_spines,1,30);
        obj.sep_mask_fixed = repmat(lb>0,1,1,size(sep,3));
        obj.sep_mask_fixed_label = repmat(lb,1,1,size(sep,3));
    end
    
    function calculate_sepintensities(obj)
        spinelabel = obj.sep_mask_fixed_label;
        sepim = dip_image(obj.sep);
        maxframe = size(obj.sep,3);
        msrarray = cell(1,maxframe);
        wb = waitbar(0);
        for ll = 0:(maxframe-1)
            msr = measure(spinelabel(:,:,ll),sepim(:,:,ll),{'size','sum'});
            msrarray{ll+1} = msr;
            waitbar((ll+1)/maxframe,wb);
        end
        close(wb)
        obj.sep_sums = zeros(max(spinelabel),maxframe);
        obj.sep_sizes = zeros(max(spinelabel),1);
        %reshape array
        for ll = 1:max(spinelabel)
            for tt = 1:maxframe
                obj.sep_sums(ll,tt) = msrarray{tt}(ll).sum;
            end
            obj.sep_sizes(ll) = msrarray{1}(ll).size;
        end
    end
    
    function [hm, f] = calculate_sepIHeatMap(obj)
        f = figure;
        trace = obj.sep_sums./obj.sep_sums(:,1);
        allsum = sum(trace,2);
        [~, ordx] = sort(allsum, 'descend');
        ord_trace = trace(ordx,:);
        obj.hm = heatmap(ord_trace);
        obj.hm.GridVisible = 'off';
        obj.hm.Colormap = jet(50);
        obj.hm.ColorLimits = [0 4];
    end
    
    function loadIJROIs1(obj)
        strFilename = fullfile(filepath,[filename(1:end-4) '_ROI1.zip']);
        obj.oldROIs1 = ReadImageJROI(strFilename);
    end
    
    function loadIJROIs2(obj)
        strFilename = fullfile(filepath,[filename(1:end-4) '_ROI2.zip']);
        obj.oldROIs2 = ReadImageJROI(strFilename);
    end
    
    function plotoldROIs(obj)
        sROI = obj.oldROIs1;
        for ii = 1:numel(sROI)
            radii = (sROI{ii}.vnRectBounds(3)-sROI{ii}.vnRectBounds(1))/2;
            centers = [sROI{ii}.vnRectBounds(2) + radii, sROI{ii}.vnRectBounds(1) + radii];
            viscircles(centers,radii,'Color','g','LineWidth',1.5);
        end
        sROI = obj.oldROIs2;
        for ii = 1:numel(sROI)
            radii = (sROI{ii}.vnRectBounds(3)-sROI{ii}.vnRectBounds(1))/2;
            centers = [sROI{ii}.vnRectBounds(2) + radii, sROI{ii}.vnRectBounds(1) + radii];
            viscircles(centers,radii,'Color',[1 0 1],'LineWidth',1);
        end
    end
        
    function calculateSpineShaftIntensities(obj)
        sepim = dip_image(obj.sep);
        % ---- fixed mask sizes ----
        % create mask of shaft - area of spines/sep intensity
        shaftmask = obj.msk_cellfill.*(-obj.sep_mask_fixed);
        obj.shafttrace_intensity_fixedsize = sum(sepim,shaftmask,[1 2]);
        obj.shafttrace_size_fixedsize = sum(shaftmask(:,:,1),[],[1 2]);
        % create mask of spines/sep intensity
        obj.spinetrace_intensity_fixedsize = sum(sepim,obj.sep_mask_fixed,[1 2]);
        obj.spinetrace_size_fixedsize = sum(obj.sep_mask_fixed(:,:,1),[],[1 2]);
        
               
        % ---- variable mask sizes ----
        % create mask of shaft - area of spines/sep intensity
        % create mask of spines/sep intensity
        clear shaftmask;
        maskedspines = bdilation(obj.sep_mask,1);
        shaftmask = msk_cellfill.*(-maskedspines);
        obj.shafttrace_intensity_varsize = sum(sepim,shaftmask,[1 2]);
        obj.shafttrace_size_varsize = sum(shaftmask,[],[1 2]);
        
        obj.spinetrace_intensity_varsize = sum(sepim,maskedspines,[1 2]);
        obj.spinetrace_size_varsize = sum(maskedspines,[],[1 2]);
    end
    
    function saveSEP(obj,savepath)
        if nargin<2
            savepath = obj.filepath;
        end
        save(fullfile(savepath,obj.filename),'obj');
    end
    
    function h = viewSEPMaskFixed(obj)
        h = ga.viewMaskOverlayPerimStatic(dip_image(obj.sep),obj.sep_mask_fixed);
        dipmapping(h,'colormap',jet);
        dipmapping(h,[0 2500]);
        diptruesize(h,100);
    end
    function h = viewSEPMask(obj)
         h = ga.viewMaskOverlayPerimStatic(dip_image(obj.sep),obj.sep_mask);
        dipmapping(h,'colormap',jet);
        dipmapping(h,[0 2500]);
        diptruesize(h,100);   
    end
    function h = viewSEPMaskFixedwOldRois(obj)
        h = ga.viewMaskOverlayPerimStatic(dip_image(obj.sep),obj.sep_mask_fixed);
        dipmapping(h,'colormap',jet);
        dipmapping(h,[0 2500]);
        diptruesize(h,100);
        plotoldROIs(obj); 
    end
    function h = viewHeatMap(obj)
       h = figure;
       obj.hm;
    end
end
end