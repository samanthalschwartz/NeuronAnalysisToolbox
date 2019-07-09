function checkViewInput(TreeNodeText)
cellfill_im = single(aa.cellFill.image);
cellfill_mask = single(aa.cellFill.mask);
surfcarg_im = single(aa.surfaceCargo.image);
surfcarg_mask = single(aa.surfaceCargo.mask);

switch TreeNodeText
    case 'Cell Marker'
        if isempty(aa.cellFill.rawimage)
            f = errordlg(['Need to Load in a Cell Fill image first (Use Image Loader Button)']);
            return;
        end
        if ~isempty(cellfill_mask)
            maskval = prctile(cellfill_im,0.75);
            CellFill = cat(4,cellfill_im,single(cellfill_mask.*maskval));
        else
            CellFill = cellfill_im;
        end
        imp = MatIJ.showImage(CellFill);
        
    case 'Cell Image'
        if isempty(aa.cellFill.rawimage)
            f = errordlg(['Need to Load in a Cell Fill image first (Use Image Loader Button)']);
            return;
        end
        CellFillImage = cellfill_im;
        imp = MatIJ.showImage(CellFillImage);
        
    case 'Cell Mask'
        if isempty(cellfill_mask)
            f = errordlg(['Need to create a Cell Image Mask First (Use Masking Button)']);
            return;
        end
        CellFillMask = cellfill_mask;
        imp = MatIJ.showImage(CellFillMask);
        
    case 'Surface Signal'
        if isempty(aa.surfaceCargo.rawimage)
            f = errordlg(['Need to Load in a Surface Cargo image first (Use Image Loader Button)']);
            return;
        end
        if ~isempty(surfcarg_mask)
            maskval = prctile(surfcarg_im,0.75);
            SurfaceSignal = cat(4,surfcarg_im,single(surfcarg_mask.*maskval));
        else
            SurfaceSignal = surfcarg_im;
        end
        imp = MatIJ.showImage(SurfaceSignal);
        
    case 'Surface Signal Image'
         if isempty(aa.surfaceCargo.rawimage)
            f = errordlg(['Need to Load in a Surface Signal Image first (Use Image Loader Button)']);
            return;
        end
        SurfaceImage = surfcarg_im;
        imp = MatIJ.showImage(SurfaceImage);
        
    case 'Surface Signal Mask'
         if isempty(surfcarg_mask)
            f = errordlg(['Need to create a Surface Signal Mask First (Use Masking Button)']);
            return;
        end
        SurfaceMask = surfcarg_mask;
        imp = MatIJ.showImage(SurfaceMask);
        
    case 'Distance Mask'
    case 'AIS Region'
    case 'Soma Region'
        
    case 'Temporal Heatmap'
        
end

end