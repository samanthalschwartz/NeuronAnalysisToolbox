function imagelist = updateImageList()
        try
            IJfilenames = cell(MIJ.getListImages);
        catch
            IJfilenames = [];
        end
        figHandles = findobj('Type', 'image');
        Mfignames = arrayfun(@(x) ['MatFig Window: ' num2str(x.Parent.Parent.Number)],figHandles,'UniformOutput',false);
        imagelist.numIJfiles = numel(IJfilenames);
        imagelist.numMATfiles = numel(Mfignames);
        imagelist.imagehandles = [IJfilenames;arrayfun(@(x) {x}, figHandles)];
        imagelist.list = [IJfilenames;Mfignames];
        if ~isempty(imagelist.list)
            imagelist.list = ['  ';imagelist.list];
        else
            imagelist.list = {'  ';'Nothing Open'};
        end
    end
    function allvars = getMatvars()
        allvars = arrayfun(@(x) x.name,whos,'UniformOutput',false);
        for ii = 1:numel(allvars)
            if ndims(allvars{ii})<2 %only include possible images
                allvars{ii} = [];
            end
        end
        if ~isempty(allvars)
            allvars = ['  ';allvars];
        else
            allvars = {'  ';'Nothing Open'};
        end
    end
    