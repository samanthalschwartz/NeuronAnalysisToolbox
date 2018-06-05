% Mark J. Olah (mjo@cs.unm.edu)
% 04/2015

function outpath=collapsepath(inpath)
    % Removes internal . and .. paths appropriately
    % [IN] 
    % inpath - string: reltaive or absolute path to a file
    % [OUT]
    % outpath - string: modified path with the internal charactors removed
    if isempty(inpath)
        outpath=[];
        return
    end
    indirs = strsplit(inpath,filesep());
    indirs(cellfun(@(d) strcmp(d,'.') | isempty(d), indirs)) = [];
    %remove and save leading ..
    notdoubledots = find(~strcmp(indirs,'..'),1,'first');
    initialdoubledots={};
    if isempty(notdoubledots)
        %all dirs are doubledots
        if inpath(1)==filesep()
            error('collapsepath:input','bad path structure');
        end
        outpath = fullfile(indirs{:});
        return;
    elseif notdoubledots>1
        if inpath(1)==filesep()
            error('collapsepath:input','bad path structure');
        end
        %first several are double dots
        initialdoubledots=indirs(1:notdoubledots-1);
        indirs=indirs(notdoubledots:end);
    end
    doubledots = find(strcmp(indirs,'..'),1,'first');
    while ~isempty(doubledots)
        if doubledots==1
            indirs(doubledots:doubledots) = [];
        else
            indirs(doubledots-1:doubledots) = [];
        end
        doubledots = find(strcmp(indirs,'..'),1,'first');
    end
    if isempty(indirs)
        outpath='.';
    elseif isunix() && inpath(1)==filesep()        
        outpath = fullfile(filesep(),indirs{:});
    else
        outpath = fullfile(initialdoubledots{:},indirs{:});
    end
end
