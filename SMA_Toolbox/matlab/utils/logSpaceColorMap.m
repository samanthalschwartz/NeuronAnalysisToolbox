function logmapfun = logSpaceColorMap(basemap)
    function M = logmap(n)
        bM = flipud(basemap(n));
        M=flipud(bM(flipud(round(logspace(log10(1),log10(n),n))),:));
    end
    logmapfun = @logmap;
end
