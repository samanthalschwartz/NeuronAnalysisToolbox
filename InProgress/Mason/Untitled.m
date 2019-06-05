fold = 'C:\Users\schwsama\Dropbox\Sam Kennedy Lab\Mason';
uiopen('C:\Users\schwsama\Dropbox\Sam Kennedy Lab\Mason\mjk59_mch_sp488_gm130640_1_7soma_w1561.TIF',1)
soma = image;
uiopen('C:\Users\schwsama\Dropbox\Sam Kennedy Lab\Mason\mjk59_mch_sp488_gm130640_1_7den_w1561.TIF',1)
den = image;

close all
ccpeak = {2055,940};
[stitchim, ccpeak]= GeneralAnalysis.stitch2images(maxden,maxsoma,ccpeak);
dipshow(stitchim,'log');
ccpeaks = cell(size(den,3),1);
ccpeaks(:) = {ccpeak};
[stitchstack, ccpeaks]= GeneralAnalysis.stitch2movies(den,soma,ccpeaks);
dipshow(stitchstack,'log');
LibTiff(fullfile(fold,'mjk59_mch_sp488_gm130640_1_w1561'),'stitchstack')

%%
uiopen('C:\Users\schwsama\Dropbox\Sam Kennedy Lab\Mason\mjk59_mch_sp488_gm130640_1_7soma_w3640.TIF',1)
soma = image;

uiopen('C:\Users\schwsama\Dropbox\Sam Kennedy Lab\Mason\mjk59_mch_sp488_gm130640_1_7den_w3640.TIF',1)
den = image;

[stitchstack, ccpeaks]= GeneralAnalysis.stitch2movies(den,soma,ccpeaks);
dipshow(stitchstack,'log');
GeneralAnalysis.LibTiff(stitchstack,fullfile(fold,'mjk59_mch_sp488_gm130640_1_w3640'))
%%
%make bins for calculating densities

h = histogram(sinkdist,100);
edges = h.BinEdges;
density = zeros(size(edges,2),1);

wb = waitbar(0);
for ii = 1:size(edges,2)

    testmask = sinkdist<edges(ii+1) & sinkdist>=edges(ii);
    inmask = testmask.*corr_golgi;
    ints = sum(inmask(:));
    numpixels = sum(testmask(:));
    density(ii) = ints/numpixels;
    waitbar(ii/size(edges,2),wb);
end
close(wb)
figure; plot(edges',density)