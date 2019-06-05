% load images

[filename, pathname] = uigetfile( ...
       {'*.tiff;*.tif;*.TIFF;*.TIF', 'All TIF Files (*..tiff, *.tif, *.TIFF, *.TIF)'; ...
        '*.*',                   'All Files (*.*)'}, ...
        'Load Composite Image');
uiopen(fullfile(pathname,filename)); close all;
%%
[cellfillid, golgiid] =  Mason_selectChannels(image)
cellfill = image(:,:,:,cellfillid);
golgi= image(:,:,:,golgiid);

%%
numslices= size(cellfill,3);
% mask cellfill channel
cellfill_l = GeneralAnalysis.imgLaplaceCutoff(cellfill); %laplace cuttoff filter
[cellfill_lm ,threshval_l,C]= GeneralAnalysis.imgThreshold_fixedUserInput(cellfill_l(:,:,ceil(numslices/2))); %user select
cellfill_g = gaussf(cellfill); % gaussfilter
reg = GeneralAnalysis.crop(cellfill_g,C); %apply gaussfilter cuttoff in region selected above
threshval_g = max(reg);
cellfill_gm = threshold(cellfill_g,'fixed',threshval_g);
cellfillmask = cellfill_lm | cellfill_gm; % or both laplace and gauss masks
cellfillmask = bdilation(cellfillmask,3);
cellfillmask = GeneralAnalysis.bwmorph_timeseries(cellfillmask,'bridge',20);

% make background corrected golgi image
    % - mask golgi (everything)
    % - calculate background filled image 
    % - subtract from golgi image
golgi_g = gaussf(golgi);
[golgi_gm ,threshval_golgi,C_golgi]= GeneralAnalysis.imgThreshold_fixedUserInput(golgi_g(:,:,ceil(numslices/2))); %user select
golgibg = golgi.*~bdilation(golgi_gm,4);
bgim_out = GeneralAnalysis.regionfill_timeseries(golgibg,dip_image(logical(golgi_gm)));
% select soma region
    
% calculate geo distance
sinkdist = bwdistgeodesic(logical(mask_out),logical(soma_mask),'quasi-euclidean');
% make distance bins

% calculate density per bin
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