% load images

[filename, pathname] = uigetfile( ...
       {'*.tiff;*.tif;*.TIFF;*.TIF', 'All TIF Files (*..tiff, *.tif, *.TIFF, *.TIF)'; ...
        '*.*',                   'All Files (*.*)'}, ...
        'Load Composite Image');
uiopen(fullfile(pathname,filename)); close all;
%%
[cellfillid, golgiid] =  Mason_selectChannels(image);
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
[golgi_gm ,threshval_golgi,C_golgi]= GeneralAnalysis.imgThreshold_fixedUserInput(golgi_g); %user select
golgi_gm = bdilation(golgi_gm,4);
golgibg = golgi.*~golgi_gm;
bgim_out = GeneralAnalysis.regionfill_timeseries(golgibg,golgi_gm);
golgi_bgim = gaussf(bgim_out);
golgi_bgsubtract = golgi - golgi_bgim;

% select soma region
h = dipshow(cellfill(:,:,ceil(numslices/2)))
diptruesize(h,50);
[roi, v] = diproi(h);
s_mask = repmat(roi,[1 1 size(cellfill,3)]);
soma_vertices = v;
close(h);
soma_mask = s_mask.*cellfillmask;
% calculate geo distance
sinkdist = bwdistgeodesic(logical(cellfillmask),logical(soma_mask),'quasi-euclidean');
h = dipshow(sinkdist,jetblack,[0 max(sinkdist(sinkdist~=Inf))*4]);
diptruesize(h,50);

% save the distance mask
savefilename = GeneralAnalysis.filename_addon(fullfile(pathname,filename),'_distancemask');
GeneralAnalysis.LibTiff(sinkdist,savefilename);
% make distance bins

% calculate density per bin
[density,edges] = GeneralAnalysis.calculate_DensityPerDistace(golgi_bgsubtract,sinkdist);
f = figure; axe= axes(f);
plot(edges',density,'Linewidth',2);
axe.FontSize = 16;
xlabel('Distance in pixels','FontSize',16);
ylabel('Density (AUs)','FontSize',16);

savefilename = GeneralAnalysis.filename_addon(fullfile(pathname,filename),'_distanceValues');
xlswrite(savefilename,[edges',density]);
