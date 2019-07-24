%%
zscale = 4.2
4.2*50
%% load images
startdir = 'C:\Users\sammy\Dropbox\Sam Kennedy Lab\Mason';

[filename, pathname] = uigetfile( ...
       {'*.tiff;*.tif;*.TIFF;*.TIF', 'All TIF Files (*..tiff, *.tif, *.TIFF, *.TIF)'; ...
        '*.*',                   'All Files (*.*)'}, ...
        'Load CellFill Image',startdir);
uiopen(fullfile(pathname,filename)); close all;
cellfill = image;
[filename, pathname] = uigetfile( ...
       {'*.tiff;*.tif;*.TIFF;*.TIF', 'All TIF Files (*..tiff, *.tif, *.TIFF, *.TIF)'; ...
        '*.*',                   'All Files (*.*)'}, ...
        'Load Golgi Image',startdir);
uiopen(fullfile(pathname,filename)); close all;
golgi = image;

%%
numslices= size(cellfill,3);
% mask cellfill channel
cellfill_l = GeneralAnalysis.imgLaplaceCutoff(cellfill); %laplace cuttoff filter
[cellfill_lm ,threshval_l,C]= GeneralAnalysis.imgThreshold_fixedUserInput(cellfill_l(:,:,ceil(numslices/2))); %user select
cellfill_g = gaussf(cellfill); % gaussfilter
reg = GeneralAnalysis.crop(cellfill_g,C); %apply gaussfilter cuttoff in region selected above
threshval_g = max(reg);
%%
cellfill_gm = threshold(cellfill_g,'fixed',threshval_g);
cellfillmask = cellfill_lm | cellfill_gm; % or both laplace and gauss masks
% this is where you can increase the dilation and bridge parameter to make
% sure mask is connected
cellfillmask = bdilation(cellfillmask,5);
cellfillmask = GeneralAnalysis.bwmorph_timeseries(cellfillmask,'bridge',50);
dipshow(cellfillmask,'log')
%%
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
uiwait(msgbox(['Select a region to represent the soma'],'Title','modal'));

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

% save the distance mask, soma mask, golgi mask
savefilename = GeneralAnalysis.filename_addon(fullfile(pathname,filename),'_distancemask');
GeneralAnalysis.LibTiff(sinkdist,savefilename);

somafilename = GeneralAnalysis.filename_addon(fullfile(pathname,filename),'_somamask');
GeneralAnalysis.LibTiff(sinkdist,somafilename);

golgibgfilename = GeneralAnalysis.filename_addon(fullfile(pathname,filename),'_golgiBgSubtract');
GeneralAnalysis.LibTiff(golgi_bgsubtract,golgibgfilename);

golgimaskfilename = GeneralAnalysis.filename_addon(fullfile(pathname,filename),'_golgiMASK');
GeneralAnalysis.LibTiff(golgi_gm,golgimaskfilename);


% make distance bins

% calculate density per bin
[density,edges] = GeneralAnalysis.calculate_DensityPerDistace(golgi_bgsubtract,sinkdist);
f = figure; axe= axes(f);
plot(edges',density,'Linewidth',2);
axe.FontSize = 16;
xlabel('Distance in pixels','FontSize',16);
ylabel('Density (AUs)','FontSize',16);

% save excel
[FILEPATH,newfilename,EXT] = fileparts(filename);
xlswrite(fullfile(FILEPATH,newfilename),[edges',density]);
