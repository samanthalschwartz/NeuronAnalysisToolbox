% load images
startdir = '/Users/masonkleinjan/Desktop/Data analysis/101618_Giantin organelle coloc';
[filename_C] = uipickfiles('FilterSpec',startdir,'Prompt',...
    'Load CellFill Image',...
    'REFilter','*.*');
uiopen(filename_C{1}); close all;
cellfill = image;
[filename] = uipickfiles('FilterSpec',startdir,'Prompt',...
    'Load Golgi Image',...
    'REFilter','*.*');
uiopen(filename{1}); close all;
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
dipshow(cellfillmask,'log');
%%
% make background corrected golgi image
% - mask golgi (everything)
% - calculate background filled image
% - subtract from golgi image
golgi_g = gaussf(golgi);
[golgi_gm ,threshval_golgi,C_golgi]= GeneralAnalysis.imgThreshold_fixedUserInput(golgi_g); %user select
cleaned_golgimask = GeneralAnalysis.cleanUpMask_manual_square(dip_image(golgi),logical(golgi_gm),25);
golgi_gm = bdilation(cleaned_golgimask,4);
golgibg = golgi.*~golgi_gm;
bgim_out = GeneralAnalysis.regionfill_timeseries(golgibg,golgi_gm);
golgi_bgim = gaussf(bgim_out);
golgi_bgsubtract = golgi - golgi_bgim;
cleaned_golgimask = GeneralAnalysis.cleanUpMask_manual_square(dip_image(golgi),logical(golgi_gm),25);


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
savefilename = GeneralAnalysis.filename_addon(filename_C{1},'_distancemask');
GeneralAnalysis.LibTiff(sinkdist,savefilename);

somafilename = GeneralAnalysis.filename_addon(filename_C{1},'_somamask');
GeneralAnalysis.LibTiff(sinkdist,somafilename);

golgibgfilename = GeneralAnalysis.filename_addon(filename{1},'_golgiBgSubtract');
GeneralAnalysis.LibTiff(golgi_bgsubtract,golgibgfilename);

golgimaskfilename = GeneralAnalysis.filename_addon(filename{1},'_golgiMASK');
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
[FILEPATH,newfilename,EXT] = fileparts(filename{1});
xlswrite(fullfile(FILEPATH,newfilename),[edges',density]);
saveas(f,fullfile(FILEPATH,newfilename),'png');
saveas(f,fullfile(FILEPATH,newfilename),'fig');







%%------- function for stitching images ----

%% first load the images

startdir = '/Users/masonkleinjan/Desktop/Data analysis/052819_golgilabeling_SP';

uiwait(msgbox('Select the 2 files you would like to stitch (remember the order)','Title','modal'));

[filenames] = uipickfiles('FilterSpec',startdir,'Prompt',...
    
'Load files to stitch',...
    
'REFilter','*.*');

if numel(filenames)==2
    
    uiopen(filenames{1});
    
    im1 = image;
    
    uiopen(filenames{2});
    
    im2 = image;
    
else
    
    errh = errordlg('You must select exactly 2 files','Error');
    
end

%% then try setting the offset : run this cell until you like the overlay

close all;

ccpeak = {2040,936}; % <------ this is what you change: {up/down,rightleft}

%                                          - {larger moves 'left' image down,larger moves 'left' image right}

im1max = single(max(im1,[],3));

im2max = single(max(im2,[],3));

[stitchimage, ccpeak] = GeneralAnalysis.stitch2images(im1max,im2max,ccpeak,1);

h = dipshow(stitchimage,'log');

%% now stitch this image and select the corresponding other channels you would like also stitch

ccpeaks = cell(size(im1,3),1);

ccpeaks(:) = {ccpeak};

[stitched1, ccpeaks] = GeneralAnalysis.stitch2movies(im1,im2,ccpeaks);

[FILEPATH,saveNAME,EXT] = fileparts(filenames{1});

GeneralAnalysis.LibTiff(stitched1,fullfile(FILEPATH,[saveNAME '_stitched']));

uiwait(msgbox(['Select 2 files from the corresponding channel you would like to stitch '...
    
'(select the file corresponding to ' saveNAME ' first)'],'Title','modal'));

[filenames2] = uipickfiles('FilterSpec',startdir,'Prompt',...
    
'Load files to stitch',...
    
'REFilter','*.*');



if numel(filenames2)==2
    
    uiopen(filenames2{1});
    
    im1 = image;
    
    uiopen(filenames2{2});
    
    im2 = image;
    
else
    
    errh = errordlg('You must select exactly 2 files','Error');
    
end

[stitched2, ccpeaks] = GeneralAnalysis.stitch2movies(im1,im2,ccpeaks);

h = dipshow(stitched2,'log');

[FILEPATH,saveNAME,EXT] = fileparts(filenames2{1});

GeneralAnalysis.LibTiff(stitched2,fullfile(FILEPATH,[saveNAME '_stitched']));

