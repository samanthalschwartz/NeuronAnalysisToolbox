%%------- function for stitching images ----
%% first load the images
startdir = 'C:\Users\sammy\Dropbox\Sam Kennedy Lab\Mason';
uiwait(msgbox('Select the 2 files you would like to stitch (remember the order)','Title','modal'));
[filename, pathname] = uigetfile( ...
    {'*.tiff;*.tif;*.TIFF;*.TIF', 'All TIF Files (*..tiff, *.tif, *.TIFF, *.TIF)'; ...
    '*.*',                   'All Files (*.*)'}, ...
    'Load files to stitch',startdir,'Multiselect','on');
if numel(filename)==2
    
    uiopen(fullfile(pathname,filename{1}));
    im1 = image;
    uiopen(fullfile(pathname,filename{2}));
    im2 = image;
else
    errh = errordlg('You must select exactly 2 files','Error');
end
%% then try setting the offset : run this cell until you like the overlay
close all;
ccpeak = {2060,950}; % <------ this is what you change: {up/down,rightleft} 
%                                          - {larger moves 'left' image down,larger moves 'left' image right}
im1max = single(max(im1,[],3));
im2max = single(max(im2,[],3));
[stitchimage, ccpeak] = GeneralAnalysis.stitch2images(im1max,im2max,ccpeak,1);
h = dipshow(stitchimage,'log');
%% now stitch this image and select the corresponding other channels you would like also stitch
ccpeaks = cell(size(im1,3),1);
ccpeaks(:) = {ccpeak};
[stitched1, ccpeaks] = GeneralAnalysis.stitch2movies(im1,im2,ccpeaks);
[FILEPATH,saveNAME,EXT] = fileparts(filename{1});
GeneralAnalysis.LibTiff(stitched1,fullfile(pathname,[saveNAME '_stitched']));
uiwait(msgbox(['Select 2 files from the corresponding channel you would like to stitch '...
    '(select the file corresponding to ' saveNAME ' first)'],'Title','modal'));
[filename, pathname] = uigetfile( ...
    {'*.tiff;*.tif;*.TIFF;*.TIF', 'All TIF Files (*..tiff, *.tif, *.TIFF, *.TIF)'; ...
    '*.*',                   'All Files (*.*)'}, ...
    ['Load files to stitch, (select the file corresponding to ' saveNAME ' first)'],startdir,'Multiselect','on');

if numel(filename)==2 
    uiopen(fullfile(pathname,filename{1}));
    im1 = image;
    uiopen(fullfile(pathname,filename{2}));
    im2 = image;
else
    errh = errordlg('You must select exactly 2 files','Error');
end
[stitched2, ccpeaks] = GeneralAnalysis.stitch2movies(im1,im2,ccpeaks);
h = dipshow(stitched2,'log');
[FILEPATH,saveNAME,EXT] = fileparts(filename{1});
GeneralAnalysis.LibTiff(stitched2,fullfile(pathname,[saveNAME '_stitched']));



