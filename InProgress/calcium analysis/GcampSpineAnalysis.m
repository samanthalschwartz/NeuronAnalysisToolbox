%% load in files and all that (cell fill and align etc)
ga = GeneralAnalysis(); %make just an empty version of the class to be able to use ga as a shortcut
% -- get pre file
[FILENAME_pre, PATHNAME_pre] = uigetfile(fullfile(pwd,'*.*'),'Select a pre condition calcium image');
prompt_pre = 'Select the Unique File Identifier: use * for wildcard';
name = 'Calcium Data Selection';
defaultanswer = {FILENAME_pre};
numlines = 1;
flstring_pre=inputdlg(prompt_pre,name,numlines,defaultanswer);
flstring_pre = flstring_pre{1};
% -- get post file
[FILENAME_post, PATHNAME_post] = uigetfile(fullfile(pwd,'*.*'),'Select a post condition calcium image');
prompt_post = 'Select the Unique File Identifier: use * for wildcard';
name = 'Calcium Data Selection';
defaultanswer = {FILENAME_post};
numlines = 1;
flstring_post=inputdlg(prompt_post,name,numlines,defaultanswer);
flstring_post = flstring_post{1};
% -- get cell fill file
[FILENAME_cf, PATHNAME_cf] = uigetfile(fullfile(PATHNAME,'*.*'),'Select the corresponding cell fill image');

img_pre = dip_image(ga.loadtiffseries(PATHNAME_pre,flstring_pre));
img_post = dip_image(ga.loadtiffseries(PATHNAME_post,flstring_post));
cellfill_pre = ga.loadtiff_1ch(fullfile(PATHNAME_cf,FILENAME_cf));

% align post image with last frame of pre-image
pre_lastframe = img_pre(:,:,end);
post_firstframe = img_post(:,:,end);
img_post_foralign = cat(3,pre_lastframe,post_firstframe);
[out,sv_arr] = ga.timedriftCorrect(img_post_foralign);
test = ga.applyshift2series(img_post,sv_arr);
test = ga.applydriftCorrect(post_firstframe,sv_arr);
test1 = cat(3,pre_lastframe,test);
img_post_align = ga.applydriftCorrect(img_post,repmat(sv_arr,1,size(img_post,3)));
img_post = img_post_aligned(:,:,1:end);

test = cat(3,img_pre,img_post);
% align cellfill to img
imgsum = sum(img,[],3);
corrimg = cat(3,imgsum,cellfill_pre);
dc_corrimg = ga.timedriftCorrect(corrimg);
cellfill = dc_corrimg(:,:,end);

%% try and find the calcium transients
% if last frame is empty uncomment this line
threshold_parameter = 1.5;
ca = img;
ca = ca(:,:,0:end-1);
% first smooth image
ca_g = gaussf(ca,[1 1 0]);
% now normalize to median
ca_gnorm = ca_g;%- medif(ca_g);
% now take difference
ca_dzz = dzz(ca_gnorm,[1 1 1]);
% get absolute values
ab_ca_dzz = abs(ca_dzz);
gab_ca_dzz = gaussf(ab_ca_dzz,[1 1 0]);
tt = threshold(ab_ca_dzz^threshold_parameter,'otsu');
lbl = label(tt,1,1,10^10);
% this is to view the identified masks over the image
[h,overlayim] = GeneralAnalysis.viewMaskOverlayPerimStatic(ca,bdilation(lbl>0,2));

%% try another way
thresh_param = 1;
imgg = gaussf(ca,[1 1 0]);
diffimage = dzz(imgg);
adiff = abs(diffimage);
lp = ga.imgLaplaceCutoff(adiff,[1 1 0],[1 1 0]); %this helps smooth image (especially enhancing round objects)
lp_ = lp.^thresh_param;
tt = threshold(lp_,'otsu');
lb = label(tt);
[h,overlayim] = GeneralAnalysis.viewMaskOverlayPerimStatic(ca,bdilation(lb>0,2));