topdir = 'C:\Users\sammy\Dropbox\Sam Kennedy Lab\Emma';
basename = 'Slice 2 RS 1_XY1556137811_Z00_T0_';
ext = '.tif';
ch_cellfillstr = 'C2';
ch_tag1 = 'C0';
ch_tag2 = 'C1';

ch_cellfillpath = fullfile(topdir,[basename ch_cellfillstr ext]);
ch_tag1path = fullfile(topdir,[basename ch_tag1 ext]);
ch_tag2path = fullfile(topdir,[basename ch_tag2 ext]);
uiopen(ch_cellfillpath); cellfill = image;
uiopen(ch_tag1path); tag1 = image;
uiopen(ch_tag2path); tag2 = image;
joinchannels('rgb',stretch(tag1),stretch(tag2),stretch(cellfill))

numslices= size(cellfill,3);
T = adaptthresh(cellfill,0.1,'ForegroundPolarity','dark');
cellfillb = single(cellfill) - 4*single(T);
cellfill_l = GeneralAnalysis.imgLaplaceCutoff(cellfillb); %laplace cuttoff filter
cellfill_g = gaussf(cellfillb); % gaussfilter

% mask cellfill channel
[~ ,threshval_l,C]= GeneralAnalysis.imgThreshold_fixedUserInput(cellfill_l(:,:,ceil(numslices/2))); %user select
cellfill_lm = cellfill_l>threshval_l;


 
cellmask = cellfill_lm | cellfill_gm;

% newmask = GeneralAnalysis.cleanUpMaskSuper_byframe_square(dip_image(image),cellmask,100);

h = dipshow(max(cellfill,[],3),'lin')
diptruesize(h,50);
[roi, v] = diproi(h);
s_mask = repmat(roi,[1 1 size(cellfill,3)]);
soma_vertices = v;
close(h);
soma_mask = s_mask.*cellmask;

