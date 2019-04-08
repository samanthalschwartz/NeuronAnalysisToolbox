% identify files
filepath = 'Z:\Matt\091018_v2_bassoon_ratio\Dark_V2_basoon5_w2640.tif';
files = uipickfiles('Prompt','Pick Files','FilterSpec',filepath);
% load files
image = loadtiff(files{1});
% max project
maxproj = @(x)(max(image,[],3));
mimage = maxproj(image);
% filter image
gmimage = gaussf(mimage);
lgmimage = GeneralAnalysis.imgLaplaceCutoff(gmimage);
% threshold
masks = threshold(lgmimage,'isodata');
%%

test= GeneralAnalysis.overlay(dip_image(mimage),masks)




%%
