% load('C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\111617\psd95ibGFP_NR1561_Abeta647_003_Reconstructed_SIM.mat')
clear all
startdir  = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files';
[FILENAME, PATHNAME, FILTERINDEX] = uigetfile(startdir);
% load in image and look at image and mask
% 'n' and 'p' for next and previous frame
% 'f' and 'b' for toggling the mask on and off
load(fullfile(PATHNAME,FILENAME));
% just re-run abeta masking step in case the new threshold wasn't properly
% run
tic
obj.make_maskchAB;
toc
% look at ch1 overlay to find good planes
GeneralAnalysis.viewMaskOverlay(obj.ch1.image,obj.ch1.mask);
% look at ch2 overlay to find good planes
%GeneralAnalysis.viewMaskOverlay(obj.ch2.image,obj.ch2.mask);
GeneralAnalysis.viewMaskOverlay(obj.abeta.image,obj.abeta.mask);
%% set the plane values 
opts.WindowStyle='normal';
prompt = {'Bottom Plane','Top Plane'};
title = 'select the planes';
dims = [1 35];
definput = {'',''};
answer = inputdlg(prompt,title,dims,definput,opts);
obj.planeBOTTOM = answer{1};
obj.planeTOP = answer{2};
 % save image
save(fullfile(PATHNAME,FILENAME),'obj');
close all;