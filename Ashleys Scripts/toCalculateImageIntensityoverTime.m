clear all
close all
load('E:\Zapalog project\GluA1_norelease_control\slip1\2merge_AshleyFile.mat'); % load the file
% call dipimage sum 
% this takes in an image, a mask to use, the dimensions of summation
% squeeze just reduces the dimensionality
% repmat is being used to make cellfillmask the same dimensions as the
% surfaceCargo image
cellfillmask = repmat(max(aa.cellFill.mask,[],3),1,1,size(aa.surfaceCargo.image,3));
surfaceCimage = aa.surfaceCargo.image;
surfaceCmask = aa.surfaceCargo.mask;
rawtrace = squeeze(sum(surfaceCimage,cellfillmask,[1 2]));
%rawtrace_s is a new trace that has different dimensions then rawtrace (dipimage) so
%that it can be plotted in a nice figure
rawtrace_s = single(rawtrace);
figure; plot(rawtrace_s)
%%
% to correct for background
% bdilation is making the cellfillmask larger
mask4bg  = bdilation(cellfillmask|surfaceCmask,10);
bgim_out = GeneralAnalysis.regionfill_timeseries(surfaceCimage,mask4bg);
%%
% now get same trace as above but for bgim_out and subtract this from hi
% trace to get background subtracted/normalized image intensity trace
% also set all negative values to equal 0
corr_image=surfaceCimage-bgim_out;
corr_image(corr_image < 0) = 0;

trace_bgcorr = squeeze(single(sum(corr_image,aa.cellFill.mask,[1 2])));
figure; plot(trace_bgcorr)
% to normalize the background trace, run the line of code below instead of
% the line above
% figure; plot(trace_bgcorr./trace_bgcorr(1))

% xlswrite to export into Excel
savedir = 'C:\Users\Ashley\Documents\MATLAB\NeuronAnalysisToolbox\Ashleys Scripts';
savename = 'NoReleaseControls_GluA1_ExportedfromMATLAB.xls';
xlswrite(fullfile(savedir,savename),trace_bgcorr)