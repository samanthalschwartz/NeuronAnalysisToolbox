% Grouping Binned Information
% Make Page for each protein
% - FileName
% - Bins
% - NumAbeta Raw
% - NumObj
% - NormNumAbeta

close all; clear all;
load('G:\Hannah Dropbox SIM data\SIM_Files\071117\PSD95488_Abeta561_Bassoon647_001_Reconstructed_SIM.mat');
xs = obj.ch1.results.bins;
ys = obj.ch1.results.numabeta/obj.ch1.results.numobj;
