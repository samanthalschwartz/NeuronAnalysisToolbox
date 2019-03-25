%%-- script for automatically finding calcium transients and determining the extent of calcium
% this will prompt you to select a file to analyze and an excel sheet name
% to save it as
%  Make sure the excel file is closed to enable writing to it. Otherwise,
%  there will be an error.
%% run this cell to test out thresholding parameter using 'Run Section' button.
% then once you're happy with the parameter, 'Run' the entire script
threshold_parameter = 1.1; % smaller number is more sensitive, larger number is less sensitive

% Select the file you want to analyze: for now this is a 1 color time
% series tif file
excelfileinfo_maxextent = 'E:\Bill\Geph_IB-Cry2Olig-GFP plus blue light\calcium imaging\0 Mg QCTs and dendritic spikes\MaxExtensions_REDO.xlsx';
excelfileinfo_Fvalues = 'E:\Bill\Geph_IB-Cry2Olig-GFP plus blue light\calcium imaging\0 Mg QCTs and dendritic spikes\ROIintensities_REDO.xlsx';
excelfileinfo_Fnorm = 'E:\Bill\Geph_IB-Cry2Olig-GFP plus blue light\calcium imaging\0 Mg QCTs and dendritic spikes\MaxIntensities_REDO.xlsx';


startdir = 'E:\Bill\Geph_IB-Cry2Olig-GFP plus blue light\calcium imaging';
[filename, pathname, filterindex] = uigetfile(fullfile(startdir,'*.tif'));
prompt = {'Select a sheet name for this cell/dataset'};
title = 'Input';
dims = [1 35];
definput = {filename(1:end-4)};
sheetname = inputdlg(prompt,title,dims,definput);
ca = GeneralAnalysis.loadtiff_1ch(fullfile(pathname,filename));
% if last frame is empty uncomment this line
ca = ca(:,:,0:end-1);
% first smooth image
ca_g = gaussf(ca,[1 1 0]);
% now normalize to median
ca_gnorm = ca_g;%- medif(ca_g);
% now take difference
ca_dzz = dzz(ca_gnorm,[2 2 3]);
% get absolute values
ab_ca_dzz = abs(ca_dzz);
gab_ca_dzz = gaussf(ab_ca_dzz,[1 1 0]);
tt = threshold(ab_ca_dzz^threshold_parameter,'otsu');
lbl = label(tt,1,20,10^10);
% this is to view the identified masks over the image
[h,overlayim] = GeneralAnalysis.viewMaskOverlayPerimStatic(ca,bdilation(lbl>0,2),bone,[0 0 0]);
%%
close(h);
tic;
msr = measure(lbl,ca,({'size','DimensionsCube'}));
toc;
maxsizes  = max(msr.DimensionsCube([1 2],:));
% now get this setup to write to excel file
if isempty(sheetname)
sheetname = filename(1:end-4);
else 
    sheetname = sheetname{1};
end
%% plot heatmap of intensity within ROI over time

%label sum projection
slbl = sum(lbl,[],3);
sumlbl = label(slbl>0);
trace_raw = zeros(max(sumlbl),size(ca,3));
wb = waitbar(0,'Quantifying Calcium Change in ROIs...');
for ll = 1:max(sumlbl)
%     tic
currmask = sumlbl==ll;
mask2use = repmat(currmask,1,1,size(ca,3));
sumval = sum(ca,mask2use,[1 2]);
trace_raw(ll,:) = sumval;
% toc
waitbar(ll/max(sumlbl),wb);
end
close(wb)

% plot in order of roi max extension
tic;
msrAll = measure(sumlbl,ca(:,:,1),({'P2A','size','DimensionsCube'}));
toc;
spanval = max(msrAll.DimensionsCube([1 2],:));
[~, ordx] = sort(spanval, 'ascend');
ord_trace = trace_raw(ordx,:);
%%
Fval = zeros(max(lbl),1);
Fnorm = zeros(max(lbl),1);
% loop through each label
wb = waitbar(0,'waiting to calculate max intensities....');
for ll = 1:max(lbl)
    % now calculate sum intensity inside mask
    currmask = lbl==ll;
    sumvaltrace = sum(ca,currmask,[1 2]);
    [Y,I] = max(single(sumvaltrace),[],3);
    maxmask = currmask(:,:,I-1);
    imagetrace = sum(ca,repmat(maxmask,1,1,size(ca,3)),[1 2]);
    baseline = median(imagetrace);
    Fval(ll) = Y;
    Fnorm(ll) = Y/baseline;
    waitbar(ll/max(lbl),wb);
end
close(wb);
% sheetname = [foldernames{ff} '-' savenames{pf}];
goodids = Fnorm>=1.05;
goodFnorms = Fnorm.*goodids;
goodmaxs = maxsizes'.*goodids;
goodmaxs(goodmaxs == 0) = NaN;

goodFnorms(goodFnorms == 0) = NaN;
xlswrite(excelfileinfo_Fnorm,[Fval,Fnorm,goodFnorms],sheetname);
%% write to excel file
% save max extents
newmaxsizes = [maxsizes',goodmaxs];
xlswrite(excelfileinfo_maxextent,newmaxsizes,sheetname);
% save F values
xlswrite(excelfileinfo_Fvalues,ord_trace,sheetname);
%% save some of the file info
save(fullfile(pathname,[filename(1:end-4) '_analysisinfo']),'threshold_parameter','ord_trace','trace_raw','ordx','sumlbl','maxsizes','lbl','tt');
% saveas(g,fullfile(pathname,[filename(1:end-4) '_rawIntensityHeatMap']),'png');

