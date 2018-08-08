%%-- script for automatically finding calcium transients and determining the extent of calcium
% this will prompt you to select a file to analyze and an excel sheet name
% to save it as
%  Make sure the excel file is closed to enable writing to it. Otherwise,
%  there will be an error.
%% run this cell to test out thresholding parameter using 'Run Section' button.
% then once you're happy with the parameter, 'Run' the entire script
threshold_parameter = 1.3; % smaller number is more sensitive, larger number is less sensitive

% Select the file you want to analyze: for now this is a 1 color time
% series tif file
excelfileinfo_maxextent = 'E:\Bill\Geph_IB-Cry2Olig-GFP plus blue light\calcium imaging\0 Mg QCTs and dendritic spikes\MaxExtensions.xlsx';
excelfileinfo_Fvalues = 'E:\Bill\Geph_IB-Cry2Olig-GFP plus blue light\calcium imaging\0 Mg QCTs and dendritic spikes\ROIintensities.xlsx';

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
[h,overlayim] = GeneralAnalysis.viewMaskOverlayPerimStatic(ca,bdilation(lbl>0,2));
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
% tt = table(ord_trace(1:10,1:10))
%% now plot
% h = msgbox('Plotting HeatMap....');
% times = 1:size(ord_trace,2);
% ord_trace_norm = ord_trace./median(ord_trace,2);
% g = figure;
% % wb = waitbar(0,'Plotting some things...');
% % cnt = 0;
% for ii = 1:size(ord_trace_norm,1)
%     for jj = 1:size(ord_trace_norm,2)
%         p = patch([times(jj),times(jj)+1,times(jj)+1,times(jj)],[ii-1, ii-1, ii, ii],ord_trace(ii,jj));
%         set(p,'FaceColor','flat','EdgeColor','none');
%         
%     end
% end
% title('rGecko Intensity (AU)','FontSize',16)
% xlabel('Frame','FontSize',16);
% yticks(1:size(ord_trace,1));
% yticklabels(num2str(round(spanval(ordx),2)'))
% ylabel('Span Value (pixels)','FontSize',16)
% xlim([0 size(ord_trace,2)]);
% ylim([0 size(ord_trace,1)]);
% c = colorbar;
% close(h);
%% write to excel file
% save max extents
xlswrite(excelfileinfo_maxextent,maxsizes',sheetname);
% save F values
xlswrite(excelfileinfo_Fvalues,ord_trace,sheetname);
%% save some of the file info
save(fullfile(pathname,[filename(1:end-4) '_analysisinfo']),'threshold_parameter','ord_trace','trace_raw','ordx','sumlbl','maxsizes','lbl','tt');
% saveas(g,fullfile(pathname,[filename(1:end-4) '_rawIntensityHeatMap']),'png');

