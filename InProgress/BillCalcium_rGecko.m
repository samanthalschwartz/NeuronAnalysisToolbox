%% for automatically finding calcium transients and determining the extent of calcium

% uiopen('E:\GephIntrabody Project\Bill rGecko\Composite.tif',1);
% uiopen('Z:\Bill\2018\8.6.18_ibOLIG_QCTs_well1\baseline_movie.tif',1);
uiopen('F:\GephIntrabody Project\Bill rGecko\180806\cell1\postZ_movie.tif',1);

ca = dip_image(image(:,:,1:end-1,1));

% uiopen('C:\Users\schwsama\Documents\Data\Hannah Calcium\Abeta1Post.tif',1)
% ca = image;

% first smooth image
ca_g = gaussf(ca,[1 1 0]);
% now normalize to median
ca_gnorm = ca_g;%- medif(ca_g);
% now take difference
ca_dzz = dzz(ca_gnorm,[2 2 3]);
% get absolute values
ab_ca_dzz = abs(ca_dzz);
gab_ca_dzz = gaussf(ab_ca_dzz,[1 1 0]);
% tt = threshold(ab_ca_dzz,'otsu');

tt = threshold(ab_ca_dzz^1.3,'otsu');
% lbl = label(tt,1,200,10^10);
lbl = label(tt,1,20,10^10);
tic;
msr = measure(lbl,ca,({'size','DimensionsCube'}));
toc;
GeneralAnalysis.viewMaskOverlayPerimStatic(ca,bdilation(lbl>0,2));


% tic
% circularities = msr.P2A;
% toc

%% plot CDF of max extension
maxsizes  = max(msr.DimensionsCube([1 2],:));
[F,X] = ecdf(maxsizes);
figure; plot(X,F); set(gca,'FontSize',16); 
xlabel('Max Extension (pixels) of Calcium Propagation','FontSize',16);
ylabel('Cumulative Distribution','FontSize',16);

%% plot heatmap of intensity within ROI over time

%label sum projection
slbl = sum(lbl,[],3);
sumlbl = label(slbl>0);
trace_raw = zeros(max(lbl),size(ca,3));
wb = waitbar(0,'Quantifying Calcium Change in ROIs...');
for ll = 1:max(sumlbl)%[1:8,10:max(lbl)]
%     tic
currmask = sumlbl==ll;
% bcurrmask = bdilation(currmask,1);
mask2use = repmat(currmask,1,1,size(ca,3));
sumval = sum(ca,mask2use,[1 2]);
% imgINmask = image*mask2use;
% sumval = sum(imgINmask,[],[1 2]);
% sizeval = sum(bcurrmask,[],[1,2]);
% tracemed = median(sumval);
% trace(ll,:) = sumval./tracemed;
trace_raw(ll,:) = sumval;
% toc
waitbar(ll/max(lbl),wb);
end
close(wb)


% plot in order of roi max extension
tic;
msrAll = measure(sumlbl,ca(:,:,1),({'P2A','size','DimensionsCube'}));
toc;

spanval = max(msrAll.DimensionsCube([1 2],:));
[~, ordx] = sort(spanval, 'ascend');
ord_trace = trace_raw(ordx,:);
% now plot
h = msgbox('Plotting HeatMap....');
times = 1:size(ord_trace,2);
ord_trace_norm = ord_trace./median(ord_trace,2);
figure;
% wb = waitbar(0,'Plotting some things...');
% cnt = 0;
for ii = 1:size(ord_trace_norm,1)
    for jj = 1:size(ord_trace_norm,2)
        p = patch([times(jj),times(jj)+1,times(jj)+1,times(jj)],[ii-1, ii-1, ii, ii],ord_trace_norm(ii,jj));
        set(p,'FaceColor','flat','EdgeColor','none');
        
    end
end
title('rGecko Intensity (AU)','FontSize',16)
xlabel('Frame','FontSize',16);
yticklabels(num2str(round(spanval(ordx),2)'))
ylabel('Span Value (pixels)','FontSize',16)
xlim([0 size(ord_trace,2)]);
ylim([0 size(ord_trace,1)]);
c = colorbar;
close(h);


%%
h = dipshow(sum(ca,[],3));
[roi, v] = diproi(h)
