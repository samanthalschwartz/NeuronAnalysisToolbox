
uiopen('F:\Hiester et al., 2017 Frontiers Submission\Matt Becker Data (For Review)\SEPGlua1_mch\20150609_SEPGlua1_mch_4_2.tif',1);
% uiopen('E:\Matt Becker Data (For Review)\SEPGlua1_mch\20150609_SEPGlua1_mch_4_2.tif',1);
sep = image(:,:,:,3);
cellfill = image(:,:,:,1);
ga = GeneralAnalysis;
dggsep = ga.imgDggCutoff(gaussf(sep));
image_out = varif(dggsep,7,'elliptic')
threshep = ga.imgThreshold_fixedUserInput(image_out)

spines = bdilation(threshep,2);
wshed = ga.watershed_timeseries(-gaussf(image_out).^1.5,1);
spines(wshed==1) = 0;
h = ga.viewMaskOverlayPerimStatic(dip_image(sep),spines);
dipmapping(h,'colormap',jet)
dipmapping(h,[0 5000]);
diptruesize(h,100);
% 
lb = label(spines,1);
newlb = zeros(size(lb,1),size(lb,2));
wb = waitbar(0);
for ll = 1:max(lb)
s = sum(spines,lb == ll,3);
newlb(s>0) = ll;
waitbar(ll/max(lb),wb);
end
close(wb)

%% -- 
% max project threshold
% watershed on max projection of image
% subtract watershed from max proj of threshold 
% filter labeling on size
% good to go?

%%
% sm = sum(spines,[],3)
% slb = label(sm>0)
msr = measure(lb,sep,{'sum','size'})
spineFo_array = zeros(max(lb),size(sep,3));
sepch = dip_image(sep);
for ll = 1:max(lb)
    spineFo_array(ll,:) = sum(sepch,lb == ll,[1 2]);   
end
spineFo_array = zeros(max(lb),size(sep,3));
sepch = dip_image(sep);
for ll = 1:max(lb)
    currmask = lb == ll;
    forsum = sepch.*currmask;
    spineFo_array(ll,:) = sum(forsum,3);
end
firstids = zeros(1,size(spineFo_array,1));
baselines = zeros(1,size(spineFo_array,1));
for aa = 1:size(spineFo_array,1)
firstids(aa) = find(spineFo_array(aa,:),1,'first'); 
if firstids(aa)< (size(spineFo_array,2)-2)
baselines(aa) = mean(spineFo_array(aa,firstids(aa):firstids(aa)+2));
else
   baselines(aa) = mean(spineFo_array(aa,firstids(aa)));
end
end
forplot = spineFo_array./baselines';
% forplot(forplot(:)>10) = 10;
figure;
hm = heatmap(forplot)
hm.GridVisible = 'off';
hm.Colormap = jet(50);
hm.ColorLimits = [0 3]
