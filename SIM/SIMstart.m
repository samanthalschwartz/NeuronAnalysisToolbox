% filename = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\071117\Gephyrin_Abeta_Bassoon\Geph488_Abeta561_Bassoon647_003_Reconstructed.nd2';
%% load the image
filepath = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\071117\Gephyrin_Abeta_Bassoon';
file = uipickfiles('Prompt','Pick Files','FilterSpec',filepath);
image = ndFileloader(file{1});
ch1 = image(:,:,:,1);
ch2 = image(:,:,:,2);
ch3 = image(:,:,:,3);
%% pick the channel
currch = ch1;
abetach = ch2;
%% mask the channel
gch = gaussf(currch);
out = GeneralAnalysis.imgLaplaceCutoff(gch,[2 2 1],[2 2 1]);
thr = multithresh(single(out),2);
mask = bdilation(out>thr(1),1);
GeneralAnalysis.viewMaskOverlay(dip_image(currch),mask)
%% find COM
labelim = label(mask,Inf,4,0);
msr = measure(labelim,currch,{'Size','Sum','Gravity'});
%%
% % - optional for displaying image;
f = dipshow(currch); 
for ff = 1:size(msr,1)
    x = msr(ff).Gravity(1);
    y = msr(ff).Gravity(2);
    rectangle('Position',[x,y,1,1],'EdgeColor','r')
end
%% calculate density as a function of distance from COM
%-- try to find the cell from the labels:
 sumim = ch1+ch2+ch3;
 maxsumim = max(sumim,[],3);
 gmaxsumim = gaussf(maxsumim,10);
%  threshs = multithresh(single(gmaxsumim),2);
%  cellmask = gmaxsumim>threshs(1);
 openeninimage_out = opening(gmaxsumim,20,'rectangular'); 
 threshs = multithresh(single(openeninimage_out),2);
 cellmask = openeninimage_out>threshs(1);
 %%
 cellmask3 = repmat(cellmask,1,1,size(sumim,3));
 h = GeneralAnalysis.viewMaskOverlay(dip_image(sumim),cellmask3)
%% 
% -- calculate 3D distance transform inside cell
%----something isn't right here!@-----
distmask = dt(~mask);
distcellmask = distmask.*cellmask;

nbins = 20;
inc = max(distcellmask)/nbins;

abetach = ch2;
bins = 0:10;
density = zeros(size(bins));
for nn = bins
    if nn == 1
        curr_distmask = (distcellmask<=nn).*cellmask3;
    else
        curr_distmask = ((distcellmask>(nn-1)) & distcellmask<=nn).*cellmask3;
    end
    cnt = nn+1;
   density(nn+1) = sum(abetach.*curr_distmask)./sum(curr_distmask(:));
end
avgdensity = sum(abetach.*cellmask3)./sum(cellmask3(:));
figure; plot(bins,density./avgdensity);
ylim([0.90 1.2])

xlabel('Distance from Mask (pixels)');
ylabel('Relative Intensity Density')

%%
distmask = dt(~mask);
distcellmask = distmask.*cellmask;

nbins = 20;
inc = max(distcellmask)/nbins;
density = zeros(nbins,1);
for nn = 1:nbins
    if nn == 1
        curr_distmask = (distcellmask<=(inc*nn)).*cellmask;
    else
        curr_distmask = (distcellmask>(inc*(nn-1)) & distcellmask<=(inc*nn)).*cellmask;
    end
   density(nn) = sum(abetach.*curr_distmask)./sum(curr_distmask);
end
avgdensity = sum(abetach.*cellmask)./sum(cellmask);
xs = inc:inc:max(distcellmask);
figure; plot(xs,density/avgdensity)
%%
distmask = dt(~mask);
distcellmask = distmask.*cellmask;

nbins = 20;
inc = max(distcellmask)/nbins;
density = zeros(nbins,1);
for nn = 1:10
    if nn == 1
        curr_distmask = (distcellmask<=nn).*cellmask;
    else
        curr_distmask = (distcellmask>nn) & distcellmask<=nn.*cellmask;
    end
   density(nn) = sum(abetach.*curr_distmask)./sum(curr_distmask);
end
avgdensity = sum(abetach.*cellmask)./sum(cellmask);
xs = inc:inc:max(distcellmask);
xs = 1:10;
figure; plot(xs,density/avgdensity)


%% old stuff

filename = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\071117\Gephyrin_Abeta_Bassoon\Geph488_Abeta561_Bassoon647_003_Reconstructed.nd2';
% filename = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\071117\GluA1_Abeta_bassoon\GluA1488_Abeta561_Bassoon647_002_Reconstructed.nd2';

image = ndFileloader(filename);
% image = ndFileloader();
geph = image(:,:,:,1);
abeta = image(:,:,:,2);
bassoon = image(:,:,:,3);

gbassoon = gaussf(bassoon);
out = GeneralAnalysis.imgLaplaceCutoff(gbassoon,[1 1 1],[1 1 1]);
bassoonmask = threshold(out);
label_bassoonmask = label(bassoonmask,Inf,2,0);
msr = measure(label_bassoonmask,bassoon,{'SurfaceArea','Size','Sum','Gravity'});
msr = measure(label_bassoonmask,bassoon,{'Gravity'});
f = dipshow(bassoon); 
for ff = 1:size(msr,1)
    x = msr(ff).Gravity(1);
    y = msr(ff).Gravity(2);
    rectangle('Position',[x,y,1,1],'EdgeColor','r')
end
%%
ab = gaussf(ch2);
out = GeneralAnalysis.imgLaplaceCutoff(ab,[2 2 1],[2 2 1]);
thr = multithresh(single(out),2);
mask = out>thr(1);
GeneralAnalysis.viewMaskOverlay(dip_image(ch2),mask)
sch2 = sum(dip_image(ch2),[],3);
sch2_series = repmat(sch2,1,1,size(ch2,3));
wshed = GeneralAnalysis.watershed_timeseries(-gaussf(sch2_series,1),1);
lbl = label(mask.*~wshed,1,8)
msr = measure(lbl,dip_image(ch2),{'Size','Sum'});
figure; histogram(msr.Sum./msr.Size)
%%
omeMeta = test{1, 4};
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
voxelSizeXdefaultValue = omeMeta.getPixelsPhysicalSizeX(0).value();           % returns value in default unit
voxelSizeXdefaultUnit = omeMeta.getPixelsPhysicalSizeX(0).unit().getSymbol(); 
voxelSizeZdefaultValue = omeMeta.getPixelsPhysicalSizeZ(0).value();           % returns value in default unit
voxelSizeZdefaultUnit = omeMeta.getPixelsPhysicalSizeZ(0).unit().getSymbol(); 