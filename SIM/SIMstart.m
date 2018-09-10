% filename = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\071117\Gephyrin_Abeta_Bassoon\Geph488_Abeta561_Bassoon647_003_Reconstructed.nd2';
filepath = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\071117\Gephyrin_Abeta_Bassoon';
file = uipickfiles('Prompt','Pick Files','FilterSpec',filepath);
image = ndFileloader(file{1});
ch1 = image(:,:,:,1);
ch2 = image(:,:,:,2);
ch3 = image(:,:,:,3);

currch = ch1;

gch = gaussf(currch);
out = GeneralAnalysis.imgLaplaceCutoff(gch,[2 2 1],[2 2 1]);
thr = multithresh(single(out),2);
mask = out>thr(1);
a = overlay(currch,mask);
a{1} = mask.*max(currch);
a{2} = dip_image(currch).*~mask;
a{3} = dip_image(currch).*~mask;
a
GeneralAnalysis.viewMaskOverlayPerimStatic(dip_image(currch),mask)
overlay(dip_image(currch),mask)
labelim = label(mask,Inf,2,0);
msr = measure(labelim,currch,{'Size','Sum','Gravity'});
% % - optional for displaying image;
f = dipshow(currch); 
for ff = 1:size(msr,1)
    x = msr(ff).Gravity(1);
    y = msr(ff).Gravity(2);
    rectangle('Position',[x,y,1,1],'EdgeColor','r')
end
%






%%

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


