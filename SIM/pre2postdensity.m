
%%
ids = [2 1 3];
channelorderingstr = {'chABeta','PSD95','Bassoon'}; % channel abeta, channel 1, channel2
filepath = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\071917\PSD95488_Abeta561_Bassoon647_001_Reconstructed.nd2';
% ids = [2 1 3];
% channelorderingstr = {'chABeta','PSD95','Bassoon'}; % channel abeta, channel 1, channel2
% filepath = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\071117\PSD95_Abeta_Bassoon\PSD95488_Abeta561_Bassoon647_003_Reconstructed.nd2';

sor  = SIM();
sor.channelordering = ids;
sor.channelorderingstr = channelorderingstr;
sor.loadNDfile(filepath);

%-- measure ab com and plot with overlay
sor.make_cellmask;
sor.make_maskchAB;
sor.measure_AB;
% GeneralAnalysis.viewMaskOverlay(sor.abeta.image,sor.abeta.mask)

alllabels = sor.abeta.msr.ID;
abim = dip_image(0*sor.abeta.image);
for nn = alllabels
    xval = round(sor.abeta.msr.Gravity(1,nn));
    yval = round(sor.abeta.msr.Gravity(2,nn));
    zval = round(sor.abeta.msr.Gravity(3,nn));
    zmin = max(0,zval-1);
    zmax = min(size(sor.abeta.mask,3)-1,zval+1);
    abim(xval-1:xval+1,yval-1:yval+1,zmin:zmax) = 1;
end
joinchannels('rgb',stretch(abim),stretch(sor.ch1.image),stretch(sor.ch2.image))

%%
[B,C] = dipcrop(gcf);
% C =
%    468   416
%     11    10
% from movie C = 
% 827   140
%     21    20
s  = SIM();
s.channelordering = ids;
s.channelorderingstr = channelorderingstr;
s.loadNDfile(filepath);
s.ch1.image = dip_image(s.ch1.image);
s.ch1.image = s.ch1.image(C(1,1):(C(1,1)+C(2,1)),C(1,2):(C(1,2)+C(2,2)),:);
s.ch2.image = dip_image(s.ch2.image);
s.ch2.image = s.ch2.image(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2),:);
s.abeta.image = dip_image(s.abeta.image);
s.abeta.image = s.abeta.image(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2),:);
% joinchannels('rgb',stretch(s.abeta.image),stretch(s.ch1.image),stretch(s.ch2.image))
figure; 
camlight; lighting phong
scaleval = 2.8;

ch1 = dip_image(s.ch1.image);
szch1 = size(ch1);
colormatrix = ones(szch1(2),szch1(1),szch1(3));
% value = (max(ch1)-min(ch1))/1.5;
value = (max(ch1)+min(ch1))/scaleval
patch(isosurface(0:szch1(1)-1,0:szch1(2)-1,0:szch1(3)-1,double(ch1),value,colormatrix.*double(ch1)),'parent',gca,...
'facelighting','phong','facecolor','g','edgecolor','none');
hold on;

ch2 = dip_image(s.ch2.image);
szch2 = size(ch2);
colormatrix = ones(szch2(2),szch2(1),szch2(3));
value = (max(ch2)+min(ch2))/scaleval*1.4
patch(isosurface(0:szch2(1)-1,0:szch2(2)-1,0:szch2(3)-1,double(ch2),value,colormatrix.*double(ch2)),'parent',gca,...
'facelighting','phong','facecolor','b','edgecolor','none');
hold on;

ab = dip_image(s.abeta.image);
szab = size(ab);
colormatrix = ones(szab(2),szab(1),szab(3));
value = (max(ab)+min(ab))/scaleval;
patch(isosurface(0:szab(1)-1,0:szab(2)-1,0:szab(3)-1,double(ab),value,colormatrix.*double(ab)),'parent',gca,...
'facelighting','phong','facecolor','r','edgecolor','none');
hold on;
pbaspect([1 1 s.XYpxsize/(s.Zpxsize*.7)])

%%
s.make_cellmask;
s.make_maskch1;
s.make_maskch2;
s.measure_ch1;
s.measure_ch2;
p1 = round(s.ch1.thisCh.msr.Gravity(:,1));
p2 = round(s.ch2.thisCh.msr.Gravity(:,1));
pts = s.BresenhamPoints(p1',p2');
img = dip_image(s.abeta.image);
imcurve = zeros(1,size(pts,1));
thickness = 5;rang = floor(thickness/2);
for ii = 1:size(pts,1)
    imx = pts(ii,1);
    imy = pts(ii,2);
    imz = pts(ii,3);
    xmin = max(0,imx-rang);
    xmax = min(size(img,1)-1,imx+rang);
    ymin = max(0,imy-rang);
    ymax = min(size(img,2)-1,imy+rang);
    zmin = max(0,imz-rang);
    zmax = min(size(img,3)-1,imz+rang);
    imcurve(ii) = sum(img(xmin:xmax,ymin:ymax,zmin:zmax));
end  
figure;plot(imcurve);
% iso plotting 
s.make_maskchAB;
s.measure_AB;
