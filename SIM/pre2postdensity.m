ids = [2 1 3];
channelorderingstr = {'chABeta','Gephyrin','Bassoon'}; % channel abeta, channel 1, channel2
filepath = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\071917\Gephyrin_Bassoon\Gephyrin488_Abeta561_Bassoon647_004_Reconstructed.nd2';


s  = SIM();
s.channelordering = ids;
s.channelorderingstr = channelorderingstr;
s.loadNDfile(filepath);

s.make_cellmask;
s.make_maskch1;
s.measurements = [s.measurements,'Gravity','Center'];
s.measure_ch1;

p1 = round(s.ch1.thisCh.msr.Gravity(:,1));
p2 = round(s.ch1.thisCh.msr.Gravity(:,2));

pts = s.BresenhamPoints(p1',p2');

dispimg = dip_image(false(size(s.ch1.image)));
for ii = 1:size(pts,1)
    imx = pts(ii,1);
    imy = pts(ii,2);
    imz = pts(ii,3);
    dispimg(imx,imy,imz) = true;
end    
out = slice_op('bdilation',dispimg,2);

img = dip_image(s.ch1.image);
imcurve = zeros(1,size(pts,1));
for ii = 1:size(pts,1)
    imx = pts(ii,1);
    imy = pts(ii,2);
    imz = pts(ii,3);
    imcurve(ii) = img(imx,imy,imz);
end  

%%
ids = [2 1 3];
channelorderingstr = {'chABeta','Gephyrin','Bassoon'}; % channel abeta, channel 1, channel2
filepath = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\071917\Gephyrin_Bassoon\Gephyrin488_Abeta561_Bassoon647_004_Reconstructed.nd2';


s  = SIM();
s.channelordering = ids;
s.channelorderingstr = channelorderingstr;
s.loadNDfile(filepath);

% joinchannels('rgb',s.ch1.image,s.ch2.image)
[B,C] = dipcrop(gcf)
% C =
%    468   416
%     11    10
s.ch1.image = dip_image(s.ch1.image);
s.ch1.image = s.ch1.image(C(1,1):(C(1,1)+C(2,1)),C(1,2):(C(1,2)+C(2,2)),:);
s.ch2.image = dip_image(s.ch2.image);
s.ch2.image = s.ch2.image(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2),:);
s.abeta.image = dip_image(s.abeta.image);
s.abeta.image = s.abeta.image(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2),:);
joinchannels('rgb',stretch(s.abeta.image),stretch(s.ch1.image),stretch(s.ch2.image))


s.make_cellmask;
s.make_maskch1;
s.make_maskch2;
s.measurements = [s.measurements,'Gravity','Center'];
s.measure_ch1;
s.measure_ch2;
p1 = round(s.ch1.thisCh.msr.Gravity(:,1));
p2 = round(s.ch2.thisCh.msr.Gravity(:,1));
pts = s.BresenhamPoints(p1',p2');
img = dip_image(s.abeta.image);
imcurve = zeros(1,size(pts,1));
thickness = 5;
for ii = 1:size(pts,1)
    imx = pts(ii,1);
    imy = pts(ii,2);
    imz = pts(ii,3);
    rang = floor(thickness/2);
    imcurve(ii) = sum(img(imx-rang:imx+rang,imy-rang:imy+rang,repmat(imz,1,5)));
end  
figure;plot(imcurve);
