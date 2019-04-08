load('Z:\Sam\Hannah SIM_Files\071117\Geph488_Abeta561_Bassoon647_001_Reconstructed_SIM.mat')
joinchannels('rgb',obj.abeta.mask,obj.ch1.mask,obj.ch2.mask)
[B,C] = dipcrop(gcf);

ch1.image = obj.ch1.image(C(1,1):(C(1,1)+C(2,1)),C(1,2):(C(1,2)+C(2,2)),:);
ch2.image = obj.ch2.image(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2),:);
abeta.image = obj.abeta.image(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2),:);

ch1.mask = obj.ch1.mask(C(1,1):(C(1,1)+C(2,1)),C(1,2):(C(1,2)+C(2,2)),:);
ch2.mask = obj.ch2.mask(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2),:);
abeta.mask = obj.abeta.mask(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2),:);


msrch1 = measure(ch1.mask,single(ch1.mask),obj.measurements);
msrch2 = measure(ch2.mask,single(ch2.mask),obj.measurements);

testim1= 0.*obj.ch1.mask;
testim1(round(msr.Gravity(1)),round(msr.Gravity(2)),round(msr.Gravity(3))) = 1;
joinchannels('rgb',testim1,testim1*0,ch1.mask)

p1= round(msrch1.Gravity);
p2= round(msrch2.Gravity);
pts = obj.BresenhamPoints(p1',p2');
img = abeta.mask;
imcurve = zeros(1,size(pts,1));
thickness = 5;rang = floor(thickness/2);
testim1= 0.*ch1.mask;
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
    testim1(imx,imy,imz) = 1;
end  

% ch1 to ch2
u = p1'-p2';
ThetaInDegrees = atan2d(norm(cross(u,v)),dot(u,v));

