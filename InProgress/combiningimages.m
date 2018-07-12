load('C:\Users\schwsama\Documents\Data\zapERtrap\cell3_AshleyFile.mat');
body = uint16(aa.cellFill.image(:,:,1));
load('C:\Users\schwsama\Documents\Data\zapERtrap\cell3dendrites_AshleyFile.mat');
dends = uint16(aa.cellFill.image(:,:,1));
cc = normxcorr2(dends,body);
[ypeak, xpeak] = find(cc==max(cc(:)));
dipshow(cc)

diffx = size(dends,2)-xpeak;
diffy = size(dends,1)-ypeak;

newim = zeros(diffy+size(body,2),diffx+size(body,1));

newim(1:size(dends,1),1:size(dends,2)) = dends;
newim(end-size(body,1)+1:end,end-size(body,2)+1:end) = body;
dipshow(newim)