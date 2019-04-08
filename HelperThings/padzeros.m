function image_out = padzeros(image_in,borderpxs)
% add left x border
ybord = zeros(borderpxs,size(image_in,1));
newim1 = [ybord;image_in;ybord];
xbord = zeros(size(newim1,2),borderpxs);
image_out = [xbord,newim1,xbord];
end