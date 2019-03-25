testim = 0.*obj.ch1.distance_mask;
for zz = 1:size(obj.abeta.image,3)
    % 
goodmask = obj.ch1.distance_mask(:,:,zz)<distance;
% find how many abeta COM in this frame
currmask = obj.abeta.COM_image(:,:,zz-1).*goodmask;
numabeta = sum(currmask(:));
sz_x = size(goodmask,1);
sz_y = size(goodmask,2);
% simulation
simvals = zeros(numabeta,2);
cnt = 0;
while cnt < numabeta
    xval = round(rand(1,1)*(sz_x-1)) + 1;
    yval = round(rand(1,1)*(sz_y-1)) + 1;
    if goodmask(xval,yval)
        cnt=cnt+1;
        simvals(cnt,:) = [xval,yval];
        testim(xval,yval,zz) = 1; 
    end
end
end
joinchannels('rgb',bdilation(logical(obj.abeta.COM_image)),bdilation(logical(testim)))
joinchannels('rgb',obj.abeta.COM_image.*(obj.ch1.distance_mask<16),testim,obj.ch1.mask)