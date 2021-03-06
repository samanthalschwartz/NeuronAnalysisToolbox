
function dip_makemovie(varargin)
% inputs: filename,image, then remaining inputs are the same for dipshow
% call
filename = varargin{1};
image = varargin{2};
writeObj = VideoWriter(filename,'MPEG-4');
writeObj.FrameRate = 10;
open(writeObj);
if numel(varargin)==2 % no extra dipshow inputs
    for ii=1:size(image,3)
        h = dipshow(image(:,:,ii-1));
        F = getframe(h);
        writeVideo(writeObj,F);
        close(h);
    end
elseif numel(varargin)==4 % extra dipshow inputs
    for ii=1:size(image,3)
        h = dipshow(image(:,:,ii-1),varargin{3},varargin{4});
        F = getframe(h);
        writeVideo(writeObj,F);
    end
elseif numel(varargin<2)
    display('Error: Must have at least 2 inputs')
end


end