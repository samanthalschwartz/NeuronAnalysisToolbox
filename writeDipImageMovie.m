function writeDipImageMovie(h,filename,options)
%h: figure handle to window 
    % - this just allows you to set up the preferred mapping beforehand etc.
% filename: full filename and path to use for saving
% options: options.framerate = 20; %default
%          options.slices = allslices %default. also could do ex. [0:20]
% vid = VideoWriter(filename,'MPEG-4');
vid = VideoWriter(filename,'MPEG-4');
seq = dipgetimage(h);
nframes = size(seq{1},3);
if nargin<3
    options.framerate = 20;
    options.slices = 0:nframes-1;
else 
    if ~isfield(options,'framerate')
        options.framerate = 20;
    end
    if ~isfield(options,'slices')
        options.slices = 0:nframes-1;
    end
end
vid.FrameRate = options.framerate;
open(vid);
for ii=options.slices
    dipmapping(h,'slice',ii);
    frame = getframe(h);
    writeVideo(vid,frame);
end
close(vid);
end