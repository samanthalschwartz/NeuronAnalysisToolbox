% nd opener for SIM data


% hannah is a 1x4 cell array where hannah{1,1} has the useful info

% z is z dimension
% c is color dimension
% plane is counter for entire stack
% image comes in as size(planes,1) = max(z) * max(c);
% runs as z1,c1
%           z1,c2
%           z1,c3
%           z2,c1
%           z2,c2
%           z2,c3


%% example to start. just get the input for z size and c size from user.
% this is a cell array hannahimage(:,1) is the image stack
%                                                 hannahimage(:,2) has the information for how to put the stack together





function im_out = loadnd2file(filepath,zstacksize,numcolors)
hannah = bfopen(filepath);
hannahimage = hannah{1,1};
xsize = size(hannahimage{1,1},1);
ysize = size(hannahimage{1,1},2);

im_out = zeros(xsize,ysize,zstacksize,numcolors);





end