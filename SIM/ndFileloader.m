function image = ndFileloader(filename)
if nargin<1
    [name, path] = uigetfile('*.nd2');
    filename = fullfile(path,name);
end
hannah = bfopen(filename);
% use first input to determine #channels, #z-planes
imsize = size(hannah{1,1}{1,1}); %size of single image frame
infostr = hannah{1,1}{1,2};
info = parseNDtext(infostr);
image = zeros([imsize,info.totalZ,info.totalCol]);
for ii = 1:size(hannah{1,1},1)
    infostr = hannah{1,1}{ii,2};
    info = parseNDtext(infostr);
    image(:,:,info.currZ,info.currCol) = hannah{1,1}{ii,1};
end
end
function info = parseNDtext(infostr)

splstr1 = strsplit(infostr,{';'});

planebool = cellfun(@(x) contains(x,'plane '),splstr1);
planestr = splstr1{planebool};
planeparts = strsplit(planestr,{' plane ','/'});
info.currPlane = str2double(planeparts{2});
info.totalPlane = str2double(planeparts{3});

zbool = cellfun(@(x) contains(x,'Z='),splstr1);
zstr = splstr1{zbool};
zparts = strsplit(zstr,{'=','/'});
info.currZ = str2double(zparts{2});
info.totalZ = str2double(zparts{3});

colbool = cellfun(@(x) contains(x,'C='),splstr1);
colstr = splstr1{colbool};
colparts = strsplit(colstr,{'=','/'});
info.currCol = str2double(colparts{2});
info.totalCol = str2double(colparts{3});
end