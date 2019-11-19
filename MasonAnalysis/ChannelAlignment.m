datadir = 'Z:\Mason\082819_dualslice';
%--- find filebasenames
ndfiles = dir(fullfile(datadir,'*.nd'));
filebasearr = cell(1,numel(ndfiles));
for nn = 1:numel(ndfiles)
   filebasearr{nn} = ndfiles(nn).name(1:end-3);  
end
for ff = 1:numel(filebasearr)
filebase = filebasearr{ff};
chstr_red = '_w1561';
chstr_green = '_w2488';
chstr_blue = '_w3640';
disp('loading files...');
loadtiff(fullfile(datadir,[filebase,chstr_red,'.TIF']));
ch_red_z = ans;
ch_red = max(ch_red_z,[],3);

loadtiff(fullfile(datadir,[filebase,chstr_green,'.TIF']));
ch_green_z = ans;
ch_green = max(ch_green_z,[],3);

loadtiff(fullfile(datadir,[filebase,chstr_blue,'.TIF']));
ch_blue_z = ans;
ch_blue = max(ch_blue_z,[],3);
disp('shifting files...');

% find shift from green to red
 shift_green2red = findshift(ch_red,ch_green,'iter');
 new_ch_green = shift(ch_green,shift_green2red);

% find shift from blue to red
shift_blue2red = findshift(ch_red,ch_blue,'iter');
new_ch_blue = shift(ch_blue,shift_blue2red);
 
new_green = zeros(size(ch_red_z));
new_blue = zeros(size(ch_red_z));
for zz = 1:size(ch_red_z,3)
    new_green(:,:,zz) = shift(ch_green_z(:,:,zz),shift_green2red);
     new_blue(:,:,zz) = shift(ch_blue_z(:,:,zz),shift_blue2red);
end
% joinchannels('rgb',max(ch_red_z,[],3),max(new_green,[],3),max(new_blue,[],3));
% joinchannels('rgb',ch_red_z*2,new_blue,ch_blue_z)
disp('saving files...');

LibTiff(new_green,fullfile(datadir,[filebase,chstr_green '_shifted']));
LibTiff(new_blue,fullfile(datadir,[filebase,chstr_blue '_shifted']));
disp('all done!');
end
disp('done,done');