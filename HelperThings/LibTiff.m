function [] = LibTiff(Vol,inputname)
%LibTiff wrote by Michael @ mic.muenter@uni-luebeck.de
% Input:uint8 Volume (x,y,z)
% INPUT: filepath

% INPUT example: Vol = 255.*ones(300,100,200); LibTiff(Vol); 

Vol = uint32(Vol);

t = Tiff([inputname,'.tiff'],'w'); % Filename by variable name
tagstruct.ImageLength = size(Vol,1); % image height
tagstruct.ImageWidth = size(Vol,2); % image width
tagstruct.Photometric = Tiff.Photometric.MinIsBlack; % https://de.mathworks.com/help/matlab/ref/tiff.html
tagstruct.BitsPerSample = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';

tic
setTag(t,tagstruct)
write(t,squeeze(Vol(:,:,1)));
for i=2:size(Vol,3) % Write image data to the file
    writeDirectory(t);
    setTag(t,tagstruct)
    write(t,squeeze(Vol(:,:,i))); % Append
end
close(t);
toc

disp('Tiff File saved')
end