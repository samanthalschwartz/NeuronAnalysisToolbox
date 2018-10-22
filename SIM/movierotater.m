clear F G H
axis vis3d

firstview = [-90 90]; view(firstview);[az,el] = view;
lastview = [-130 0.8];
loops1 = 50;
steps = (lastview-firstview)/loops1
F(loops1) = struct('cdata',[],'colormap',[]);
for j = 1:loops1
    view(az+((j-1).*steps(1)),el+((j-1).*steps(2)));
    drawnow
    F(j) = getframe(gcf);
end

firstview = [-130 0.8]; view(firstview);[az,el] = view;
lastview = [-325 0];
loops1 = 150;
steps = (lastview-firstview)/loops1
G(loops1) = struct('cdata',[],'colormap',[]);
for j = 1:loops1
    view(az+(j.*steps(1)),el+(j.*steps(2)));
    drawnow
    G(j) = getframe(gcf);
end

firstview = [-325 0]; view(firstview);[az,el] = view;
lastview = [-450 90];
loops1 = 80;
steps = (lastview-firstview)/loops1
H(loops1) = struct('cdata',[],'colormap',[]);
for j = 1:loops1
    view(az+(j.*steps(1)),el+(j.*steps(2)));
    drawnow
    H(j) = getframe(gcf);
end
% 
% firstview = [-150 10]; view(firstview);[az,el] = view;
% lastview = [-90 90];
% loops1 = 50;
% steps = (lastview-firstview)/loops1
% I(loops1) = struct('cdata',[],'colormap',[]);
% for j = 1:loops1
%     view(az+(j.*steps(1)),el+(j.*steps(2)));
%     drawnow
%     I(j) = getframe(gcf);
% end
%%
filesavename = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\SIM_files\3d images\psd95bassoonabeta_roi3'
test= cat(2,F,G,H);
myVideo = VideoWriter(filesavename,'MPEG-4');
myVideo.FrameRate = 30;  % Default 30
myVideo.Quality = 100;
myVideo.ColorChannels
open(myVideo);
writeVideo(myVideo, test);
close(myVideo);




