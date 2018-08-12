load('E:\Sam\Data\MJK_zapERtrap_for_sam\021318_NL1_cells2_3\cell2_global_AshleyFile.mat')
filesavename = 'C:\Users\schwsama\Documents\Data\myfile';
%%
dipisosurface_samchange(aa.cleanedcargomask(:,:,7:end));
%% now make it the size we want and remove buttons


ax = gca;
azend = -50;
elend = 30;
view([azend elend])
light('parent',ax,'position',get(ax,'CameraPosition'));
light('parent',ax,'position',-get(ax,'CameraPosition')); % one light on each side of the surface

%%
set(ax,'XGrid','off','YGrid','off','ZGrid','off');
view([-45 -45])
% light('parent',ax,'position',get(ax,'CameraPosition'));
% light('parent',ax,'position',-get(ax,'CameraPosition')); % one light on each side of the surface

azend = -50;
elend = 30;
view([azend elend])
% view([0 -0])
axis tight manual
ax.NextPlot = 'replaceChildren';
set(gca,'cameraviewanglemode','manual')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'ZTickLabel',[]);



fig = gcf;
view([0 -90])
% light('parent',ax,'position',get(ax,'CameraPosition'));
% light('parent',ax,'position',-get(ax,'CameraPosition')); % one light on each side of the surface

cb = colorbar;
set(gca,'clim',[3 30])
set(get(cb,'title'),'string','Frame After Release','FontSize',16)

loops1 = 100;
F(loops1) = struct('cdata',[],'colormap',[]);
for j = 1:loops1
    view([-j+1 -90+j-1])
    drawnow
    F(j) = getframe(fig);
end
[az,el] = view;
loops2 = 50;
G(loops2) = struct('cdata',[],'colormap',[]);
azend = -50;
elend = 30;
azinc = (az-azend)/loops2;
elinc = (el-elend)/loops2;
for j = 1:loops2
   view([az - j*azinc el-j*elinc])
    drawnow
    G(j) = getframe(fig);
end

test= cat(2,F,G);
% fig2 = figure;
% set(fig2,'Position',fig.Position)
% movie(test)

myVideo = VideoWriter(filesavename,'MPEG-4');
myVideo.FrameRate = 25;  % Default 30
myVideo.Quality = 100;
open(myVideo);
writeVideo(myVideo, test);
close(myVideo);