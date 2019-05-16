vd = viewer3D;
filename = 'Y:\Theresa\Images\0325slice_AAVSphHaloTent_gephIBgfp_tdTom_n3_pre_3drecon.tif';
% vd.filepath = filename;
vd.loadimagefile(filename);
vd.initialize3Dimage();
vd.makepanel();
vd.dimension = [ 1 1 0.25];

syn = vd.image{3};
scaleval = 2.5;
ab = gaussf(syn,[1 1 0]);
value = (max(ab)+min(ab))/scaleval;
abmask = ab>value;
lb = label(abmask);
lb(lb == 2) = 0;
vd.image{3} = gaussf(lb,[1 1 0]);
vd.initialize3Dimage();

vd.patch{1}.FaceAlpha = 0.2
vd.endpatch{1}.FaceAlpha = 0.2

camlight; material dull;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'ZTickLabel',[]);
%%
axis vis3d
clear F G
firstview = [-130 -90]; view(firstview);[az,el] = view;
lastview = [-130 60];
loops1 = 50;
steps = (lastview-firstview)/loops1
F(loops1) = struct('cdata',[],'colormap',[]);
for j = 1:loops1
    view(az+((j-1).*steps(1)),el+((j-1).*steps(2)));
    drawnow
    F(j) = getframe(gcf);
end
G(loops1) = struct('cdata',[],'colormap',[]);
steps = (lastview-firstview)/loops1;
[az,el] = view;
for j = 1:loops1
    view(az-((j-1).*steps(1)),el-((j-1).*steps(2)));
    drawnow
    G(j) = getframe(gcf);
end
filesavename = 'Y:\Theresa\Images\0325slice_AAVSphHaloTent_gephIBgfp_tdTom_n3_pre_3drecon_movie_try1'
myVideo = VideoWriter(filesavename,'MPEG-4');
myVideo.FrameRate = 15;  % Default 30
myVideo.Quality = 100;
open(myVideo);
test= cat(2,F,G);
writeVideo(myVideo, test);
close(myVideo);

