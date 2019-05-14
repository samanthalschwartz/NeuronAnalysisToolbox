
        %%
        
class viewer3D
uiopen('Y:\Theresa\Images\0325slice_AAVSphHaloTent_gephIBgfp_tdTom_n3_pre_3drecon.tif',1)
cellfill = image(:,:,:,1);
gephyrin = image(:,:,:,2);
syn = image(:,:,:,3);
close all




camlight; lighting phong
whitebg('w')

scaleval = 5.5;
ch1 = gaussf(cellfill,[1 1 0]);

function p = makepatch(data,scaleval)
szd = size(data);
colormatrix = ones(szd(2),szd(1),szd(3));
value = (max(data)+min(data))/scaleval;
p = patch(isosurface(0:szd(1)-1,0:szd(2)-1,1:szd(3),double(data),value,colormatrix.*double(data)),'parent',gca...
,'facelighting','phong','facecolor','b','edgecolor','none','FaceAlpha',.5);
end

function h = plotchannels(ax,channels,scalevals)
p1 = makepatch(ax,channels.ch1,scalevals.ch1); hold on;
p2 = makepatch(ax,channels.ch2,scalevals.ch2); hold on;
p3 = makepatch(ax,channels.ch3,scalevals.ch3); hold on;
end

newfig = figure;
a = uicontrol('Parent',newfig,'Style','slider','Position',[121,100,410,23],...
    'value',zeta, 'min',0, 'max',1);
b = uicontrol('Parent',newfig,'Style','slider','Position',[121,50,410,23],...
    'value',zeta, 'min',0, 'max',1);
c = uicontrol('Parent',newfig,'Style','slider','Position',[121,150,410,23],...
    'value',zeta, 'min',0, 'max',1);

b.Callback = @slider_callback;
function slider_callback(hObject, callbackdata)
    A_t = num2str(hObject.Value);
    plotchannels(a.Value,b.
end



%%
hold on;
scaleval = 5;
ch2 = gaussf(gephyrin,[1 1 0]);
szch2 = size(ch2);
colormatrix = ones(szch2(2),szch2(1),szch2(3));
value = (max(ch2)+min(ch2))/scaleval*1.4;
p1 = patch(isosurface(0:szch2(1)-1,0:szch2(2)-1,0:szch2(3)-1,double(ch2),value,colormatrix.*double(ch2)),'parent',gca,...
'facelighting','phong','facecolor','g','edgecolor','none','FaceAlpha',.5);
hold on;


scaleval = 2.5;
ab = gaussf(geph,[1 1 0]);
value = (max(ab)+min(ab))/scaleval;
abmask = ab>value;
lb = label(abmask);
lb(lb == 2) = 0;
szab = size(lb);
colormatrix = ones(szab(2),szab(1),szab(3));
patch(isosurface(0:szab(1)-1,0:szab(2)-1,0:szab(3)-1,double(lb),1,colormatrix.*double(lb)),'parent',gca,...
'facelighting','phong','facecolor','m','edgecolor','none','FaceAlpha',.5);
pbaspect([1 1 .25]);

set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'ZTickLabel',[]);
%
newfig = figure;
a = uicontrol('Parent',newfig,'Style','slider','Position',[121,100,410,23],...
    'value',zeta, 'min',0, 'max',1);
b.Callback = @(es,ed) updateSystem(h,tf(wn^2,[1,2*(es.Value)*wn,wn^2])); 
b = uicontrol('Parent',newfig,'Style','slider','Position',[121,50,410,23],...
    'value',zeta, 'min',0, 'max',1);
c = uicontrol('Parent',newfig,'Style','slider','Position',[121,150,410,23],...
    'value',zeta, 'min',0, 'max',1);

          
          
          
bl1 = uicontrol('Parent',newfig,'Style','text','Position',[50,54,23,23],...
                'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',newfig,'Style','text','Position',[500,54,23,23],...
                'String','1','BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',newfig,'Style','text','Position',[240,25,100,23],...
                'String','Damping Ratio','BackgroundColor',bgcolor);
%             b.Callback = @(es,ed) updateSystem(h,tf(wn^2,[1,2*(es.Value)*wn,wn^2])); 



%%
cellfill = image(:,:,:,1);
gephyrin = image(:,:,:,2);
syn = image(:,:,:,3);
figure; 
camlight; lighting phong
scaleval = 5.5;
ch1 = gaussf(cellfill,[1 1 0]);
szch1 = size(ch1);
colormatrix = ones(szch1(2),szch1(1),szch1(3));
% value = (max(ch1)-min(ch1))/1.5;
value = (max(ch1)+min(ch1))/scaleval;
patch(isosurface(0:szch1(1)-1,0:szch1(2)-1,1:szch1(3),double(ch1),value,colormatrix.*double(ch1)),'parent',gca...
,'facelighting','phong','facecolor','b','edgecolor','none','FaceAlpha',.5);
hold on;

% figure;
scaleval = 5;
ch2 = gaussf(gephyrin,[1 1 0]);
szch2 = size(ch2);
colormatrix = ones(szch2(2),szch2(1),szch2(3));
value = (max(ch2)+min(ch2))/scaleval*1.4;
patch(isosurface(0:szch2(1)-1,0:szch2(2)-1,0:szch2(3)-1,double(ch2),value,colormatrix.*double(ch2)),'parent',gca,...
'facelighting','phong','facecolor','g','edgecolor','none','FaceAlpha',.5);
hold on;


scaleval = 2.5;
ab = gaussf(syn,[1 1 0]);
value = (max(ab)+min(ab))/scaleval;
abmask = ab>value;
lb = label(abmask);
lb(lb == 2) = 0;
szab = size(lb);
colormatrix = ones(szab(2),szab(1),szab(3));
patch(isosurface(0:szab(1)-1,0:szab(2)-1,0:szab(3)-1,double(lb),1,colormatrix.*double(lb)),'parent',gca,...
'facelighting','phong','facecolor','m','edgecolor','none','FaceAlpha',.5);
pbaspect([1 1 .25]);

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
