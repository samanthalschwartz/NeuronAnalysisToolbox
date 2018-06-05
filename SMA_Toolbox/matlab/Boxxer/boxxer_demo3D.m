%% Make objects


%% Make sim image
sigma=[1.3, 2.0, 2.8;    %L
       1.1, 1.4, 1.5;    %Y
       1.0, 1.2, 1.3  ]; %X
szL=125;
szY=64;
szX=32;
N=50;
Nframes=2;
Nscales = size(sigma,2);
imsize = [szL szY szX];
tic;
im=zeros([imsize Nframes],'single');
for i=1:Nscales
    fprintf('Scale: %.2f, %.2f %.2f\n',sigma(1,i), sigma(2,i), sigma(3,i));
    hss=HSSim(imsize([2 1 3]), sigma(:,i));
    [nim,npos]=hss.simulate3DImage(Nparticles,Nframes);
    im=im+nim;
end
simT=toc;
fprintf('Size [LYXT] %i x %i x %i x %i \nSimulation time: %.3f s\n', szL, szY, szX, Nframes, simT);
%% Get FAP image
% load('/home/mjo/LidkeLab/Data/FAP/FAP_3-3V-2013-7-3-12-33-9.mat')
% sz=size(sequence,1);
% psf=[1.3, 1.3];
% im=single(sequence(:,:,1:100));
% gain=0.04;
% bg=118;
% im=gain*(im-bg);

%% Filter image
DIPim=dip_image(im)
boxxer=Boxxer3D(imsize, sigma(:,1));
logfim=boxxer.filterDoG(im,sigma(:,1));
[maxima, max_vals]=boxxer.enumerateImageMaxima(logfim,5);
[fmaxima, fmax_vals]=boxxer.filterMaxima(maxima, max_vals);
boxxer.plotMaxima(logfim, fmaxima);
boxxer.plotMaxima(im, fmaxima);
boxxer.checkMaxima(logfim, maxima, max_vals);
boxxer.checkMaxima(logfim, fmaxima, fmax_vals);

bc=boxxer.generateBoxCoords(fmaxima);
boxim=boxxer.plotBoxCoordsDIP(logfim, bc);
