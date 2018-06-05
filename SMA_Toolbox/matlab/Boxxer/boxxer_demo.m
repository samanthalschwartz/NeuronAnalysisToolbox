%% Make objects


sigma=[1.0 1.6 2.0 2.6 1.0 1.0 2.0 4.0;  %sigma_x in pixels
       1.0 1.6 2.0 2.6 2.0 4.0 1.0 1.0]; %sigma_y in pixels 
imsize=[256, 512]; %[Y X]


%% Make sim image
Nparticles=100;
Nframes=4;
Nscales=size(sigma,2);
im=zeros([imsize Nframes],'single');
tic;
for i=1:Nscales
    hss=HSSim(flip(imsize), sigma(:,i)');
    [nim,npos]=hss.simulate2DImage(Nparticles,Nframes);
    im=im+nim;
    bd=Nparticles*Nframes;
    orig_pos(1:3,(i-1)*bd+1:i*bd)=npos';
    orig_pos(4,(i-1)*bd+1:i*bd)=i;
end
simT=toc;
fprintf('Simulation time: %.3f s',simT);

%% Get FAP image
% load('/home/mjo/LidkeLab/Data/FAP/FAP_3-3V-2013-7-3-12-33-9.mat')
% sz=size(sequence,1);
% psf=[1.3, 1.3];
% im=single(sequence(:,:,1:100));
% gain=0.04;
% bg=118;
% im=gain*(im-bg);

%% Filter image
dip_image(im)
diptruesize(400);
tic;
boxxer=Boxxer2D(imsize, sigma); %put this is backwards since our images will be [Y X]
logfim=boxxer.filterLoG(im);
dogfim=boxxer.filterDoG(im);
gfim=boxxer.filterGauss(im);
slogfim=boxxer.filterScaledLoG(im);
sdogfim=boxxer.filterScaledDoG(im);
[smaxima, smax_vals] = boxxer.scaleSpaceDoGMaxima(im, 5, 3);
[fsmaxima, fsmax_vals, filter] = Boxxer2D.filterMaxima(smaxima, smax_vals);
overlayIm= boxxer.plotScaleMaxima(sdogfim,smaxima,filter);
toc
dip_image(overlayIm)
% [maxima, max_vals]=boxxer.enumerateMaxima(logfim,5);
% [fmaxima, fmax_vals, filter]=Boxxer2D.filterMaxima(maxima, max_vals);
% boxxer.plotMaxima(logfim, maxima, max_vals, filter)
% boxxer.plotMaxima(im, maxima, max_vals, filter)
% 
% bc=boxxer.generateBoxCoords(fmaxima,Nframes);
% closeboxim=boxxer.plotBoxCoordsDIP(logfim, bc, fmaxima, fmax_vals)
