imageSize=[256,512];
psfSigma=[1.1, 1.1];
nParticles=100;
nFrames=10;

hss=HSSim(imageSize, psfSigma);
[im2d, loc2d]=hss.simulate2DImage(nParticles, nFrames);
