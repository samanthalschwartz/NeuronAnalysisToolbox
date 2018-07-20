s = load('C:\Users\schwsama\Documents\Data\zapERtrap\cell3_AshleyFile.mat');
soma = s.aa.cellFill.image;
d = load('C:\Users\schwsama\Documents\Data\zapERtrap\cell3dendrites_AshleyFile.mat');
dendrites = d.aa.cellFill.image;

soma1 = squeeze(single(soma(:,:,1)));
dend1 = squeeze(single(dendrites(:,:,1)));

tic
[newimage, ccpeak] = GeneralAnalysis.stitch2images(soma1,dend1);
toc

soma1 = single(soma);
dend1 = single(dendrites);

tic
[stitchmovie1,ccpeak] = GeneralAnalysis.stitch2movies(soma1,dend1);
toc


soma = single(s.aa.TfR.image);
dendrites = single(d.aa.TfR.image);

tic
[stitchmovie2,ccpeak] = GeneralAnalysis.stitch2movies(soma,dendrites,ccpeak);
toc


soma2 = single(s.aa.surfaceCargo.image);
dendrites2 = single(d.aa.surfaceCargo.image);

tic
[stitchmovie3,ccpeak] = GeneralAnalysis.stitch2movies(soma2,dendrites2,ccpeak);
toc
