s = load('C:\Users\schwsama\Documents\Data\zapERtrap\cell3_AshleyFile.mat');
soma = s.aa.cellFill.image;
d = load('C:\Users\schwsama\Documents\Data\zapERtrap\cell3dendrites_AshleyFile.mat');
dendrites = d.aa.cellFill.image;

soma1 = squeeze(single(soma(:,:,1)));
dend1 = squeeze(single(dendrites(:,:,1)));

tic
GeneralAnalysis.stitch2images(soma1,dend1);
toc

soma1 = single(soma);
dend1 = single(dendrites);

tic
GeneralAnalysis.stitch2movies(soma1,dend1);
toc