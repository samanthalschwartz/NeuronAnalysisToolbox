filename = 'F:\SIM\180704 Coverslip3\Coverslip3_Reconstructed-cell1pos1.nd2';
image = ndFileloader(filename);
% image = ndFileloader();
geph = image(:,:,:,1);
gfp = image(:,:,:,2);
gaba = image(:,:,:,3);

out = GeneralAnalysis.imgLaplaceCutoff(gaba,[1 1 1],[1 1 1]);