filename = 'D:\SIM\180704 Coverslip3\Coverslip3_Reconstructed-cell1pos1.nd2';
image = ndFileloader(filename);

geph = image(:,:,:,1);
gfp = image(:,:,:,2);
gaba = image(:,:,:,3);