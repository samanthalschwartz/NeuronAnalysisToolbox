filepath = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\071117\Gephyrin_Abeta_Bassoon\Geph488_Abeta561_Bassoon647_003_Reconstructed.nd2';
s  = SIM();
s.loadNDfile(filepath);
s.make_masks();
s.make_distancemasks();
s.calculateDensity_abetaINch1();
s.calculateDensity_abetaINch2();