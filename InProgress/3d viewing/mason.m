filename = 'E:\SIM Worked and 3D\92217\092217SynpoTfRmChSIM\Fixed\fix092217_GFPSynpo_TfRHT_mCh_1_Reconstructed.nd2';
vd = viewer3D();
vd.loadNDFile(filename);
vd.dimension = [1 1 0.1];
vd.initialize3Dimage;
vd.makepanel;