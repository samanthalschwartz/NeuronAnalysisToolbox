datadir = "D:/WDKennedyLabHDDBackup/Projects/Project Cry2Olig-Gephyrin/HT GABA/190424 GabaA2-HaloTag/cov2_Cry2_HTGaba+JF565_+488/";
fullfile = datadir+"cov2_Cry2_HTGaba+JF565_+488_w2488_s1_t10.TIF";
strmatch = "cov2_Cry2_HTGaba+JF565_+488_w248_s1_t<1,6,11,16,21,26,31,36,41>.TIF";
savename = "cov2_Cry2_HTGaba+JF565_+488_ch488_s1_timeseries";
fullstrmatch = datadir+strmatch;
run("Bio-Formats Importer", "open=["+fullfile+"] color_mode=Default group_files rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT name=["+fullstrmatch+"]");
saveAs("Tiff", datadir+savename+".TIF");
run("Z Project...", "projection=[Max Intensity] all");
saveAs("Tiff", datadir+savename+"_zproj.TIF");
close(); 
close();