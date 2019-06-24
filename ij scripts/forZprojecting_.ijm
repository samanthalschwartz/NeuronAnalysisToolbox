for (i=8; i>0; i--) {
filename= "gephOlig-488light-mchtimelapse";
saveAs("Tiff", "I:/Projects/Project Cry2Olig-Gephyrin/180724 GephIntraOlig-GabamCh/coverslip1-olig-mChGaba/process/"+filename+"_s"+i+".tif");
run("Z Project...", "projection=[Max Intensity] all");
saveAs("Tiff", "I:/Projects/Project Cry2Olig-Gephyrin/180724 GephIntraOlig-GabamCh/coverslip1-olig-mChGaba/process/"+filename+"_Zproj_s"+i+".tif");
close();
close();
}