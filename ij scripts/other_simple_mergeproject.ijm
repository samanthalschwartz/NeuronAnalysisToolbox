topdir= "I:/Projects/Project Cry2Olig-Gephyrin/180724 GephIntraOlig-GabamCh/coverslip1-olig-mChGaba/process/";
dir1 = "gephOlig-488light__2colors/";
dir2 = "gephOlig-488lighagain__2colors/";
dir3 = "gephOlig-488lighagain2__2colors/";
savedir= topdir + "combined/";
filename1 = "gephOlig-488light_Zproj_s";
filename2 = "gephOlig-488lighagain_Zproj_s";
filename3 = "gephOlig-488lighagain2_Zproj_s";
for (i=1; i<9; i++) {
open(topdir+dir1+filename1+i+".tif");
open(topdir+dir2+filename2+i+".tif");
open(topdir+dir3+filename3+i+".tif");
run("Concatenate...", "  title=gephOlig+488_2color_s"+i+" open image1="+filename1+i+".tif image2="+filename2+i+".tif image3="+filename3+i+".tif");
saveAs("Tiff", savedir+"gephOlig-488light_2color_Zproj_s"+i);
close();
}