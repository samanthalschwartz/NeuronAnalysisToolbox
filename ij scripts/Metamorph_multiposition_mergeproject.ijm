id = File.openDialog("Choose an .nd file");
savedir = getDirectory("Choose a Save Directory");
savestr= "cell1cry2olig_round1_with488_postTimeSeries";
numstagepos = 18;

for (s=1; s<(numstagepos+1); s++){
run("Bio-Formats Importer", "open=["+id+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+d2s(s,0));
saveAs("Tiff", savedir+savestr+"_position"+d2s(s,0)+".tif");
print("saving "+savedir+savestr+"_position"+d2s(s,0)+".tif");
run("Z Project...", "projection=[Max Intensity] all");
saveAs("Tiff", savedir+savestr+"_position"+d2s(s,0)+"_zproject.tif");
close();
close();
}