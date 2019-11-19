Dialog.create("Crop ROIs");
Dialog.addChoice("Type:", newArray("Load image and ROIs", "Run from open Image and ROIs")) ;
Dialog.show();
type = Dialog.getChoice();

if (type == "Load image and ROIs") {
id_roi = File.openDialog("Choose an ImageJ ROI file");
id_file = File.openDialog("Choose an Image to Crop ROIs");
run("Open",id);
roiManager("Open",id);
}
savedir = getDirectory("Choose Save Directory");

for (i=0; i<roiManager("count"); ++i) {
    	roiManager("Select", i);
    	run("Duplicate...", "title=tempim duplicate");
		selectWindow("tempim");
        saveAs("Tiff", savedir + "roi_" + d2s(i,0) + ".tif");
        close();
    }


