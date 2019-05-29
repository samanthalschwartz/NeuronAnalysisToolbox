//Collect basic image information and set global parameters
run("Set Measurements...", "mean decimal=3");
slices=nSlices();
title=getTitle();
//Create ROI
setSlice(nSlices);
run("Duplicate...", "title=&title");
run("Gaussian Blur...", "sigma=2 slice");
setAutoThreshold("Mean dark");
run("Create Selection");
setBackgroundColor(0,0,0);
setForegroundColor(255, 255, 255);
run("Fill", "slice");
run("Clear Outside");
setMinAndMax(0, 0);
run("8-bit");
run("Select None");
run("3D Objects Counter", "threshold=128 slice=0 min.=1000 max.=1048576 objects");
selectImage(2);
close();
selectImage(2);
setThreshold(1,255);
run("Create Selection");
roiManager("Add");
roiManager("Select", 0);
close();