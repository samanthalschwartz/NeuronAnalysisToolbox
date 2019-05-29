setBatchMode(false);
cumulativeIntDenArray=newArray();
directory=getDirectory("Choose a Directory");
files=getFileList(directory);
treatment=File.getParent(files[0]);
for(a=0; a<files.length; a++){
	open(directory+files[a]);
	
//Collect basic image information and set global parameters
	run("Set Measurements...", "area mean standard modal min integrated median display redirect=None decimal=3");
	slices=nSlices();
	title=getTitle();
	//Create ROI
	setSlice(nSlices);
	run("Duplicate...", "duplicate channels=1");
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
	run("3D Objects Counter", "threshold=128 slice=0 min.=20000 max.=1048576 objects");
	selectImage(2);
	close();
	selectImage(2);
	setThreshold(1,255);
	run("Create Selection");
	roiManager("Add");
	roiManager("Select", 0);
	close();
	selectImage(title);
	makeOval(25, 320, 97, 97);
	waitForUser("Find background");
	roiManager("Add");
	
	
//Measure intensity for each slice
	roiManager("Select", newArray(0,1));
	roiManager("Multi Measure");
	intDenArray=newArray(nSlices);
	intDenArraybackground=newArray(nSlices);
	frameArray=newArray(nSlices);
	baseline=getResult("RawIntDen1", 0);
	for(b=0; b<nResults; b++){
		intDen=getResult("RawIntDen1", b);
		intDenArray[b]=intDen;
		intDen=getResult("RawIntDen2", b);
		intDenArraybackground[b]=intDen;
	}
		intDenArray=Array.concat(intDenArray, intDenArraybackground);
		AreaArray = newArray(2);
		Area = getResult("Area1", 0);
		AreaArray[0] = Area;
		Area = getResult("Area2",1);
		AreaArray[1]=Area;
		intDenArray=Array.concat(intDenArray, AreaArray);
		intDenArray=Array.concat(title, intDenArray);

		cumulativeIntDenArray=Array.concat(intDenArray, cumulativeIntDenArray);
		
	selectWindow("Results");
	run("Close");
	selectImage(1);
	run("Close");
	roiManager("Deselect");
	roiManager("Delete");
	wait(1000);

}
	
for(a=0; a<files.length; a++){
	tempArray=newArray(intDenArray.length);
	tempArray = Array.slice(cumulativeIntDenArray, (a*intDenArray.length), (a+1)*(intDenArray.length));
	setResult("File", a, tempArray[0]);
	setResult("Int Den 561", a, tempArray[1]);
	setResult("Int Den 488", a, tempArray[2]);
	setResult("Int Den 640", a, tempArray[3]);
	setResult("Bckgrd Int Den 561", a, tempArray[4]);
	setResult("Bckgrd Int Den 488", a, tempArray[5]);
	setResult("Bckgrd Int Den 640", a, tempArray[6]);
	setResult("Area Neuron", a, tempArray[7]);
	setResult("Area Background", a, tempArray[8]);
}
