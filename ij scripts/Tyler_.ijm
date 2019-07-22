//Begin with a two-channel projection or single-z-plane 16-bit .tif (no z-stacks) where c:1 is the cell fill and c:2 is a punctate protein.
//There can only be one image open and the Roi Manager should be empty.


/*  It helps to have the scale set beforehand.
 *  Ideally, the scale (pixels/micron) is contained within the OME-XML metadata of the .tif.
 *  If not, uncomment the two lines below and in the prompt...
 *  Set the known "Distance in pixels:" (i.e. the number of pixels per micron).
 *  Set the "Known distance:" to "1.00".
 *  Set the "Pixel aspect ratio:" to "1.0".
 *  Set the "Unit of length:" to "micron".
 *  Click "Ok".
 */

selectImage(1)
run("Set Scale...")


//Alternatively, modify and uncomment the line below to accomplish this without any user input.

//run("Set Scale...", "distance=15.6279 known=1 pixel=1 unit=micron");


//Runs "Analyze > Set Measurements..." and ensures that only "Mean gray value" is selected.
run("Set Measurements...", "mean redirect=None decimal=3");

//Runs "Edit > Options > Input/Output..." and ensures that "Copy column headers" and "Copy row numbers" are NOT selected.
run("Input/Output...", "jpeg=85 gif=-1 file=.csv save_column save_row");


//Duplicates the channels separately, runs "Plugins > Mosaic > Utility > Background Subtractor" on each, then merges the background-subtracted channels.
//Note: Requires selection of "MOSAIC ToolSuite" in list of update sites. Also, see MosaicSuite documentation for how to set appropriate slinding-window length.
selectImage(1);
run("Duplicate...", "title=c1 duplicate channels=1-1");
selectImage(2);
run("Background Subtractor", "length=40");
selectImage(1);
run("Duplicate...", "title=c4 duplicate channels=2-2");
selectImage(3);
run("Background Subtractor", "length=25");
run("Merge Channels...", "c1=[c1] c4=[c4] create");


//Duplicates the background-subtracted cell fill channel (c:1), applies an Otsu threshold, and manipulates the mask to generate the "Cell Fill" roi.
selectImage("Composite");
run("Duplicate...", "duplicate channels=1-1");
selectImage(3);
run("Auto Threshold", "method=Otsu ignore_black ignore_white white");
run("Create Selection");
//Note: If the entire cell is NOT neatly contained within the edges of the image, you must run the next line.
run("Make Inverse");
run("Create Mask");
run("Dilate");
run("Close-");
//Note: "Analyze > Analyze Particles..." is used here simply to select only regions of the mask that meet a minimum and/or maximum pixel area and to exclude any regions that don't.
run("Analyze Particles...", "size=4-Infinity pixel show=Masks");
run("Create Selection");
roiManager("Add");
roiManager("Select", 0);
roiManager("Rename", "Cell Fill");
//Note: The following three lines are very useful whenever you finish interacting with the Roi Manager...
roiManager("Deselect");
roiManager("Show All");
roiManager("Show None");
selectWindow("Mask of Mask");
close();
selectWindow("Mask");
close();
selectImage(3);
close();


//Duplicates the background-subtracted puncta channel (c:2) and applies an Otsu threshold to generate the "Total Puncta" roi.
selectImage("Composite");
run("Duplicate...", "duplicate channels=2-2");
selectImage(3);
run("Auto Threshold", "method=Otsu ignore_black ignore_white white");
run("Create Selection");
//Note: As above, if all of the puncta are NOT neatly contained within the edges of the image, you must run the next line.
run("Make Inverse");
run("Create Mask");
//Note: As above, "Analyze > Analyze Particles..." is used here to effectively exclude any regions of the mask that fail to meet a minimum and/or maximum pixel area.
run("Analyze Particles...", "size=4-Infinity pixel show=Masks");
run("Create Selection");
roiManager("Add");
roiManager("Select", 1);
roiManager("Rename", "Total Puncta");
roiManager("Deselect");
roiManager("Show All");
roiManager("Show None");
selectWindow("Mask of Mask");
close();
selectWindow("Mask");
close();
selectImage(3);
close();


//Below is an alternative way to generate the "Total Puncta" roi that employs a more arbitrary threshold (2.5 * the mean intensity of c:2 within the Cell Fill roi).
/*
roiManager("Select", 0);
run("Clear Results");
print("\\Clear");
run("Measure");
selectWindow("Results");
c2Mean = getResult("Mean",0);
c2Thresh = (2.5*round(c2Mean));
//print(c2Mean);
//print(c2Thresh);
selectWindow("Results");
run("Close");
selectWindow("Log");
run("Close");
selectImage(3);
setThreshold(c2Thresh, 65535);
setOption("BlackBackground", true);
run("Convert to Mask");
run("Analyze Particles...", "size=4-Infinity pixel show=Masks");
run("Create Selection");
roiManager("Add");
roiManager("Select", 1);
roiManager("Rename", "Total Puncta");
roiManager("Deselect");
roiManager("Show All");
roiManager("Show None");
selectWindow("Mask of Composite-1");
run("Close");
selectImage(3);
close();
*/


selectImage("Composite");
run("Split Channels");


//Creates "Mask1" from the "Cell Fill" roi on c:1.
selectWindow("C1-Composite");
roiManager("Select", 0);
setBackgroundColor(0, 0, 0);
run("Clear Outside");
run("Create Mask");
roiManager("Deselect");
roiManager("Show All");
roiManager("Show None");
selectWindow("C1-Composite");
close();
selectWindow("Mask");
rename("Mask1");


//Creates "Mask2" from the "Total Puncta" roi on c:2.
selectWindow("C2-Composite");
roiManager("Select", 1);
setBackgroundColor(0, 0, 0);
run("Clear Outside");
run("Create Mask");
roiManager("Deselect");
roiManager("Show All");
roiManager("Show None");
selectWindow("C2-Composite");
close();
selectWindow("Mask");
rename("Mask2");


//Labels all objects in "Mask2" ("Plugins > MorphoLibJ > Binary Images > Connected Components Labeling"), generating "Mask2-lbl".
//For a given object, each pixel is assigned the grayscale value of that object's label.
//i.e. for an 8-bit output objects are labeled 1-255, and for a 16-bit output objects are labeled 1-65535.
//Note: Requires selection of "IJPB-plugins" in list of update sites.
selectWindow("Mask2");
run("Connected Components Labeling", "connectivity=8 type=[16 bits]");


//Assigns all nonzero pixels in "Mask1" the grayscale value of 1.
selectWindow("Mask1");
roiManager("Select", 0);
run("Set...", "value=1");


//Multiplies "Mask2-lbl" by "Mask1" to generate a new mask of their intersection ("Result of Mask2-lbl").
imageCalculator("Multiply create", "Mask2-lbl","Mask1");


//Generates a string array of the unique labels present in "Result of Mask2-lbl", prints the contents of this array to the "Log" window, and then copies this output to the clipboard.
//Note: Requires "IJPB-plugins"...
selectWindow("Result of Mask2-lbl");
run("Region Morphometry");
IJ.renameResults("Result-Morphometry", "Results");
A = newArray(nResults);
for (i=0; i<nResults; i++){
  A[i] = getResultString("Label",i);
}
Array.print(A);
selectWindow("Result of Mask2-lbl");
close();
selectWindow("Mask1");
close();
selectWindow("Mask2");
close();
selectWindow("Results");
run("Close");
String.copy(getInfo("log"));
selectWindow("Log");
run("Close");


//Returns to "Mask2-lbl", prompts user to enter a subset of labels/objects to keep, and generates a new mask ("Mask2-lbl-keepLabels").
//Note: Requires "IJPB-plugins"...
//Note: The script hangs here waiting for user input. Paste and click "Ok" (Ctrl + v, Return/Enter).
selectWindow("Mask2-lbl");
run("Select Label(s)");


//Prompts user to enter a subset of labels/objects to modify and a "Final Value" to assign to their pixels.
//Note: Requires "IJPB-plugins"...
//Note: The script hangs here again. Paste, manually enter "65535", and click "Ok" (Ctrl + v, Tab, 65535, Return/Enter).
selectWindow("Mask2-lbl-keepLabels");
run("Replace/Remove Label(s)");


//Generates the "Puncta" roi. This is the subset of puncta from "Total Puncta" that overlap with "Cell Fill".
selectWindow("Mask2-lbl-keepLabels");
run("8-bit");
run("Create Selection");
run("Make Inverse");
roiManager("Add");
roiManager("Select", 2);
roiManager("Rename", "Puncta");
roiManager("Deselect");
roiManager("Show All");
roiManager("Show None");
selectWindow("Mask2-lbl-keepLabels");
close();
selectWindow("Mask2-lbl");
close();


//Returns to the original two-channel projection or single-z-plane 16-bit .tif, duplicates c:2, applies a LUT, applies the "Puncta" roi, and clears outside.
selectImage(1);
run("Duplicate...", "title=c4 duplicate channels=2-2");
selectImage(2);
run("Grays");
roiManager("Select", 2);
run("Clear Outside");
roiManager("Deselect");
roiManager("Show All");
roiManager("Show None");


//Creates a mask from "Puncta" to analyze, then returns a summary table with puncta "Count" and "Average Size" (average area).
selectImage(2)
roiManager("Select", 2);
run("Create Mask");
run("Analyze Particles...", "size=4-Infinity pixel show=Nothing summarize");
selectWindow("Mask");
close();
roiManager("Deselect");
roiManager("Show All");
roiManager("Show None");


//Returns to the original two-channel projection or single-z-plane 16-bit .tif, duplicates c:1, applies a LUT, applies the "Cell Fill" roi, and clears outside.
selectImage(1);
run("Duplicate...", "title=c1 duplicate channels=1-1");
selectImage(3);
run("Red");
roiManager("Select", 0);
run("Clear Outside");
roiManager("Deselect");
roiManager("Show All");
roiManager("Show None");


//Creates a new composite image from the segmented channels for visual comparison / validation against the original.
run("Merge Channels...", "c1=c1 c4=c4 create");
selectWindow("Composite");
Stack.setSlice(1);
run("Enhance Contrast", "saturated=0.35");
run("Next Slice [>]");
run("Enhance Contrast", "saturated=0.35");


selectImage(1);
Stack.setDisplayMode("composite");
Stack.setSlice(1);
run("Previous Slice [<]");
run("Red");
run("Next Slice [>]");
run("Grays");
Stack.setSlice(1);
run("Previous Slice [<]");
run("Enhance Contrast", "saturated=0.35");
run("Next Slice [>]");
run("Enhance Contrast", "saturated=0.35");


selectImage(2);