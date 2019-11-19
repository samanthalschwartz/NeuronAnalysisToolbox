dir = getDirectory("Choose Source Directory");
list = getFileList(dir);
for (g=0; g<list.length; g++) {
	path = dir+list[g];
	open(path);
	run("Split Channels");
	selectWindow("C1-" + list[g]);
	run("Duplicate...", "title=duplicateCh1");
	run("Combine...", "stack1=C1-"+list[g]+" stack2=C2-"+list[g]);
	rename("Ch1-Ch2"); 
	saveAs("Tiff",dir+list[g]+"c1c2");
	//run("Combine...", "stack1=duplicateCh1 stack2=C3-"+list[g]);A
	//rename("Ch1-Ch3"); 
	//saveAs("Tiff",dir+list[g]+"c1c3");
	close();
	close();
}

