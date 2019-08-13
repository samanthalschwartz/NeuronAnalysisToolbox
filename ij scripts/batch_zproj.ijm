run("Bio-Formats Macro Extensions");
 dir = getDirectory("Choose Source Directory");
      savedir = getDirectory("Choose a Save Directory");
      print("Data Directory: " + dir);
      print("Save Directory: " + savedir);
      list = getFileList(dir);
      alllen = list.length;
      filelist = newArray(0);
      for (f=0; f<alllen; f++){
		if (endsWith(list[f],".tif")) {
		filelist = Array.concat(filelist,list[f]);   
		}   	
      }
		len = filelist.length;
      print(len);
      print(filelist[0]);
		for (g=0; g<len; g++) {
      	id = dir + filelist[g];
      	makeZproj(id,savedir);
      	close();
      	}
      	print("Finished the thing");

function makeZproj(id,savedir){
	run("Bio-Formats Importer", "open=["+id+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
	basename = File.getName(id);
	basename = substring(basename,0,lastIndexOf(basename, "."));
	currZprojfilename = savedir+basename+"_Zproj_position";
	run("Z Project...", "projection=[Max Intensity] all");
	run("Make Composite");
	Ext.setId(id);
	run("Green");
	run("Next Slice [>]");
	run("Red");
	run("Next Slice [>]");
	run("Grays");
	saveAs("Tiff",currZprojfilename);
	print("Making File: "+currZprojfilename);
	close();
}