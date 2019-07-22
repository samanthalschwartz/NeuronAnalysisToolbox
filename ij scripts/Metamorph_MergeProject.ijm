Dialog.create("Metamorph Data Batch Convert");
Dialog.addChoice("Type:", newArray("Process Individual .nd file", "Process All .nd files in Folder")) ;
Dialog.show();
type = Dialog.getChoice();
run("Bio-Formats Macro Extensions");
print(type);
if (type == "Process All .nd files in Folder") {
      dir = getDirectory("Choose Source Directory ");
      print(dir);
      list = getFileList(dir);
      alllen = list.length;
      filelist = newArray(0);
      for (f=0; f<alllen; f++){
		if (endsWith(list[f],".nd")) {
		filelist = Array.concat(filelist,list[f]);   
		}   	
      }
      len = filelist.length;
      print(len);
      for (f=0; f<len; f++) {
      	id = filelist[f];
      	print(id);
      	processNDfiles(id);
      }
}
else {
      id = File.openDialog("Choose a file");
      processNDfiles(id);
};
     
	function processNDfiles(id){
		Ext.setId(id);
		Ext.getSeriesCount(seriesCount);
		print("Num Series = " +seriesCount);
		basename = File.nameWithoutExtension;
		basepath = File.directory;
		print(basepath);
		for (s=0; s<seriesCount; s++) {
			run("Bio-Formats Importer", "open=["+id+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series"+d2s(s,0)); 
			inputId = getImageID();
			print("id = " + inputId);
			currfilename = basepath+basename+"_position_"+d2s(s+1,0);
			print(currfilename);
			saveAs("Tiff", currfilename);
			currZprojfilename = basepath+basename+"_Zproj_position_"+d2s(s+1,0);
			run("Z Project...", "projection=[Max Intensity] all");
			saveAs("Tiff",currZprojfilename);
			close();
			close();
		}
	}