Dialog.create("RGB Batch Convert");
  Dialog.addString("Red Suffix:", "561");
  Dialog.addString("Green Suffix:", "488");
  Dialog.addString("Blue Suffix:", "640");
  Dialog.addCheckbox("Open as Stack", false);
  Dialog.show();
  redSuffix = Dialog.getString() + ".";
  greenSuffix = Dialog.getString() + ".";
  blueSuffix = Dialog.getString() + ".";
 openAsStack = Dialog.getCheckbox();
  if (openAsStack)
      openImagesAsStack();
  else
      batchConvert();
  exit;
  
    function openImagesAsStack() {
      dir = getDirectory("Choose Source Directory ");
      dir2 = getDirectory("Choose Target Directory ");
      list = getFileList(dir);
      setBatchMode(true);
      n = list.length;
      if ((n%3)!=0)
         exit("The number of files must be a multiple of 3");
      stack = 0;
      first = 0;
      for (i=0; i<n/3; i++) {
          showProgress(i+1, n/3);
          red="?"; green="?"; blue="?";
          for (j=first; j<first+3; j++) {
              if (indexOf(list[j], redSuffix)!=-1)
                  red = list[j];
              if (indexOf(list[j], greenSuffix)!=-1)
                  green = list[j];
              if (indexOf(list[j], blueSuffix)!=-1)
                  blue = list[j];
          }
          open(dir+red);
          open(dir+green);
          open(dir+blue);
        run("Merge Channels...", "c1=["+red+"] c2=["+green+"] c4=["+blue+"]");
        run("Z Project...", "projection=[Max Intensity]");
		save(dir2+red);
		close()
      }
        
  function batchConvert() {
      dir1 = getDirectory("Choose Source Directory ");
      dir2 = getDirectory("Choose Destination Directory ");
      list = getFileList(dir1);
      setBatchMode(true);
      n = list.length;
      if ((n%3)!=0)
         exit("The number of files must be a multiple of 3");
      stack = 0;
      first = 0;
      for (i=0; i<n/3; i++) {
          showProgress(i+1, n/3);
          red="?"; green="?"; blue="?";
          for (j=first; j<first+3; j++) {
              if (indexOf(list[j], redSuffix)!=-1)
                  red = list[j];
              if (indexOf(list[j], greenSuffix)!=-1)
                  green = list[j];
              if (indexOf(list[j], blueSuffix)!=-1)
                  blue = list[j];
          }
          open(dir1 +red);
          open(dir1 +green);
          open(dir1 +blue);
         run("Merge Channels...", "c1=["+red+"] c2=["+green+"] c4=["+blue+"]");
         run("Z Project...", "projection=[Max Intensity]");
         setSlice(1);
         run("Enhance Contrast", "saturated=0.35");
         setSlice(2);
		run("Enhance Contrast", "saturated=0.35");
		setSlice(3);
		run("Enhance Contrast", "saturated=0.35");
          index = indexOf(red, redSuffix);
          name = substring(red, 0, index);
          saveAs("tiff", dir2+name);
          first += 3;
      }
  }    
