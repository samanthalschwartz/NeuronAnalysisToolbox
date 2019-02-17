dir = getDirectory("Choose a Directory ");
list = getFileList(dir);
startstring = "open=";
for (i=0; i<list.length; i++) {
	v = dir+list[i];
	if (endsWith(v, ".csv")) {
run("XY Coordinates... ","open="+v);
roiManager("add")
	}
      }