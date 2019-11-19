//load the data
path = File.openDialog("Select a SR Results .csv file");
 dir = File.getParent(path);
  name = File.getName(path);
  print("Path:", path);
  print("Name:", name);
  print("Directory:", dir);

open(path);
outputarry = Table.getColumn("frame");
newarray = newArray(outputarry.length);
Array.getStatistics(outputarry, min, max, mean, stdDev);
for (g=0; g<outputarry.length; g++) {
	curval = outputarry[g];
	Table.set("frame", g, max-curval+1);
	//newarray[g] = max-curval+1;
}
//Table.setColumn("frame",newarry);
//save the  data
saveAs("Measurements", dir + "/flipped"+name);