outputarry = Table.getColumn("frame");
newarray = newArray(outputarry.length);
Array.getStatistics(outputarry, min, max, mean, stdDev);
for (g=0; g<outputarry.length; g++) {
	curval = outputarry[g];
	newval = max-curval+1;
	setResult("frame",g,newval);
}