dir = getDirectory("Choose Source Directory ");
print(dir)
base = "nolight_Cry2Olig488_RIM561_Gaba647_cell2";
ch2 = "_w1488";
ch1 = "_w2561";
ch3 = "_w3640";
stagestr = "_s";
for (i=1; i<3; i++) {
ch1str = base + ch1 + stagestr + i + ".TIF";
ch2str = base + ch2 + stagestr + i + ".TIF";
ch3str = base + ch3 + stagestr + i + ".TIF";
open(dir+ch1str);
open(dir+ch2str);
open(dir+ch3str);
run("Merge Channels...", "c1=["+ch1str+"] c2=["+ch2str+"] c4=["+ch3str+"]");
save(dir+base+i);
run("Z Project...", "projection=[Max Intensity]");
save(dir+base+"_"+i+"maxZproject");
close();
close();
}