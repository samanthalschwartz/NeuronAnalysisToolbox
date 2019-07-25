dir = getDirectory("Choose Source Directory ");
      dir2 = getDirectory("Chsoose Target Directory ");
      list = getFileList(dir);
      setBatchMode(true);
      n = list.length;
      basename = "GFP-gephOlig_star568antiGabag2-AF647antiGFP-nolight";
for (i=1; i<n; i++) {
green = basename+"_w1488_s"+i+".TIF";
red = basename+"_w2561_s"+i+".TIF";
white = basename+"_w3640_s"+i+".TIF";
savename = dir2 + basename+"_s"+i+".TIF";
savename_zproj = dir2 + basename+"_Zprojected_s"+i+".TIF";
open(dir+green);
open(dir+red);
open(dir+white);
run("Merge Channels...", "c1=["+red+"] c2=["+green+"] c4=["+white+"]");
selectWindow("Composite");
saveAs("Tiff", savename);
run("Z Project...", "projection=[Max Intensity]");
saveAs("Tiff", savename_zproj);
close();
}