namestr = "NL1-halo_130minpostrelease_50nM660olddy+newdye"
for (y=1; y<=9; y++){
run("Merge Channels...", "c2="+namestr+"_w1488_s"+y+"_maxproj.tiff c1="+namestr+"_w2561_s"+y+"_maxproj.tiff c4="+namestr+"_w3640_s"+y+"_maxproj.tiff create");
saveAs("Tiff", "G:/FromMicroscopeComputer/190108 NL1Halo_2color early-late/overlays/"+namestr+"_s"+y+".tif");
close();
}