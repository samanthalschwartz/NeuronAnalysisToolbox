uiopen('C:\Users\sammy\Dropbox\Sam Kennedy Lab\Projects\GephyrinIntrabody-Olig Clustering\GABA-Geph Quant\+light_Cry2Olig488_RIM561_Geph647\+light_Cry2Olig488_RIM561_Geph647_cells8-11.tif',1);
%make GephyrinIB channel mask
ch2 = squeeze(image(:,:,:,2));
ch2l = GeneralAnalysis.imgLaplaceCutoff(ch2);
ch2lm = GeneralAnalysis.imgThreshold_fixedUserInput(ch2l);
ch2m = ch2lm.*max(ch2(:));
ch2mi = single(cat(4,ch2m,ch2));
MIJ.createImage('ch2overlay.tif',ch2mi,true)
MIJ.run("Make Composite", "display=Composite");
MIJ.run('Red');
MIJ.run("Next Slice [>]");
MIJ.run("Grays");
MIJ.close

ch2_lesssens_lm = GeneralAnalysis.imgThreshold_fixedUserInput(ch2l);
ch2_lesssens_m = ch2_lesssens_lm.*max(ch2(:))*2;
ch2_lesssens_mi = single(cat(4,ch2_lesssens_m,ch2));
MIJ.createImage('ch2overlay_lesssens.tif',ch2_lesssens_mi,true)
MIJ.run("Make Composite", "display=Composite");
MIJ.run('Red');
MIJ.run("Next Slice [>]");
MIJ.run("Grays");
MIJ.close




% make anti-Gephyrin ab channel mask
ch3 = squeeze(image(:,:,:,3));
ch3l = GeneralAnalysis.imgLaplaceCutoff(ch3);
ch3lm = GeneralAnalysis.imgThreshold_fixedUserInput(ch3l);
ch3m = ch3lm.*max(ch3(:));
ch3mi = single(cat(4,ch3m,ch3));
MIJ.createImage('ch3overlay.tif',ch3mi,true)
MIJ.run("Make Composite", "display=Composite");
MIJ.run('Red');
MIJ.run("Next Slice [>]");
MIJ.run("Grays");

newovlay = single(cat(4,ch2,ch3,ch2_lesssens_m,ch3m));
MIJ.createImage('ch23maskoverlay.tif',newovlay,true);
MIJ.run("Make Composite", "display=Composite");
MIJ.run('Green');
MIJ.run("Next Slice [>]");
MIJ.run("Grays");
MIJ.run("Next Slice [>]");
MIJ.run("Red");
MIJ.run("Next Slice [>]");
MIJ.run("Blue");




