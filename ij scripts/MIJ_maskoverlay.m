function MIJ_maskoverlay(mask,greyim,filename)
if nargin<3
    filename = 'ch_overlay.tif';
end
mask = mask.*max(greyim(:)*10);
ovly = double(cat(4,mask,greyim,mask.*0));
MIJ.createColor(ovly,true)
MIJ.run("Make Composite", "display=Composite");
MIJ.run('Red');
MIJ.run("Next Slice [>]");
MIJ.run("Grays");
end