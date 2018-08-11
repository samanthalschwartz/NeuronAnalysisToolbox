% stitching images together
% Z:\Sam\MJK_zapERtrap_for_sam\040318_globalrelease_NL1 data
ga = GeneralAnalysis();
datadir = 'E:\Sam\Data\MJK_zapERtrap_for_sam\040318_globalrelease_NL1';
cargostr = 'cell1*_fr.tif';
anchorstr = 'cell1*_g.tif';
cellfillstr = 'cell1*_r.tif';

cellfillfiles = dir(fullfile(datadir,cellfillstr));
anchorfiles = dir(fullfile(datadir,anchorstr));
cargofiles = dir(fullfile(datadir,cargostr));



cellfills = cell(3,1);
anchors = cell(3,1);
cargos = cell(3,1);
tic
for ff = 1:numel(cellfillfiles)
    temp = ga.loadtiff_1ch(fullfile(datadir,cellfillfiles(ff).name));
    [img_out,sv_arr] = ga.timedriftCorrect(temp);
    cellfills{ff} = img_out;
    clear temp; clear img_out;
    
    temp = ga.loadtiff_1ch(fullfile(datadir,anchorfiles(ff).name));
    img_out = ga.applydriftCorrect(temp,sv_arr);
    anchors{ff} = img_out;
    clear temp; clear img_out;
    
    temp = ga.loadtiff_1ch(fullfile(datadir,cargofiles(ff).name));
    img_out = ga.applydriftCorrect(temp,sv_arr);
    cargos{ff} = img_out;
    clear temp; clear img_out;
end
toc
cellbody = 1;
dend1 = 2;
dend2 = 3;
[stitchmovie,ccpeak_out] = ga.stitch2movies(cellfills{dend1},cellfills{cellbody});