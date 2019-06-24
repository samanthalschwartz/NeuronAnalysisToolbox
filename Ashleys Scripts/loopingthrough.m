close all; clear all;
[files, pathname]  = uigetfile('G:\zapERtrap\Raw Data\GLOBAL RELEASE\GluA1\060118','multiselect','on');

for ff = 1:numel(files)
    load(fullfile(pathname,files));
    aa.cellFill.ROI_trim = [];
    aa.cellFill.trim_rawimage();
    aa.surfaceCargo.ROI_trim = aa.cellFill.ROI_trim;
    aa.surfaceCargo.trim_rawimage;
    aa.cleanedcargomask = []; %just make sure to reset cleaned cargomask
    save(fullfile(pathname,files),'aa', '-v7.3');
    clear aa;
end




wb = waitbar(0);
for ff = 2:numel(files)
    load(fullfile(pathname,files{ff}));
    [img_out_cf,sv_arr] = GeneralAnalysis.timedriftCorrect(aa.cellFill.image);
    aa.cellFill.setimage(img_out_cf);
    img_out_sc = GeneralAnalysis.applydriftCorrect(aa.surfaceCargo.image,sv_arr);
    aa.surfaceCargo.setimage(img_out_sc);
    
    save(fullfile(pathname,files{ff}),'aa', '-v7.3');
    waitbar(ff/numel(files),wb)
    clear aa;
end
