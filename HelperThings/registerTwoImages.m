function [reg_im,shift_im] = registerTwoImages()
[filename_reg, pathname_reg] = uigetfile('*.*','Select file to register');
[filename_shift, pathname_shift] = uigetfile(fullfile(pathname_reg,'*.*'),'Select file to shift');
ch1im = loadtiff(fullfile(pathname_reg,filename_reg));
ch2im = loadtiff(fullfile(pathname_shift,filename_shift));
[reg_im,sv_arr] = GeneralAnalysis.timedriftCorrect_parfor(ch1im);
shift_im =  GeneralAnalysis.applydriftCorrect(dip_image(ch2im),sv_arr);
end