% convert gain tiff
gainfilepath = 'G:\Sam\Data\180522 SRFF pilot\180522 01_gainAssesment_XY1527001456_Z0_T00_C1.tif';
sequence = GeneralAnalysis.loadtiff_1ch(gainfilepath);
save(gainfilepath(1:end-4),'sequence');
% convert background tiff
bkfilepath = 'G:\Sam\Data\180522 SRFF pilot\180522 03_background_XY1527002044_Z0_T00_C0.tif';
sequence = GeneralAnalysis.loadtiff_1ch(bkfilepath);
save([bkfilepath(1:end-4),'.mat'],'sequence');

%convert bead files

beadch1file = 'G:\Sam\Data\180522 SRFF pilot\180522 06_beadChrom - 4_XY1527005934_Z0_T0_C1.tif';
sequence = GeneralAnalysis.loadtiff_1ch(beadch1file);
save([beadch1file(1:end-4) '.mat'],'sequence');


beadch2file = 'G:\Sam\Data\180522 SRFF pilot\180522 06_beadChrom - 4_XY1527005934_Z0_T0_C2.tif';
sequence = GeneralAnalysis.loadtiff_1ch(beadch2file);
save([beadch2file(1:end-4) '.mat'],'sequence');

ch1 = load([beadch1file(1:end-4) '.mat']);
ch2 = load([beadch2file(1:end-4) '.mat']);

savedir = 'G:\Sam\Data\180522 SRFF pilot\beadcomposite';
sequence = [ch1.sequence,ch2.sequence];
save(savedir,'sequence')
%%

ra = RegistrationAnalysis();
ra.gui