day1.datadir = 'F:\Pharm Department Channel Alignment\180625 channel alignment test';
day1.filenamec0 = 'zstack_XY1529948111_Z0_T0_C0.tif';
day1.filenamec1 = 'zstack_XY1529948111_Z0_T0_C1.tif';
day1.filenamec2 = 'zstack_XY1529948111_Z0_T0_C2.tif';
day1.filenamec3 = 'zstack_XY1529948111_Z0_T0_C3.tif';

day1.c0 = GeneralAnalysis.loadtiff_1ch(fullfile(day1.datadir,day1.filenamec0));
day1.c0 = squeeze(max(day1.c0,[],3));
day1.c1 = GeneralAnalysis.loadtiff_1ch(fullfile(day1.datadir,day1.filenamec1));
day1.c1 = squeeze(max(day1.c1,[],3));
day1.c2 = GeneralAnalysis.loadtiff_1ch(fullfile(day1.datadir,day1.filenamec2));
day1.c2 = squeeze(max(day1.c2,[],3));
day1.c3 = GeneralAnalysis.loadtiff_1ch(fullfile(day1.datadir,day1.filenamec3));
day1.c3 = squeeze(max(day1.c3,[],3));

day1.sv10 = findshift(day1.c1,day1.c0,'iter');
day1.sv12 = findshift(day1.c1,day1.c2,'iter');
day1.sv13 = findshift(day1.c1,day1.c3,'iter');

%%
day2.datadir = 'F:\Pharm Department Channel Alignment\20180531_alignment_test';
day2.filenamec0 = 'Capture 2_XY1527802513_Z0_T0_C0.tiff';
day2.filenamec1 = 'Capture 2_XY1527802513_Z0_T0_C1.tiff';
day2.filenamec2 = 'Capture 2_XY1527802513_Z0_T0_C2.tiff';
day2.filenamec3 = 'Capture 2_XY1527802513_Z0_T0_C3.tiff';

day2.c0 = GeneralAnalysis.loadtiff_1ch(fullfile(day2.datadir,day2.filenamec0));
% day2.c0 = max(day2.c0,[],3);
day2.c1 = GeneralAnalysis.loadtiff_1ch(fullfile(day2.datadir,day2.filenamec1));
% day2.c1 = max(day2.c1,[],3);
day2.c2 = GeneralAnalysis.loadtiff_1ch(fullfile(day2.datadir,day2.filenamec2));
% day2.c2 = max(day2.c2,[],3);
day2.c3 = GeneralAnalysis.loadtiff_1ch(fullfile(day2.datadir,day2.filenamec3));
% day2.c3 = max(day2.c3,[],3);

day2.sv10 = findshift(day2.c1,day2.c0,'iter');
day2.sv12 = findshift(day2.c1,day2.c2,'iter');
day2.sv13 = findshift(day2.c1,day2.c3,'iter');
%%
day3.datadir = 'F:\Pharm Department Channel Alignment\20180531_alignment_test';
day3.filenamec0 = 'Capture 1_XY1527802393_Z0_T0_C0.tiff';
day3.filenamec1 = 'Capture 1_XY1527802393_Z0_T0_C1.tiff';
day3.filenamec2 = 'Capture 1_XY1527802393_Z0_T0_C2.tiff';
day3.filenamec3 = 'Capture 1_XY1527802393_Z0_T0_C3.tiff';

day3.c0 = GeneralAnalysis.loadtiff_1ch(fullfile(day3.datadir,day3.filenamec0));
% day3.c0 = max(day3.c0,[],3);
day3.c1 = GeneralAnalysis.loadtiff_1ch(fullfile(day3.datadir,day3.filenamec1));
% day3.c1 = max(day3.c1,[],3);
day3.c2 = GeneralAnalysis.loadtiff_1ch(fullfile(day3.datadir,day3.filenamec2));
% day3.c2 = max(day3.c2,[],3);
day3.c3 = GeneralAnalysis.loadtiff_1ch(fullfile(day3.datadir,day3.filenamec3));
% day3.c3 = max(day3.c3,[],3);

day3.sv10 = findshift(day3.c1,day3.c0,'iter');
day3.sv12 = findshift(day3.c1,day3.c2,'iter');
day3.sv13 = findshift(day3.c1,day3.c3,'iter');
%%
shifttest = shift(day1.c3,day2.sv13);
jc = joinchannels('rgb',day1.c1,shifttest)
