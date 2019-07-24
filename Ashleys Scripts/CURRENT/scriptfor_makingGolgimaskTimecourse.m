close all; clear all
%filename = 'E:\Zapalog project\090817 COS7 Local Release\movies\Movies For Fig 1d\cell1.tif';
%[out] = loadtiff(filename);
%chGFP = dip_image(out(:,:,:,2));
%chSurface = dip_image(out(:,:,:,3));
%chER = dip_image(out(:,:,:,1));

%Golgi Accumulation
clear all;
close all;
chGFP = GeneralAnalysis.loadtiff_1ch('Z:\Matt\022018_nl1_zapalog\MOVIES\slip2_negcntrl_withzapalog_releaseatend\cell5_g_stabilized.tif');
maxprojGFP = max(chGFP,[],3);
[maskgolgi,threshval_golgi] = GeneralAnalysis.imgThreshold_fixedUserInput(maxprojGFP); 
golgimasktimeseries = repmat(maskgolgi,1,1,size(chGFP,3));
overlay(chGFP,golgimasktimeseries)
% clean up the image:
uiwait(msgbox('Select regions in the mask to remove. Once you are satisfied, close the window.','Clean UP','modal'));
cargomask = GeneralAnalysis.cleanUpMask_manual_square(chGFP,golgimasktimeseries);
output = sum(cargomask.*chGFP,[],[1 2]);
output = single(squeeze(output));
golgi_norm_baseline = sum(output(:,[1 2 3 4 5 6 7 8]))
golgi_norm_baseline = golgi_norm_baseline/8
golgi_intensity = output./golgi_norm_baseline;
figure; plot(golgi_intensity);
figure;uitable('Data',golgi_intensity')


%Surface Accumulation
chSurface = GeneralAnalysis.loadtiff_1ch('E:\Zapalog project\090817 COS7 Local Release\movies\Movies For Fig 1d\cell7_640.tif');
chER = GeneralAnalysis.loadtiff_1ch('E:\Zapalog project\090817 COS7 Local Release\movies\Movies For Fig 1d\cell7_561.tif');
[maskER,threshval_ER] = GeneralAnalysis.imgThreshold_fixedUserInput(chER);
[maskSurface,threshval_Surface] = GeneralAnalysis.imgThreshold_fixedUserInput(chSurface);
ERmask = dip_image(logical(sum(maskER,[],(3))));
ERmask_timeseries = repmat(ERmask,1,1,size(chSurface,3));
overlay(chSurface,ERmask_timeseries)
output = sum(chSurface.*ERmask_timeseries,[],[1 2]);
output = single(squeeze(output));
norm_baseline = sum(output(:,[16 17 18 19 20]))
norm_baseline = norm_baseline/5
%norm_baseline2 = output(:,[1 2 3 4])/norm_baseline
%norm_baseline3 = mean(norm_baseline2)
surface_intensity = output./norm_baseline
figure; plot(surface_intensity);
figure;uitable('Data',surface_intensity')
