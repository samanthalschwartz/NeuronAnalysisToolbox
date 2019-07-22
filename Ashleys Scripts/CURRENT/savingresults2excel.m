filess = uipickfiles('FilterSpec','Y:\Lab Projects\zapERtrap\Raw Data\LOCAL RELEASE\050918_local_NL1')
load(filess{2})
 [FILEPATH,NAME,EXT] = fileparts(filess{2});

%plot intensity density as a function of time norm to the max somatic intensity
%density
figure;
plot(aa.M.areanormintensity'./aa.M.areanormintensity(1,:)')
M1 = aa.M.areanormintensity'./aa.M.areanormintensity(1,:)'

% % plot the intensity density per time norm to the max intensity
% density for each distance
figure;
plot(aa.M.areanormintensity'./aa.M.areanormintensity(:,71)')
M2 = aa.M.areanormintensity'./aa.M.areanormintensity(:,71)';

xlswrite(fullfile(FILEPATH,[NAME '_results_M1']),M1);
xlswrite(fullfile(FILEPATH,[NAME '_results_M2']),M2)
aa.save(filess{1})
close all
