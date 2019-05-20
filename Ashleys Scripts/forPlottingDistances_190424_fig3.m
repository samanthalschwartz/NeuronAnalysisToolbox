close all;
clear all;
% filename = 'G:\zapERtrap\Raw Data\GLOBAL RELEASE\GluA1\051718\slip2_7_merge_stitch-AshleyFile.mat';
pxsize = 0.114374*2;
d1 = 5/pxsize;
d2 = 40/pxsize;
d3 = 200/pxsize;
distances = [d1 d2 d3];
files = uipickfiles('FilterSpec','G:\zapERtrap\Raw Data\GLOBAL RELEASE\GluA1\051718');
for ff = 1:numel(files)
    try
        load(files{ff})
        M = aa.plotDensityperTime([distances]);
        
        % show distance 1 area
        % currmask = aa.distmask<=d1;
        % maskd1 = dipshow(currmask)
        %
        % % show distance 2 area
        % currmask = aa.distmask>d1 & aa.distmask<=d2;
        % maskd2 = dipshow(currmask)
        %
        % % show distance 3 area
        % currmask = aa.distmask>d2 & aa.distmask<=d3;
        % maskd3 = dipshow(currmask)
        
        save(files{ff},'aa');
    catch
        continue
    end
    
end
%% can just run this section if you've already calculated everything and just want to plot
% figure;
% %plot density as a function of time
% plot(aa.M.areanormintensity')
for ff = 2:numel(files)
%plot intensity density as a function of time norm to the max somatic intensity
%density
% plot(aa.M.areanormintensity'./aa.M.areanormintensity(1,end)')
load(files{ff});
M1 = aa.M.areanormintensity'./aa.M.areanormintensity(1,end)';

% % plot the intensity density per time norm to the max intensity
% density for each distance
% figure;
% plot(aa.M.areanormintensity'./aa.M.areanormintensity(:,end)')
M2 = aa.M.areanormintensity'./aa.M.areanormintensity(:,end)';
savename = [files{ff}(1:end-4) '_results'];
xlswrite(savename,M1,'M1');
xlswrite(savename,M2,'M2');
clear aa;
end