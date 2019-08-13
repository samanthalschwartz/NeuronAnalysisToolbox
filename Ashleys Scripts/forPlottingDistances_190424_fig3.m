close all;
clear all;
pxsize = 0.114374*2;
d1 = 5/pxsize;
d2 = 40/pxsize;
d3 = 200/pxsize;
distances = [d1 d2 d3];
files = uipickfiles('FilterSpec','Y:\Lab Projects\zapERtrap\Raw Data\GLOBAL RELEASE\NL1\022018');
for ff = 2:numel(files)
    clear aa;
    try
        load(files{ff})
        aa.calcDensityperTime([distances]);
        
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
for ff = 1:numel(files)
%plot intensity density as a function of time norm to the max somatic intensity
%density

% % plot the intensity density per time norm to the max intensity
% density for each distance
