%% calculated 10 and 50 % values
close all; clear all;
excelFile = 'Y:\Lab Projects\zapERtrap\Raw Data\LOCAL RELEASE\values.xlsx';
filename = 'Y:\Lab Projects\zapERtrap\Raw Data\LOCAL RELEASE\040318_local_soma_NL1\cell4_AshleyFile.mat';
load(filename);
savename = filename(1:end-4);

%%
pxsize = 0.114374*2;
d1 = 5/pxsize;
d2 = 40/pxsize;
d3 = 200/pxsize;
distances = [d1 d2 d3];
if isempty(aa.M)
    M = aa.plotDensityperTime([distances]);
end
figure;
%plot intensity density as a function of time norm to the max somatic intensity
%density
frame_120min = 70;
plot(aa.M.areanormintensity'./aa.M.areanormintensity(1,frame_120min)')
M1 = aa.M.areanormintensity'./aa.M.areanormintensity(1,frame_120min)'

% % plot the intensity density per time norm to the max intensity
% density for each distance
figure;
M3 =  aa.M.areanormintensity';
plot(aa.M.areanormintensity'./aa.M.areanormintensity(:,frame_120min)')
M2 = aa.M.areanormintensity'./aa.M.areanormintensity(:,frame_120min)';
xlswrite(fullfile([savename '_results_areanormintensity']),M3);
xlswrite(fullfile([savename '_results_norm2soma']),M1);
xlswrite(fullfile([savename '_results_norm2each']),M2)
close all
%% 
folder = 'C:\Users\schwsama\Documents';
files = uipickfiles('FilterSpec',folder,'REFilter','MaxEachDistance.csv');
savenamestr = inputdlg('save name');
savename = fullfile(folder,savenamestr);
results.name = {};
results.table = {};
for ff = 1:numel(files)
    csvFile= files{ff};
    [FILEPATH,excelbase,EXT] = fileparts(csvFile);
    excelFile = fullfile(FILEPATH,[excelbase '.xlsx']);
    [data,txt,raw]  = xlsread(csvFile);
    xs = data(:,1);
    perc10 = [];
    perc50 = [];
    cnt = 0;
    for m = 2:4
       cnt = cnt+1;
        ys = data(:,m);
       perc10(cnt) = xs(find(ys>=0.1,1,'first'));
       perc50(cnt) = xs(find(ys>=0.5,1,'first'));
    end
    firstcelldata = raw;
    xlswrite(excelFile,firstcelldata,1);
    secondcelldata = [{'Time to 10% Max'}, {'Time to 50% Max'}; num2cell(perc10(:)),num2cell(perc50(:))];
    try
        secondcelldata = [txt',secondcelldata];
    catch
    end
    xlswrite(excelFile,secondcelldata,2);
    results.name{ff} = csvFile;
   t = table(perc10',perc50');
   t.Properties.RowNames = txt(2:end)';
   t.Properties.VariableNames = {'perc10','perc50'};
   results.table{ff} = t; 
end
perc10_5 = cellfun(@(x) x.perc10(1), results.table);
perc10_40 = cellfun(@(x) x.perc10(2), results.table);
perc10_200 = cellfun(@(x) x.perc10(3), results.table);
perc50_5 = cellfun(@(x) x.perc50(1), results.table);
perc50_40 = cellfun(@(x) x.perc50(2), results.table);
perc50_200 = cellfun(@(x) x.perc50(3), results.table);
allresults = [perc10_5',perc10_40',perc10_200',perc50_5',perc50_40',perc50_200'];
allnames = {'10% 0-5','10% 10-40','10% 10_200','50% 0-4','50% 10-40','50% 40_200'};
resultsmat = [allnames;num2cell(allresults)];
xlswrite(savename,resultsmat);