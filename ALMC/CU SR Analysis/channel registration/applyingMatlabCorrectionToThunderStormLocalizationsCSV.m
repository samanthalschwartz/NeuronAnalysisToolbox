offset = 25.5
val =[];
outval = [];
theoutval = [];
[outval] = themap([val/1000]);
theoutval = outval;
theoutval(:,1) = theoutval(:,1)+offset;
theoutval = theoutval.*1000;
%%
univals = unique(frames);
v1 = interp1(univals,flip(univals),frames);
%%
filename = 'D:\WDKennedyLabHDDBackup\Projects\Project Cry2Olig-Gephyrin\SR\analysis things\gephyrin images\geph-geph005_SR561.csv';
savename = 'D:\WDKennedyLabHDDBackup\Projects\Project Cry2Olig-Gephyrin\SR\analysis things\gephyrin images\geph-geph005_SR561_chregistered.xlsx';
themap = image.getOptimalMapMicrons;
M = importdata(filename);
xvals = M.data(:,3)/1000 + 25.6;
yvals = M.data(:,4)/1000;
newvals = themap([xvals,yvals]);
M.data(:,3:4) = newvals*1000;

newnames = cellfun(@(x) x(2:end-1),M.colheaders,'UniformOutput',false)
mytab = [M.colheaders;num2cell(M.data)];
xlswrite(savename,mytab);