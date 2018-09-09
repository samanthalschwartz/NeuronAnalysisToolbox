% first load in all the data (this is slowish)
datadir = 'E:\Matt Becker Data (For Review)\SEPfiles\-TeNT';
% files = uipickfiles('Prompt','Pick Files','FilterSpec',datadir);
files_noT = dir(fullfile(datadir,'*.mat'));
alltraces1_noT = cell(numel(files_noT),1);
alltraces2_noT = cell(numel(files_noT),1);
for ff = 1:numel(files_noT)
sep = load(fullfile(datadir,files_noT(ff).name));
trace1 = sep.obj.oldROI1_sep_sums;
trace2 = sep.obj.oldROI2_sep_sums;
alltraces1_noT{ff} = trace1./mean(trace1(:,1:3),2);
alltraces2_noT{ff} = trace2./mean(trace2(:,1:3),2);
end

datadir = 'E:\Matt Becker Data (For Review)\SEPfiles\+TeNT';
% files = uipickfiles('Prompt','Pick Files','FilterSpec',datadir);
files_T = dir(fullfile(datadir,'*.mat'));
alltraces1_T = cell(numel(files_T),1);
alltraces2_T = cell(numel(files_T),1);
for ff = 1:numel(files_T)
sep = load(fullfile(datadir,files_T(ff).name));
trace1 = sep.obj.oldROI1_sep_sums;
trace2 = sep.obj.oldROI2_sep_sums;
alltraces1_T{ff} = trace1./mean(trace1(:,1:3),2);
alltraces2_T{ff} = trace2./mean(trace2(:,1:3),2);
end


