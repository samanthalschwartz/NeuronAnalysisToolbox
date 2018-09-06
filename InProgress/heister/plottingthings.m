filepath = 'E:\Matt Becker Data (For Review)\SEPfiles\-TeNT';


files = uipickfiles('Prompt','Pick Files','FilterSpec',filepath);
shafttraces_fixedsize = nan(numel(files),60);
septraces_fixedsize = nan(numel(files),60);
for ff = 1:numel(files)
   sep = load(files{ff});
   currarr_shaft = sep.obj.shafttrace_intensity_fixedsize;
   currarr_sep = sep.obj.spinetrace_intensity_fixedsize;
   shafttraces_fixedsize(ff,1:size(currarr_shaft,1)) = currarr_shaft./currarr_shaft(1);
   septraces_fixedsize(ff,1:size(currarr_sep,1)) = currarr_sep./currarr_sep(1);
end
% septraces_fixedsize(4,:) = [];
% shafttraces_fixedsize(4,:) = [];
save(fullfile(filepath,'traces'),'shafttraces_fixedsize','septraces_fixedsize','files');
%%
filepath = 'E:\Matt Becker Data (For Review)\SEPfiles\+TeNT';
load(fullfile(filepath,'traces'));
titlestr = '+TeNT';
out_shaft = shafttraces_fixedsize(:,all(~isnan(shafttraces_fixedsize)));
out_sep = septraces_fixedsize(:,all(~isnan(septraces_fixedsize)));

fash = figure; plot(out_shaft'); ylim([0.75 1.25]); set(fash,'Position',[2109,104,1289,857]); legend(files,'Location','south'); title(titlestr);
fasp = figure; plot(out_sep'); ylim([0.75 1.25]); set(fasp,'Position',[2109,104,1289,857]); legend(files,'Location','south'); title(titlestr);
saveas(fash,fullfile(filepath,'allShaftTraces'),'png'); 
saveas(fasp,fullfile(filepath,'allSEPTraces'),'png');
close(fash);
close(fasp);

meanshaft = mean(out_shaft,1);
stdshaft = std(out_shaft,1)./sqrt(size(out_shaft,1)); 
fmsh = figure; errorbar(meanshaft,stdshaft,'CapSize',0); ylim([0.75 1.25]); set(fmsh,'Position',[2109,104,1289,857]); title(titlestr);
saveas(fmsh,fullfile(filepath,'meanShaftTraces'),'png');
saveas(fmsh,fullfile(filepath,'meanShaftTraces'),'fig');
close(fmsh)

meansep = mean(out_sep,1);
stdsep = std(out_sep,1)./sqrt(size(out_sep,1));
fmsp = figure; errorbar(meansep,stdsep,'CapSize',0); ylim([0.75 1.25]); set(fmsp,'Position',[2109,104,1289,857]);title(titlestr);
saveas(fmsp,fullfile(filepath,'meanSEPTraces'),'png');
saveas(fmsp,fullfile(filepath,'meanSEPTraces'),'fig');
close(fmsp)