filepath = 'E:\Matt Becker Data (For Review)\SEPfiles\-TeNT';
titlestr = '-TeNT';
files = uipickfiles('Prompt','Pick Files','FilterSpec',filepath);
shafttraces_fixedsize = nan(numel(files),60);
septraces_fixedsize = nan(numel(files),60);
for ff = [1:numel(files)]
   sep = load(files{ff});
   currarr_shaft = sep.obj.shafttrace_intensity_varsize./sep.obj.shafttrace_size_varsize;
   currarr_sep = sep.obj.spinetrace_intensity_varsize./sep.obj.spinetrace_size_varsize;
   shafttraces_fixedsize(ff,1:size(currarr,1)) = currarr_shaft;
   septraces_fixedsize(ff,1:size(currarr,1)) = currarr_sep;
end

save(fullfile(filepath,'traces_varsize'),'shafttraces_fixedsize','septraces_fixedsize');
out_shaft = shafttraces_fixedsize(:,all(~isnan(shafttraces_fixedsize)));
out_sep = septraces_fixedsize(:,all(~isnan(septraces_fixedsize)));

fash = figure; plot(out_shaft');
ylim([100 900]); 
set(fash,'Position',[2109,104,1289,857]); legend(files,'Location','south')
ylim([100 800]);  title(titlestr);
saveas(fash,fullfile(filepath,'allShaftTraces_varsize'),'png'); 
close(fash);

fasp = figure; plot(out_sep'); 
ylim([100 900]); 
set(fasp,'Position',[2109,104,1289,857]); legend(files,'Location','south')
legend(files,'Location','south');title(titlestr);
saveas(fasp,fullfile(filepath,'allSEPTraces_varsize'),'png');
close(fasp);

meanshaft = mean(out_shaft,1);
stdshaft = std(out_shaft,1)./sqrt(size(out_shaft,1)); 
fmsh = figure; errorbar(meanshaft,stdshaft,'CapSize',0); 
ylim([100 800]);  title(titlestr);
set(fmsh,'Position',[2109,104,1289,857]);
saveas(fmsh,fullfile(filepath,'meanShaftTraces_varsize'),'png');
saveas(fmsh,fullfile(filepath,'meanShaftTraces_varsize'),'fig');
close(fmsh);

meansep = mean(out_sep,1);
stdsep = std(out_sep,1)./sqrt(size(out_sep,1));
fmsp = figure; errorbar(meansep,stdsep,'CapSize',0); 
ylim([100 800]); title(titlestr);
set(fmsp,'Position',[2109,104,1289,857]);
saveas(fmsp,fullfile(filepath,'meanSEPTraces_varsize'),'png');
saveas(fmsp,fullfile(filepath,'meanSEPTraces_varsize'),'fig');
close(fmsp);