baselineframe_start = 1;
baselineframe_end = 6;
baselineframerate = 2;
releasetime = 1;
postreleaseframe_start = 7;
postreleaseframe_end = 'end';
postreleaseframerate = 2;
%---
filename = uipickfiles('Prompt','Pick all the AshleyFiles to Calculate Distances (can pick more than 1)');
if ~iscell(filename)
    filesize = 1;
    filename = {filename};
else
    filesize = numel(filename);
end

for ii = 1:filesize
    clear aa;
load(fullfile(filename{ii}));


aa.imagingparams.baselineframe_start =baselineframe_start;
aa.imagingparams.baselineframe_end =baselineframe_end;
aa.imagingparams.baselineframerate =baselineframerate;
aa.imagingparams.releasetime =releasetime;
aa.imagingparams.postreleaseframe_start =postreleaseframe_start;
aa.imagingparams.postreleaseframe_end =postreleaseframe_end;
aa.imagingparams.postreleaseframerate =postreleaseframerate;



save(fullfile(filename{ii}),'aa');

end