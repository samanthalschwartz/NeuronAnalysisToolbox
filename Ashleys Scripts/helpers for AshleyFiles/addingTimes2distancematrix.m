filename = uipickfiles('Prompt','Pick all the AshleyFiles to Calculate Distances (can pick more than 1)');
if ~iscell(filename)
    filesize = 1;
    filename = {filename};
else
    filesize = numel(filename);
end
wb = waitbar(0,'converting to frames and microns....');
for ii = 1:filesize
    clear aa;
    load(fullfile(filename{ii}));
    
    if ischar(aa.imagingparams.postreleaseframe_end)
        aa.imagingparams.postreleaseframe_end = size(aa.surfaceCargo.image,3);
    end
    % now get corresponding time values for frames
    baseframes = aa.imagingparams.baselineframe_start:aa.imagingparams.baselineframe_end;
    basetime_mins = (-flip(baseframes) + 1) * aa.imagingparams.baselineframerate;
    % first frame of post release occurs at time = releasetime
    % all remaining frames at postreleaseframerate intervals + releasetime
    postframes = aa.imagingparams.postreleaseframe_start:aa.imagingparams.postreleaseframe_end;
    firstframe = aa.imagingparams.releasetime;
    remainingframes = (1:(numel(postframes))-1) * aa.imagingparams.postreleaseframerate + firstframe;
    posttime_mins = [firstframe,remainingframes];
    
    timemapping_frames = [baseframes,postframes;];
    timemapping_times = [basetime_mins,posttime_mins];
    aa.distancematrix(:,5) = aa.distancematrix(:,1);
    aa.distancematrix(:,6) = aa.distancematrix(:,2)*aa.pxsize;
    for tt = 1:numel(timemapping_frames)
      idx = aa.distancematrix(:,1) == timemapping_frames(tt);
      aa.distancematrix(idx,5) = timemapping_times(tt);  
    end
    save(fullfile(filename{ii}),'aa');
    waitbar(ii/filesize,wb);
end
close(wb);