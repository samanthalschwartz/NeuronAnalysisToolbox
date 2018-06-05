% msdtemporal.m
%
% Mark J. Olah (mjo@cs.unm.edu)
% June 2015
% Compute an MSD vs time lag relatioship for a single track

function [outmsd, outlag, outcounts]=msdtemporal_multiTracks(track, timeStep)
    % [in] track - an Nx3 matrix where the first three columns are [t x y]
    % [in] (optional) timeStep - the time step in between frames in whatever units the track t's are given in.
    %       e.g. if t's are given as frame index, timeStep=1, but if ts are given in seconds, 
    %       then timeStep is also given in seconds.  If not provided, we estimate this with some extra
    %       computation that could be avoided.
    % this also includes a fixed maxdisp of 0.5 --- if you want to allow bigger jumps need to change!! 
    all_msd = {};
    all_lag = {};
    all_counts = {};
    maxdisp = 0.5;
    for ii = 1:numel(track)
    if nargin<2 || timeStep<=0
        timeStep= min(diff(track{ii}(:,1)));
    end
   
    ts = round(track{ii}(:,1)./timeStep);
    lag_steps = 1:max(ts)/4;
    nSteps = numel(lag_steps);
    lag = lag_steps*timeStep;
    msd = zeros(nSteps,1);
    counts = zeros(nSteps,1);
    for n=1:numel(lag_steps)
        step = lag_steps(n);
        t = diff(ts(1:step:end,1));
        sqDisp = diff(track{ii}(1:step:end,2)).^2 + diff(track{ii}(1:step:end,3)).^2;
        sqDisp = sqDisp(t==step);
%         sqDisp(sqDisp>maxdisp)=[];
        msd(n) = mean(sqDisp);
        counts(n) = numel(sqDisp);
    end
    msd(counts<=0) = NaN;
    lag(counts<=0) = NaN;
    counts(counts<=0) = NaN;
    all_msd{ii} = msd;
    all_lag{ii} = lag;
    all_counts{ii} = counts;
    end
    temp = cellmap(@(x) numel(x), all_msd); 
    all_msd_arrs = NaN(max(cell2mat(temp)),numel(all_msd));
    all_lag_arrs = NaN(size(all_msd_arrs));
    all_counts_arrs = NaN(size(all_msd_arrs));
    for ii = 1:numel(all_msd)
        all_msd_arrs(1:size(all_msd{ii}),ii) = all_msd{ii};
        all_lag_arrs(1:size(all_msd{ii}),ii) = all_lag{ii};
        all_counts_arrs(1:size(all_msd{ii}),ii) = all_counts{ii};
    end
    
    outcounts = sum(all_counts_arrs,2,'omitnan');
    outcounts = outcounts(outcounts>0);
    outmsd = mean(all_msd_arrs,2,'omitnan');
    outmsd = outmsd(outcounts>0);
    outlag = mean(all_lag_arrs,2,'omitnan');
    outlag = outlag(outcounts>0);
    
end
