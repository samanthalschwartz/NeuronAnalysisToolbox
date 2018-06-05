% msdtemporal.m
%
% Mark J. Olah (mjo@cs.unm.edu)
% June 2015
% Compute an MSD vs time lag relatioship for a single track

function [msd, lag, counts]=msdtemporal(track, timeStep)
    % [in] track - an Nx3 matrix where the first three columns are [t x y]
    % [in] (optional) timeStep - the time step in between frames in whatever units the track t's are given in.
    %       e.g. if t's are given as frame index, timeStep=1, but if ts are given in seconds, 
    %       then timeStep is also given in seconds.  If not provided, we estimate this with some extra
    %       computation that could be avoided.
    for ii = 1:numel(track{ii})
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
        msd(n) = mean(sqDisp);
        counts(n) = numel(sqDisp);
    end
    msd = msd(counts>0);
    lag = lag(counts>0);
    counts = counts(counts>0);
end
