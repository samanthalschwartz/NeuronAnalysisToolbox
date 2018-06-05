% Mark J. Olah (mjo@cs.unm.edu)
% June 2015
% Compute an empirical cumulative distribution of squared displacments for a
% cell-array of tracks

function [F, xs] = squaredDisplacementCDF(tracks, stepsize, thinning)
    % [in] tracks - single track or a cell array of tracks.  
    %       Each track is Nx3 with the columns being [t x y] (extra columns
    %       beyond 3 are ignored)
    % [in] stepsize -(optional) step size in units of whatever 't' column is recorded in. 
    % [in] thinning - (optional) the step size at which to sub-sample the squared displacements. [default =2]
    % [out] F - vector of empirical cdf values
    % [out] xs - vector of points the cdf is reported at
    if ~iscell(tracks); 
        tracks={tracks}; 
    end
    if nargin<3
        thinning=2;
    end
    if nargin<2 || stepsize<=0 || isempty(stepsize)
        stepsize= min(cellmatfun(@(t) diff(t(:,1)'), tracks));
    end
    disps = cellmatfun(@(ts) squaredDisplacements(ts,stepsize,thinning), tracks);
    if isempty(disps)
        F = [];
        xs = [];
        return;
    end
    [F, xs] = ecdf(disps);  %empircal CDF
end

function ds = squaredDisplacements(track, stepsize, thinning)
    % compute the squared displacements
    if size(track,2)<3
        ds = [];
        return;
    end
    ds = diff(track(:,2)).^2 + diff(track(:,3)).^2;
    keepers = abs(diff(track(:,1)) - stepsize) < (stepsize/2); %keep only jumps where result is fit in each frame
    ds = ds(keepers);
    ds = ds(1:thinning:end)';
end
