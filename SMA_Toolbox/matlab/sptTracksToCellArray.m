% Mark J. Olah (mjo@cs.unm.edu)
% 03 -2015

function tracks = sptTracksToCellArray(spt, frameT, pixelSize)
    % Given an sptobject or file return the tracks as an array
    % [IN]
    %  spt: spt object or fullpath to .spt file
    %  frameT: frame time in seconds (default = 1s)
    %  pixelSize: pixel size in microns [X Y] (default = [1 1])
    % [OUT]
    %  tracks: cell array of length = num_tracks.
    %          Each entry of the cell array is a matrix describing one track.
    %          This is the RPT format.  Rows = localizations.
    %          Cols: [t(s) x(um) y(um) I bg sigma varX(um^2) varY(um^2) varI varBG varSigma FrameIdx]
    %          FrameIdx is 1-based.  Times (t) start at 0.  Upper left corner of image is (x,y)=(0,0);
    if ischar(spt)
        spt = SPT(spt);
    end
    if nargin<2
        frameT=1;
    end
    if nargin<3
        pixelSize=[1 1];
    else
        pixelSize = pixelSize(:)' .* [1 1];
    end
    origin = spt.ParamsGeneral.Roi([1,3])+[0 1];
    lengths=arrayfun(@(t) numel(t.X), spt.Tracks); %Get length of each track
    [~,sidx]=sort(lengths,'descend'); %get indexs of sorted tracks
    tracks = cell(1,length(spt.Tracks));
    for i = sidx
        t = spt.Tracks(i);
        len = length(t.Frame);
        if isfield(t,'Bg')
            bg = t.Bg';
        else
            bg = zeros(len,1);
        end
        tracks{i} = [frameT*(t.Frame'-1), ...        %t
                     pixelSize(1)*(t.X'+origin(1)), pixelSize(2)*(t.Y'+origin(2)),... %x y
                     t.Photons', bg, ones(len,1),...  %I bg sigma         
                     pixelSize(1)^2*t.std_x', pixelSize(2)^2*t.std_y',...%var_x var_y
                     t.std_Photons'.^2, zeros(len,1), zeros(len,1), t.Frame']; %var_i var_bg var_sigma frameIdx
    end 
end
