% sample script to make movie of tracks overlayed on an image
% make sure that your mathlab path has access to the SMA_Toolbox and you
% have run 'startupSMA' on the command line

%% ----- if you have an .rpt file
savename = fullfile(savedir,rpt.saveFileBaseName); %change this to change the save name

rpt = RPT(fullfile(dirname,files(ff).name));
tm = TrackMovie(rpt);
tm.setTrackColorMethod('Sequence'); %there are a number of different coloring schemes to use. see TrackMovie class for more info.
tm.FrameLag =20; % number of frames after a given localization for a trajectory to remain visible
tm.MovieSize = [1080, 720]; % change this to change frame size for movie. this default works nicely to see tracks.
tm.VideoWriterProfile = 'MPEG-4'; %MPEG-4 works best for smaller files
tm.VideoWriterFrameRate = 20; %Set frame rate to whatever you want.
tm.TrackLineWidth = 2; % adjust this to make track lines more thick or thin.
videoOutFile = savename;
tm.initializeVideoSequence(videoOutFile);
frameSequence = tm.frameBounds(1):tm.frameBounds(2); % This is the default.  It can be any sequence of frame numbers
tm.renderActiveTracks(frameSequence);
tm.finalizeVideoSequence();


%% ---- to make from tracks and image file
% -- track: must be in the following format (matching the RPT class
%       format) % Table/Matrix columns format for Track results
        % This defines the "RPT" track format although we can covert between
        % (1) cell-array of matricies. One per track
        % (2) single table with a track index
        % (3) structure array
%         NTrackColumns = 12;
%         TrackColumnNames = {'t', 'x', 'y', 'I', 'bg', 'sigma', 'SE_x', 'SE_y', 'SE_I', 'SE_bg', 'SE_sigma','frame'};
%         TrackColumnUnits = {'s', 'um', 'um', 'photons', 'photons/px', 'um', 'um', 'um', 'photons', 'photons/px', 'um','index'};        
%         TrackColumnDescriptions = {'Time', 'x-Position estimate', 'y-Position estimate', 'Intensity estimate',...
%                                    'Mean background intensity per pixel estimate','Apparent gaussian sigma estimate',...
%                                    'Standard error of x-position estimate','Standard error of y-position estimate',...
%                                    'Standard error of intensity estimate', 'Standard error of background intensity estimate',...
%                                    'Standard error of apparent gaussian sigma estimate',...
%                                    'Frame index'};

% --- im_raw:  the image you want to display the trajectories over (need to
%               load in tiff file
savename = ''; % full path to filename you want to use for saving the movie
im_raw = [];
tracks = [];

frameBounds = [min(cellfun(@(T) min(T(:,end)),tracks)), max(cellfun(@(T) max(T(:,end)),tracks))]; %this is just to ensure that underlying image has same dimensions as trajectories
im = im_raw(:,:,frameBounds(1):frameBounds(2));
tm = TrackMovie(tracks,im);
tm.setTrackColorMethod('Sequence');
tm.FrameLag =20;
tm.MovieSize = [1080, 720];
tm.VideoWriterProfile = 'MPEG-4'; %MPEG-4 works best for smaller files
tm.VideoWriterFrameRate = 20; %Set frame rate to whatever you want.
videoOutFile = savename;
tm.initializeVideoSequence(videoOutFile);
frameSequence = tm.frameBounds(1):tm.frameBounds(2); % This is the default.  It can be any sequence of frame numbers
tm.renderActiveTracks(frameSequence);
tm.finalizeVideoSequence();


