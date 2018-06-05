%Example custom tracks movie generation with changing the frames in the background
%Mark J. Olah (mjo@cs.unm.edu)
% 10/22/15
%
% This is an example of how to use a TrackMovie object to generate a .avi movie.  We also show how to replace the
% background frames with whatever image you want.

%Location of tracked RPT. Change this to your .rpt file location
rpt_file ='/home/mjo/LidkeLab/Data/2015_5_21/RPT/Exp2_2HRG_ROI3.rpt';

%open RPT object
rpt = RPT(rpt_file);

%get the tracks.  This is a cell array of tracks.
tracks = rpt.getTracks();

%Select just one track or several to show
selected_tracks = [1, 7, 15];
my_tracks = tracks(selected_tracks);

%Make TrackMovie that includes only my_tracks
tm = TrackMovie(rpt, my_tracks);
tm.MovieSize = [1080, 720]; % This is the size of the movie in pixels [Width, Height].  
                           % The axes will be made as large as possible to fill the window but maintain the aspect ratio of the data
tm.setTrackColorMethod('Sequence'); % Set coloring method for tracks. This could also be 'Temporal' or 'Speed'

%Configure video options
tm.VideoWriterProfile = 'MPEG-4'; %MPEG-4 works best for smaller files
tm.VideoWriterProfile = 'Uncompressed AVI'; %Archival works best if you will re-encode video later
tm.VideoWriterFrameRate = 25; %Set frame rate to whatever you want.

videoOutFile = 'testMovie'; %Ouput filename

tm.initializeVideoSequence(videoOutFile);
frameSequence = tm.frameBounds(1):tm.frameBounds(2); % This is the default.  It can be any sequence of frame numbers
tm.renderActiveTracks(frameSequence);
tm.finalizeVideoSequence(); % Write file out to disk


%% Now we want to change the images shown behind the movie.  We can do this by replacing
%the tm.frames property with a new image of the same size.

frames = tm.frames;

%This is just a stupid example which makes a dummy background. 
% Replace new_frames with some real movie here
new_frames = zeros(size(frames));
new_frames(:) = 1:numel(frames); 

assert(all(size(new_frames)==size(frames))); %This checks you made the new_frames the right size
tm.ImageGlobalNormalize = true; %set to false to normalize each frame individually set to true to do a global normalization across all frames
tm.setFrames(new_frames); % Must call setFrames to propertly rescale the image

%Also we can change the colormap for the background image
%tm.ImageColorMap = @gray; %default normal grayscale image
tm.ImageColorMap = @hsv; %use any other matlab colormap http://www.mathworks.com/help/matlab/ref/colormap.html?searchHighlight=colormap#inputarg_name

%Now write this out in the new format with the new tracks behind it.
videoOutFile = 'testMovieNewFrames'; %Ouput filename
tm.initializeVideoSequence(videoOutFile);
frameSequence = tm.frameBounds(1):tm.frameBounds(2); % This is the default.  It can be any sequence of frame numbers
tm.renderActiveTracks(frameSequence);
tm.finalizeVideoSequence(); % Write file out to disk

