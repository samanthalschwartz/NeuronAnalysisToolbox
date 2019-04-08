%%
%----- specify the files to skeletonize and the save directory
filename = 'C:\Users\sammy\Downloads\Sec61-100Hz_010_SRRF.tif';
savename = 'C:\Users\sammy\Downloads\Sec61-100Hz_010_SRRF_skeleton';
minbranchlength = 5;
%---- load the file
tifim = loadtiff(filename);
tifim = gaussf(tifim);
% tifim = tifim(:,:,0:20);% -- use this to just analyze a few frames.
%----- select a background region
[threshim, threshval] = ForChris.imgThreshold_fixedUserInput(tifim);
%----- now make the skeleton
tic
skeleton_image = zeros(size(tifim));
wb = waitbar(0,'Calculating Skeleton...');
for ii = 1:size(tifim,3)
testim = logical(squeeze(threshim(:,:,ii-1)));
BW3 = bwskel(testim,'MinBranchLength',minbranchlength);
skeleton_image(:,:,ii) = BW3;
waitbar(ii/size(tifim,3),wb);
end
toc
close(wb)
%----- save results
save(savename,'skeleton_image','threshval');
LibTiff(skeleton_image,savename);
%% view results with dipimage
%--- use 'n' and 'p' keys to go next and previous frame
%--- use 'f' and 'b' keys to toggle on and off the mask
%--- this is pretty slow moving between frames but allows you to see the
            %mask
frames2view =  1:5;
ForChris.viewMaskOverlay(tifim(:,:,frames2view),skeleton_image(:,:,frames2view))
%% reshape tracks
% chris tracks are in format
%[time,trackid][x,trackid][y,trackid]
tracksfilename = 'C:\Users\sammy\Downloads\Sec61-011.mat';
rawimage = loadtiff('C:\Users\sammy\Downloads\Sec61-100Hz_010.nd2_Raw.tif');
load(tracksfilename);
matx = Tracks.matrix;
tracks=  cell(size(matx,2),1);
framerate = 1;
for tr = 1:size(matx,2)
tfr = matx(:,tr,1);
goodid = ~isnan(tfr)&tfr<5000;
fr =  matx(goodid,tr,1)+1;
x = matx(goodid,tr,3)*6.25;
y =  matx(goodid,tr,2)*6.25;
t = matx(goodid,tr,1)*framerate;
inbet = zeros(size(t,1),8);
currtr = [t,x,y,inbet,fr];
currtr(isnan(currtr(:,1)),:) = [];
tracks{tr}= currtr;
end
emptids = cellfun(@isempty,tracks);
tracks(emptids) = [];
frameBounds = [min(cellfun(@(T) min(T(:,end)),tracks)), max(cellfun(@(T) max(T(:,end)),tracks))]; %this is just to ensure that underlying image has same dimensions as trajectories
im = rawimage(:,:,frameBounds(1):frameBounds(2));

tm = TrackMovie(tracks,im);
tm.FrameLag = -1; %-1 is default. Change this to change how many frames the tracks stay on.
tm.TrackColorMap = @jet;
tm.setTrackColorMethod('Sequence');
tm.viewSequence;