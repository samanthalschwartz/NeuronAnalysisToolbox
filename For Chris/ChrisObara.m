filename = 'C:\Users\sammy\Downloads\Sec61-100Hz_010_SRRF.tif';
minbranchlength = 5;
tifim = loadtiff(filename);
tic
[threshim, threshval] = ForChris.imgThreshold_fixedUserInput(tifim);
newim = zeros(size(tifim));
for ii = 1:size(tifim,3)
testim = logical(squeeze(threshim(:,:,ii-1)));
BW3 = bwskel(testim,'MinBranchLength',minbranchlength);
newim(:,:,ii) = BW3;
end
toc
ForChris.viewMaskOverlay(tifim,newim)

%% reshape tracks
% chris tracks are in format
%[time,trackid][x,trackid][y,trackid]
tracksfilename = 'C:\Users\sammy\Downloads\Sec61-011.mat';
load(tracksfilename);
matx = Tracks.matrix;
tracks=  cell(size(matx,2),1);
frrate = 1;
for tr = 1:size(matx,2)
tfr = matx(:,tr,1);
goodid = ~isnan(tfr)&tfr<5000;
fr =  matx(goodid,tr,1)+1;
x = matx(goodid,tr,3)*6.25;
y =  matx(goodid,tr,2)*6.25;
t = matx(goodid,tr,1)*frrate;
inbet = zeros(size(t,1),8);
currtr = [t,x,y,inbet,fr];
currtr(isnan(currtr(:,1)),:) = [];
tracks{tr}= currtr;
end
emptids = cellfun(@isempty,tracks);
tracks(emptids) = [];

frameBounds = [min(cellfun(@(T) min(T(:,end)),tracks)), max(cellfun(@(T) max(T(:,end)),tracks))]; %this is just to ensure that underlying image has same dimensions as trajectories
im = image(:,:,frameBounds(1):frameBounds(2));
tm = TrackMovie(tracks,im);
tm.FrameLag = 10; %-1 is default
tm.TrackColorMap = @jet;
tm.setTrackColorMethod('Sequence');
tm.viewSequence;