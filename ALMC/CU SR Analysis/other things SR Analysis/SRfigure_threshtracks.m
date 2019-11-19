rptdir = 'C:\Users\schwsama\Documents\Data\SR\RPT';
rptfiles = dir(fullfile(rptdir,'*.rpt'));

totalemitters = [];
for rr = 1:numel(rptfiles)
   % load rptfile
   rpt = RPT(fullfile(rptdir,rptfiles(rr).name));
   
   ids = rpt.ResultsTrack.trackLengths>15;
   rpt.ResultsTrack.tracks(ids) = [];
   localizations = [];
   for tt = 1:size(rpt.ResultsTrack.tracks,2)
       localizations = cat(1,localizations,rpt.ResultsTrack.tracks{tt});
   end 
   emitters= localizations(:,[2:11,11,12,12]);
   emitters(:,[1:2,6,7]) = emitters(:,[1:2,6,7])./rpt.data.pixelSize;
   totalemitters = cat(1,totalemitters,emitters);
end

% 
rpt.ResultsFilterEmitters.emitters = totalemitters;
imSizePx = 8192;
physical_roi = rpt.ROI(1:4);
physical_roi(1) = physical_roi(1)-1;
physical_roi(3) = physical_roi(3)-1;
srr = SRRender2D(physical_roi,'single');
points = totalemitters(:,[3,1,2,6,7,13]); %I x y sigmaX sigmaY frameIdx
im = srr.renderGauss(points, imSizePx);
% im(:) = min(im(:), prctile(im(:),99.9999));
rpt.srimage = im;
h = rpt.viewEmitterSuperResGauss();

