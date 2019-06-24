% this is the super cleaner
% work on cleaning up the cell fill. to do this:
% first threshold on the sum projection of the cellfill image
% second clean up this single frame image
% then use this image as an additional mask to clean up the cellfill mask

aa.cellFill.mask_img_better; % to better mask cell

masksum = sum(aa.cellFill.mask,[],3); % sum project to bring out all dendrites etc
test = GeneralAnalysis.imgThreshold_fixedUserInput(masksum); % threshold this sum projectsion
sumproj = sum(aa.cellFill.image,[],3); % sum the original image - just for visualization purposes (with cleanup)
clean = GeneralAnalysis.cleanUpMask_manual_square(sumproj,test); % clean up the cell mask!
aa.cellFill.mask = repmat(clean,[1 1 size(aa.cellFill.image,3)]); % replicate the mask for all time points

aa.cleanSurfaceCargoMask_Manual;
aa.cleanSurfaceCargoMaskbyFrame_Manual;
