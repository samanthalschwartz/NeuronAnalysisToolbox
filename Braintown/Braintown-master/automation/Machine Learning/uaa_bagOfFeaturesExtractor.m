function [features, featureMetrics, varargout] = uaa_bagOfFeaturesExtractor(I)

%% Step 1: Preprocess the Image
% The extractor function is applied to each image, I, within the image set
% used to create the bagOfFeatures. Depending on the type of features being
% extracted, the input images may require preprocessing prior to feature
% extraction.

% Convert I to grayscale if required.
[height,width,numChannels] = size(I);
if numChannels > 1
    grayImage = rgb2gray(I);
else
    grayImage = I;
end

%% DSIFT
% tic
% im = single(I);
% binSize = 8 ;
% magnif = 6 ; %this maybe needs to change depending on zoom
% Is = vl_imsmooth(im, sqrt((binSize/magnif)^2 - .25)) ;
% 
% [f, d] = vl_dsift(Is, 'size', binSize,'Norm') ;
% features = double(d');
% featureMetrics = double(f(3,:)');
% toc
% return
% f(3,:) = binSize/magnif ;
% f(4,:) = 0 ;
% [f_, d_] = vl_sift(I, 'frames', f) ;

%% BRISK

% points = detectBRISKFeatures(grayImage);
% features = extractFeatures(grayImage, points);
% features = double(features.Features);
% featureMetrics = points.Metric;
% return

%% SURF
% 
% points = detectSURFFeatures(grayImage);
% features = extractFeatures(grayImage, points);
% featureMetrics = points.Metric;
% return

% %% MSER
% 
% regions = detectMSERFeatures(grayImage);
% features = extractFeatures(grayImage, regions);
% featureMetrics = var(features,[],2);
% 
% return

% %% MSER
% 
% regions = detectMSERFeatures(grayImage);
% features = extractFeatures(grayImage, regions);
% featureMetrics = var(features,[],2);
% 
% return
%%
%% Step 2: Select Point Locations for Feature Extraction
% Here, a regular spaced grid of point locations is created over I. This
% allows for dense SURF feature extraction. 

% Define a regular grid over I.
gridStep = 8; % in pixels
gridX = 1:gridStep:width;
gridY = 1:gridStep:height;

[x,y] = meshgrid(gridX, gridY);

gridLocations = [x(:) y(:)];

%%
% Concatenate multiple SURFPoint objects at different scales to achieve
% multiscale feature extraction.
multiscaleGridPoints = [SURFPoints(gridLocations, 'Scale', 1.6); 
                        SURFPoints(gridLocations, 'Scale', 3.2);
                        SURFPoints(gridLocations, 'Scale', 4.8);
                        SURFPoints(gridLocations, 'Scale', 6.4)];
                    
% Alternatively, you may use a feature detector such as detectSURFFeatures
% or detectMSERFeatures to select point locations. For instance:
%
% multiscaleSURFPoints = detectSURFFeatures(grayImage);
                    
%% Step 3: Extract features
% Finally, extract features from the selected point locations. By default,
% bagOfFeatures extracts upright SURF features. 
features = extractFeatures(grayImage, multiscaleGridPoints,'Upright',true);

%% Step 4: Compute the Feature Metric
% The feature metrics indicate the strength of each feature, where larger
% metric values are given to stronger features. The feature metrics are
% used to remove weak features before bagOfFeatures learns a visual
% vocabulary. You may use any metric that is suitable for your feature
% vectors.
%
% Use the variance of the SURF features as the feature metric.
featureMetrics = var(features,[],2);

% Alternatively, if a feature detector was used for point selection,
% the detection metric can be used. For example:
%
% featureMetrics = multiscaleSURFPoints.Metric;

% Optionally return the feature location information. The feature location
% information is used for image search applications. See the retrieveImages
% and indexImages functions.
if nargout > 2
    % Return feature location information
    varargout{1} = multiscaleGridPoints.Location;
end


