clearvars %Clear workspace
close all %Close all figures

%%
path = uigetdir; %Select data folder
content = dir(fullfile(path,'*.tif')); %Open images that are in tif format & assign the full path to the folder
names = {content.name}'; %Assign folder to variable "names"

%testN = names{1};
%outS = cellfun(@(x) strsplit(x,{'_','-','.'}), names, 'UniformOutput', false) ;
%outS2 = cellfun(@(x) x{12}, outS , 'UniformOutput', false);
%channels1 = outS2(ismember(outS2, 'w0000'));

channels = cellfun(@isempty,regexp(names,'w0000'))+1; % Cellfun applies function to each cell in cell array; use regexp to search for files in folder that contain "w0000" (this is the red channel)
%at symbol is anonymous function; telling matlab that following this at
%symbol is a function call
%cellfun is an important function because names is a cell array, not a vector
%Running a function using a cell array
slice_num = (regexp(names, '(?<=_t)\d+(?=_)','match')); %Search for the timepoint in each tif file
slice_num = cellfun(@(x) str2num(x{:}),slice_num)+1; %Convert string to number

%slice_num2 = cellfun(@(x) str2double(x{11}(2:end)), outS , 'UniformOutput', true);

%%
meta = imfinfo(fullfile(path,names{1})); %Get metadata
num_slices = numel(names)/2; %Number of timepoints is 24 (half of the number of elements which is 48)
% IMG = zeros(meta(1).Height, meta(1).Width,2,num_slices,'uint16');
IMG = zeros(meta(1).Height, meta(1).Width,2,num_slices);

max_max_pix = [0 0]; %Take the maximum of the max pixel intensity for constrast adjustment
for n_names = 1:numel(names)
    
    stack_name = fullfile(path,names{n_names});
    meta2 = imfinfo(stack_name); %by giving it the full path, it doesn't matter where I'm at on my computer
    % iminfo generates struct
    channel = channels(n_names);
    stack = zeros(meta(1).Height, meta(1).Width, numel(meta2),'uint16');
    %16 bit depth
    for m=1:numel(meta2) %fill stack with layer (my inner for loop)
        img = imread(fullfile(path,names{n_names}),m);
        max_pix = max(img(:));
        
        if max_pix > max_max_pix(channel) %greatest wins out
        max_max_pix(channel) = max_pix;
        end       
        stack(:,:,m) = img; %filling stack, 3rd dimension is 1st layer
    end
    
    %projection = max(stack,[],3);
    projection = mean(stack,3); % now make mean projection, convert 3D to 2D image; mean is only looking at m (the layer)
    
    IMG(:,:,channel, slice_num(n_names))=projection; %create movie 4D
end

IMG(:,:,1,:) = IMG(:,:,1,:) / max_max_pix(1); %separating channels and dividing by max pix
IMG(:,:,2,:) = IMG(:,:,2,:) / max_max_pix(2);

%% READ TEXT FILE

content = dir(fullfile(path,'*.txt')); %Open file in txt format (conveniently, there's only 1 file per dataset)

file_path = fullfile(path, content.name);

fid = fopen(file_path);
tline = fgetl(fid);
 
idx = 1 %initialize iterator
frappa = zeros(11,3); %Preallocation with a matrix larger than necessary to account for varying number of FRAPPA points
while ischar(tline) %asks if string because the end is numeric (-1); how while loop ends at end of the file
    
    %search for _FRAPPA lines
    re = regexp(tline,'(?<=_FRAPPA\tPoint\t NumberOfPoints\( \d\) : \( )\d+, \d+','match');
    
    if ~isempty(re) %if a match; if line contains FRAPPA points, re will not be empty
        lb = regexp(tline,'(?<=Label\( )\d','match'); %get label number
        
        frappa(idx,:) = [str2num(re{:}) str2num(lb{:})]; %Convert string to number
        idx = idx+1; %Plus one increments for iterator
    end
    
    tline = fgetl(fid); %read each line in .txt file %moves to next line when fgetl is called
end
fclose(fid); %always have to close a file after opening it

% Cleanup FRAPPA
max_label = max(frappa(:,3)); %only keep number of FRAPPA points corresponding to max label number
frappa(max_label+1:end,:) = []; % Remove all other FRAPPA points (they're just repeats of the original points)
frappa(:,3)=[]; %Remove label number from FRAPPA points (we just want 1st and 2nd columns for x and y coordinates)

frappa = circshift(frappa,1,2); %switch x and y with circshift

%% create mask
% MASK = false(size(IMG));

img = squeeze(IMG(:,:,1,6)); %Make 3D binary stack which only cares about 1 and 0 to create our boundaries
%squeeze out channel %6th projection
mask = false(size(img)); %initialize mask, false will also make zeros like zeros but its logical zeros

for n=1:size(frappa,1) %For loop to plot FRAPPA points on image
    mask(frappa(n,1), frappa(n,2)) = true; %now we're making them equal 1
end

dist_mask = bwdist(mask)<4; %Creating ellipse ROI with radius of 4 pixels; anything <4=1 to make ROI

MASK = repmat(dist_mask,1,1,2,num_slices); %repmat replicates matrices so ROIs are on ALL 24 projection images and the other channel
% first input argument is WHAT I want to replicate, the rest are just
% dimensions %1,1 because we don't want to change x,y dimensions

 figure; imshowpair(imadjust(img),dist_mask) %imshowpair is quicker function than subplot so we can view BOTH contrast adjusted image and ROIs

%%

RED = squeeze(IMG(:,:,1,:)); %squeeze removes elements with dimensions of 1 (squeeze out channel)
GREEN = squeeze(IMG(:,:,2,:));

%% Get average background signal
% Get the average of the 1000 lowest pixel intensities
% h = imhist(imadjust(img)) %if we want to visualize histogram of pixel
% intensity
 A = sort(img,'ascend');
 y = A(1:1000);
 meanBackgroundsignal = mean(y)
 
RED = RED - meanBackgroundsignal;

%%
% Generate table with Background Subtracted Mean Intensity calculations %
%mean_IMG = mean(RED,3); %calculate mean intensity for red channel ROIs

mean_intensity = cell(size(RED,3),1);
for mi=1:size(RED,3);
    mean_intensity{mi,1} = regionprops('table',dist_mask, RED(:,:,mi),'MeanIntensity')
end

%rp = regionprops('table',dist_mask, mean_IMG,'MeanIntensity') %create table of mean intensity values
% Ex for table indexing: mean_intensity{1,1}.MeanIntensity(1)

%% Convert tables to cell arrays & combine
Timepoint_1 = table2cell(mean_intensity{1,1});
Timepoint_2 = table2cell(mean_intensity{2,1});
Timepoint_3 = table2cell(mean_intensity{3,1});
Timepoint_4 = table2cell(mean_intensity{4,1});
Timepoint_5 = table2cell(mean_intensity{5,1});
Timepoint_6 = table2cell(mean_intensity{6,1});
Timepoint_7 = table2cell(mean_intensity{7,1});
Timepoint_8 = table2cell(mean_intensity{8,1});
Timepoint_9 = table2cell(mean_intensity{9,1});
Timepoint_10 = table2cell(mean_intensity{10,1});
Timepoint_11 = table2cell(mean_intensity{11,1});
Timepoint_12 = table2cell(mean_intensity{12,1});
Timepoint_13 = table2cell(mean_intensity{13,1});
Timepoint_14 = table2cell(mean_intensity{14,1});
Timepoint_15 = table2cell(mean_intensity{15,1});
Timepoint_16 = table2cell(mean_intensity{16,1});
Timepoint_17 = table2cell(mean_intensity{17,1});
Timepoint_18 = table2cell(mean_intensity{18,1});
Timepoint_19 = table2cell(mean_intensity{19,1});
Timepoint_20 = table2cell(mean_intensity{20,1});
Timepoint_21 = table2cell(mean_intensity{21,1});
Timepoint_22 = table2cell(mean_intensity{22,1});
Timepoint_23 = table2cell(mean_intensity{23,1});
Timepoint_24 = table2cell(mean_intensity{24,1});
Timepoint_Total = [Timepoint_1; Timepoint_2; Timepoint_3; Timepoint_4; Timepoint_5; Timepoint_6; Timepoint_7; Timepoint_8; Timepoint_9; Timepoint_10; Timepoint_11; Timepoint_12; Timepoint_13; Timepoint_14; Timepoint_15; Timepoint_16; Timepoint_17; Timepoint_18; Timepoint_19; Timepoint_20; Timepoint_21; Timepoint_22; Timepoint_23; Timepoint_24]

%% Normalize recovery to baseline
% Normalize recovery to baseline

% Calculate mean of baseline
BaselineMean = mean(cellfun(@(x) x(1), Timepoint_Total(1:12,:)))

% Convert cell to double
Timepoint_Total = cell2mat(Timepoint_Total)

%Subtract baseline mean from every point
SubtractedBaseline = Timepoint_Total / BaselineMean
Sub_Timepoint_1 = mean(SubtractedBaseline(1:4,:));
Sub_Timepoint_2 = mean(SubtractedBaseline(5:8,:));
Sub_Timepoint_3 = mean(SubtractedBaseline(9:12,:));
Sub_Timepoint_4 = mean(SubtractedBaseline(13:16,:));
Sub_Timepoint_5 = mean(SubtractedBaseline(17:20,:));
Sub_Timepoint_6 = mean(SubtractedBaseline(21:24,:));
Sub_Timepoint_7 = mean(SubtractedBaseline(25:28,:));
Sub_Timepoint_8 = mean(SubtractedBaseline(29:32,:));
Sub_Timepoint_9 = mean(SubtractedBaseline(33:36,:));
Sub_Timepoint_10 = mean(SubtractedBaseline(37:40,:));
Sub_Timepoint_11 = mean(SubtractedBaseline(41:44,:));
Sub_Timepoint_12 = mean(SubtractedBaseline(45:48,:));
Sub_Timepoint_13 = mean(SubtractedBaseline(49:52,:));
Sub_Timepoint_14 = mean(SubtractedBaseline(53:56,:));
Sub_Timepoint_15 = mean(SubtractedBaseline(57:60,:));
Sub_Timepoint_16 = mean(SubtractedBaseline(61:64,:));
Sub_Timepoint_17 = mean(SubtractedBaseline(65:68,:));
Sub_Timepoint_18 = mean(SubtractedBaseline(69:72,:));
Sub_Timepoint_19 = mean(SubtractedBaseline(73:76,:));
Sub_Timepoint_20 = mean(SubtractedBaseline(77:80,:));
Sub_Timepoint_21 = mean(SubtractedBaseline(81:84,:));
Sub_Timepoint_22 = mean(SubtractedBaseline(85:88,:));
Sub_Timepoint_23 = mean(SubtractedBaseline(89:92,:));
Sub_Timepoint_24 = mean(SubtractedBaseline(93:96,:));

% Normalize each point to FRAP
Norm_Timepoint_1 = ((Sub_Timepoint_1-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_2 = ((Sub_Timepoint_2-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_3 = ((Sub_Timepoint_3-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_4 = ((Sub_Timepoint_4-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_5 = ((Sub_Timepoint_5-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_6 = ((Sub_Timepoint_6-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_7 = ((Sub_Timepoint_7-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_8 = ((Sub_Timepoint_8-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_9 = ((Sub_Timepoint_9-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_10 = ((Sub_Timepoint_10-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_11 = ((Sub_Timepoint_11-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_12 = ((Sub_Timepoint_12-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_13 = ((Sub_Timepoint_13-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_14 = ((Sub_Timepoint_14-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_15 = ((Sub_Timepoint_15-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_16 = ((Sub_Timepoint_16-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_17 = ((Sub_Timepoint_17-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_18 = ((Sub_Timepoint_18-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_19 = ((Sub_Timepoint_19-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_20 = ((Sub_Timepoint_20-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_21 = ((Sub_Timepoint_21-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_22 = ((Sub_Timepoint_22-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_23 = ((Sub_Timepoint_23-Sub_Timepoint_4)/(1-Sub_Timepoint_4));
Norm_Timepoint_24 = ((Sub_Timepoint_24-Sub_Timepoint_4)/(1-Sub_Timepoint_4));

Norm_Timepoint_Total = [Norm_Timepoint_1; Norm_Timepoint_2; Norm_Timepoint_3; Norm_Timepoint_4; Norm_Timepoint_5; Norm_Timepoint_6; Norm_Timepoint_7; Norm_Timepoint_8; Norm_Timepoint_9; Norm_Timepoint_10; Norm_Timepoint_11; Norm_Timepoint_12; Norm_Timepoint_13; Norm_Timepoint_14; Norm_Timepoint_15; Norm_Timepoint_16; Norm_Timepoint_17; Norm_Timepoint_18; Norm_Timepoint_19; Norm_Timepoint_20; Norm_Timepoint_21; Norm_Timepoint_22; Norm_Timepoint_23; Norm_Timepoint_24]

%% Make FRAP curve showing mean recovery
figure(2);
x = 0:23,1;
y = Norm_Timepoint_Total;
f = plot(x,y,'m','MarkerSize',20,'LineWidth', 2);
xlabel('Time (min)')
ylabel('Fluorescence Recovery')
xlim([1 24])
ylim([0 1.2])

%%
%projection = squeeze(max(mat2gray(IMG),[],3));
%mask = false(size(projection));

%%
%ha = imshow(projection);
%colormap(parula)
% Re-run this code after drawing ROI to continue creating ROIs %
%he = imellipse(gca);
%position = wait(he);
%mask(createMask(he)) = true;
%imshowpair(projection,mask)


%%
% Calculate the Average Mean Intensity for all ROIs
% BS = repmat(mask,1,1,num_slices);


% mean(IMG(BS))