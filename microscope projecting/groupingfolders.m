% select the ordering of folders to find the files we want to group
folders1 = {'G:\FromMicroscopeComputer\190118 CamKII virus test\300nL_I206KTD_pre',...
    'G:\FromMicroscopeComputer\190118 CamKII virus test\300nL_I206KTD_+488'};
textstrs = {'Pre-light','+488 Light'};
numpositions = 5;
channels = {'561'};
% what space to put in between?
combinedfilesfold = 'G:\FromMicroscopeComputer\190118 CamKII virus test\CombinedMovies';
if ~exist(combinedfilesfold,'dir')
    mkdir(combinedfilesfold)
end

% loop through each posiion
% loop through each channel
% find the appropriate file from each folder and concatenate

for pp = 1:numpositions
    for cc = 1:numel(channels)
        newimage = [];
        for ff = 1:numel(folders1)
            file = dir(fullfile(folders1{ff},['*' channels{cc} '_s' num2str(pp) '_*.tiff']));
             badids = [];
            for tt = 1:numel(file) 
               if file(tt).isdir 
                   badids = [badids,tt];
               end
            end
            file(badids) = [];
            filename = fullfile(folders1{ff},file.name);
            im = loadtiff(filename);
            textframeZ = zeros(size(im,1),size(im,2));
            textframe = insertText(textframeZ,[77 107],textstrs{ff},'FontSize',18,'BoxOpacity',1,'TextColor','white');
            textframe =  textframe(:,:,3)*7000;
            if isempty(newimage)
                newimage = textframe;
            else
                newimage = cat(3,newimage,textframe);
            end
            newimage = cat(3,newimage,im);
        end
        GeneralAnalysis.LibTiff(newimage,fullfile(combinedfilesfold,[file.name(1:end-5) '_allpositionscombined']));
    end
end