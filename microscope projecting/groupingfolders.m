% select the ordering of folders to find the files we want to group
folders1 = {'G:\FromMicroscopeComputer\181220 CamKII test\HomerYFP-CamKIII205K-2b--prelight3',...
    'G:\FromMicroscopeComputer\181220 CamKII test\HomerYFP-CamKIII205K-2b--light',...
    'G:\FromMicroscopeComputer\181220 CamKII test\HomerYFP-CamKIII205K-2b-postlight',...
    'G:\FromMicroscopeComputer\181220 CamKII test\HomerYFP-CamKIII205K-2b-postlight2'};
textstrs = {'Pre-light-3','+488 Light','Post-light-1','Post-light-2'};
numpositions = 7;
channels = {'561'};
% what space to put in between?
combinedfilesfold = 'G:\FromMicroscopeComputer\181220 CamKII test\CombinedMovies';
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
            file = dir(fullfile(folders1{ff},['*' channels{cc} '_s' num2str(pp) '*']));
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
        GeneralAnalysis.LibTiff(newimage,fullfile(combinedfilesfold,[file.name(1:end-4) '_allpositionscombined']));
    end
end