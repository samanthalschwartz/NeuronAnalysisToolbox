function returnimage = registerMultiChannelImage(inputimage)
if ~exist('inputimage')
[filename, pathname] = uigetfile('*.*','Select file to register');
fileinfo = fullfile(pathname,filename);
list = {'2','3'};
[indx,~] = listdlg('ListString',list,'SelectionMode','single','PromptString', 'Number of channels in image');
numchannels = list{indx};
inputimage = [];
else
    numchannels = size(inputimage,4);
end
reglist = cellfun(@num2str,num2cell(1:str2num(numchannels)),'uni',false);
[indxR,~] = listdlg('ListString',reglist,'SelectionMode','single','PromptString', 'Select which channel to use for registration');
regchannel = reglist{indxR};

switch numchannels
    case '2'
        if isempty(inputimage)
        inputimage = GeneralAnalysis.loadtiff_2ch(fileinfo);
        end
        ch1im = inputimage(:,:,:,1);
        ch2im = inputimage(:,:,:,2);
        regchnum = str2num(regchannel);
        switch regchnum
            case 1
                [ch1im_shifted,sv_arr] = GeneralAnalysis.timedriftCorrect_parfor(ch1im);
                ch2im_shifted =  GeneralAnalysis.applydriftCorrect(dip_image(ch2im),sv_arr);
            case 2
                [ch2im_shifted,sv_arr] = GeneralAnalysis.timedriftCorrect_parfor(ch2im);
                ch1im_shifted =  GeneralAnalysis.applydriftCorrect(dip_image(ch1im),sv_arr);    
        end
        
    returnimage = cat(4,ch1im_shifted,ch2im_shifted);
        
    case '3'
        if isempty(inputimage)
        inputimage = GeneralAnalysis.loadtiff_3ch(fileinfo);
        end 
        
         ch1im = inputimage(:,:,:,1);
        ch2im = inputimage(:,:,:,2);
        ch3im = inputimage(:,:,:,3);
        regchnum = str2num(regchannel);
        switch regchnum
            case 1
                [ch1im_shifted,sv_arr] = GeneralAnalysis.timedriftCorrect_parfor(ch1im);
                ch2im_shifted =  GeneralAnalysis.applydriftCorrect(dip_image(ch2im),sv_arr);
                ch3im_shifted =  GeneralAnalysis.applydriftCorrect(dip_image(ch3im),sv_arr);
            case 2
                [ch2im_shifted,sv_arr] = GeneralAnalysis.timedriftCorrect_parfor(ch2im);
                ch1im_shifted =  GeneralAnalysis.applydriftCorrect(dip_image(ch1im),sv_arr);
                ch3im_shifted =  GeneralAnalysis.applydriftCorrect(dip_image(ch3im),sv_arr);   
            case 3
                [ch3im_shifted,sv_arr] = GeneralAnalysis.timedriftCorrect_parfor(ch3im);
                ch2im_shifted =  GeneralAnalysis.applydriftCorrect(dip_image(ch2im),sv_arr);
                ch1im_shifted =  GeneralAnalysis.applydriftCorrect(dip_image(ch1im),sv_arr);
        end
        
    returnimage = cat(4,ch1im_shifted,ch2im_shifted,ch3im_shifted);
end

end