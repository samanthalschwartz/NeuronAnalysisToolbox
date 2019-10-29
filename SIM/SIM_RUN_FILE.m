%% setup the SIM files
ids = [2 1 3];
channelorderingstr = {'chABeta','PrP'}; % channel abeta, channel 1, channel2
dirname = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\081';
savedir = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\081';
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec',dirname);
wb = waitbar(0);
for ff= 1:numel(filepath)
    obj  = SIM();
    obj.channelordering = ids;
    obj.channelorderingstr = channelorderingstr;
    obj.filepath = filepath{ff};
    [FILEPATH,NAME,EXT] = fileparts(filepath{ff});
    switch EXT
        case '.nd2'
            obj.loadNDfile();
            if isempty(obj.ch2.rawimage) % this is to handle only 2 channel data
                obj.ch2.rawimage = 0.*obj.ch1.rawimage;
            end
        case {'.tiff','.tif','.TIF','.TIFF'}
            obj.loadtiff; close all;
            obj.XYpxsize = 0.0321;           % returns value in default unit
            obj.XYpxsize_units = 'µm';
            obj.Zpxsize = 0.2000;           % returns value in default unit
            obj.Zpxsize_units = 'µm';
    end
    obj.planeBOTTOM = '3';
    obj.planeTOP = '10';
    obj.setimage();
    disp(['File ' num2str(ff)]);
    obj.save(savedir);
    waitbar(ff/numel(filepath),wb);
end
close(wb)
%% identify top/bottom planes to use
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec',savedir);
wb = waitbar(0);
for ff= 1:numel(filepath)
    load(filepath{ff});
    % -- this is in dipimage values! (one less than image j plane)
    g = dipshow(obj.ch1.rawimage,'log')
    while(ishandle(g))
        try
                w = waitforbuttonpress;
        catch
        end
    end
    opts.WindowStyle='normal';
    prompt = {'Bottom Plane','Top Plane'};
    title = 'select the planes';
    dims = [1 35];
    definput = {'',''};
    answer = inputdlg(prompt,title,dims,definput,opts);
    obj.planeBOTTOM = answer{1};
    obj.planeTOP = answer{2};
    obj.setimage();
    obj.save();
    waitbar(ff/numel(filepath),wb);
end
close(wb)

%% run all the masking, simulation etc
close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\081');
wb = waitbar(0);
for ff= 1:numel(filepath)
    tic
load(filepath{ff});
disp('Masking channels...');
obj.make_cellmask;
obj.make_maskchX;
obj.make_maskch1;
obj.make_maskch2;
disp('Finished Making Masks');
obj.make_distancemasks;
obj.measurements = [{'size'},   {'sum'}, {'Gravity'}];
obj.measure_AB;
obj.simulationAbeta(20);
obj.calculateNumberDensityCOM(0);
obj.calculateNumberDensityCOM(0,1);
obj.save(filepath{ff});
waitbar(ff/numel(filepath),wb);
clear obj;
toc
end
close(wb);

%% ch1 and ch2
close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\081');
savestr = inputdlg('File Save Name');
[FILEPATH,NAME,EXT] = fileparts(filepath{1});
savenameCh1 = fullfile(FILEPATH,[NAME,savestr{1},'_ch1']);
savenameCh2 = fullfile(FILEPATH,[NAME,savestr{1},'_ch2']);

ch2matdata = []; ch1matdata=[];
wb = waitbar(0)
for ff= 3:numel(filepath)
load(filepath{ff});
%ch1 stuff
ch1dat = movmean(obj.ch1.results.numabeta./obj.ch1.results.volume,3);
ch1sim =  movmean(obj.ch1.abetaSIM.results.numabeta./obj.ch1.abetaSIM.results.volume,10);
ch1norm = ch1dat./ch1sim;
ch1bins = obj.ch1.abetaSIM.results.bins;
ch1out_bins = [[obj.filepath ' Bins']; num2cell(ch1bins)'];
ch1out_data = [[obj.filepath ' ch1']; num2cell(ch1norm)];
if ff == 1
    ch1matdata = [ch1out_bins,ch1out_data];
else
    ch1matdata = [ch1matdata,ch1out_data];
end

if numel(obj.channelorderingstr)>2
    %ch2 stuff
    ch2dat = movmean(obj.ch2.results.numabeta./obj.ch2.results.volume,3);
    ch2sim =  movmean(obj.ch2.abetaSIM.results.numabeta./obj.ch2.abetaSIM.results.volume,10);
    ch2norm = ch2dat./ch2sim;
    ch2bins = obj.ch2.abetaSIM.results.bins;
    ch2out_bins = [[obj.filepath ' Bins']; num2cell(ch2bins)'];
    ch2out_data = [[obj.filepath ' ch2']; num2cell(ch2norm)];
    if ff == 1
        ch2matdata = [ch2out_bins,ch2out_data];
    else
        ch2matdata = [ch2matdata,ch2out_data];
    end
end
waitbar(ff/numel(filepath),wb);
end
close(wb)
xlswrite(savenameCh1,ch1matdata);
if ~isempty(ch2matdata)
    xlswrite(savenameCh2,ch2matdata);
end

