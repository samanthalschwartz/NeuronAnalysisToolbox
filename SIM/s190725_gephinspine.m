close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\');
savepath = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\geph_NOTinspines';
wb = waitbar(0);
for ff= 1:numel(filepath)
    tic
    [FILEPATH,NAME,EXT] = fileparts(filepath{ff});
    disp('Loading File...');
    
    load(filepath{ff});
    disp('Done Loading...');
    
    %psd95 obj.ch1
    %geph obj.ch2
    % make cell mask based on psd95ib only
    sumim = obj.ch1.image;
    maxsumim = max(sumim,[],3);
    gmaxsumim = gaussf(maxsumim,10);
    %  threshs = multithresh(single(gmaxsumim),2);
    %  cellmask = gmaxsumim>threshs(1);
    openeninimage_out = opening(gmaxsumim,20,'rectangular');
    threshs = multithresh(single(openeninimage_out),2);
    msk_open = openeninimage_out>threshs(1);
    msk_open = repmat(msk_open,1,1,size(obj.ch1.image,3));
    
    % identify gephyrin within certain distance of psd95
    regionsize = 6;
    psd95_region = bdilation(obj.ch1.mask,regionsize);
    geph_spine = ~psd95_region.*obj.ch2.mask;
    % MatIJ.showImage(cat(4,psd95_region,geph_spine,obj.ch2.mask));
    
    % count total gephyrin
    totgephlb = label(obj.ch2.mask);
    numel(unique(single(totgephlb(:)))); %or obj.ch2.results.numobj
    numgeph = obj.ch2.results.numobj;
    %count total geph spine
    totgephspinelb = label(geph_spine);
    numgephspine = numel(unique(single(totgephspinelb(:)))); %or obj.ch2.results.numobj
    obj.metadata.gephnum.total = numgeph;
    obj.metadata.gephnum.NOTspine = numgephspine;
    obj.metadata.gephnum.ratio = numgephspine/numgeph;
    
    %try making SIM file with new gephyrin channel, then re-run analysis
    obj.ch2.mask = geph_spine;
    % rerun analysis
    obj.make_distancemasks;
    obj.measurements = [{'size'},   {'sum'}, {'Gravity'}];
    obj.measure_AB;
    obj.simulationAbeta(20);
    obj.calculateNumberDensityCOM;
    obj.calculateNumberDensityCOM(0,1);
    disp('Saving File ...');
    obj.save(fullfile(savepath,[NAME '_' num2str(ff)]));
    waitbar(ff/numel(filepath),wb);
    clear obj;
    toc
end
close(wb);
%%
%% ch1 and ch2
close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\071_062819');
savestr = inputdlg('File Save Name');
[FILEPATH,NAME,EXT] = fileparts(filepath{1});
savenameCh1 = fullfile(FILEPATH,[NAME,savestr{1},'_ch1']);
savenameCh2 = fullfile(FILEPATH,[NAME,savestr{1},'_ch2']);

wb = waitbar(0)
for ff= 1:numel(filepath)
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



waitbar(ff/numel(filepath),wb);
end
close(wb)
xlswrite(savenameCh1,ch1matdata);
xlswrite(savenameCh2,ch2matdata);



