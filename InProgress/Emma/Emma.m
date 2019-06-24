%% select the file to analayze
% topdir = '/Volumes/Boxer/example data/Slice 2 RS';
topdir = 'C:\Users\sammy\Dropbox\Sam Kennedy Lab\Emma';
output = uigetfile(fullfile(topdir,'*.*'));
[~,NAME,EXT] = fileparts(output)
basename = NAME(1:end-2);
% savepath = '/Users/emmaboxer/Documents/PVsynapsescellfillexcelsheets';
savepath = fullfile(topdir,'test');
ext = EXT;
ch_cellfillstr = 'C2';
ch_tag1 = 'C1';
ch_tag2 = 'C0';
% and load
clear cellfill tag1 tag2
ch_cellfillpath = fullfile(topdir,[basename ch_cellfillstr ext]);
ch_tag1path = fullfile(topdir,[basename ch_tag1 ext]);
ch_tag2path = fullfile(topdir,[basename ch_tag2 ext]);
uiopen(ch_cellfillpath); cellfill.image = image;
uiopen(ch_tag1path); tag1.image = image;
uiopen(ch_tag2path); tag2.image = image; close all;
joinchannels('rgb',stretch(tag1.image),stretch(tag2.image),stretch(cellfill.image))
savename = fullfile(savepath,[strrep(basename,' ','_') 'results']);

%% Mask the Cell Fill
numslices= size(cellfill.image,3);
T = adaptthresh(cellfill.image,0.1,'ForegroundPolarity','dark');
cellfill_b = single(cellfill.image) - 4*single(T);
cellfill_l = GeneralAnalysis.imgLaplaceCutoff(cellfill_b); %laplace cuttoff filter
cellfill_g = gaussf(cellfill_b); % gaussfilter
[~ ,threshval_l,~]= GeneralAnalysis.imgThreshold_fixedUserInput(cellfill_l(:,:,ceil(numslices/2))); %user select
cellfill_lm = cellfill_l>threshval_l;
[~ ,threshval_g,~]= GeneralAnalysis.imgThreshold_fixedUserInput(cellfill_g(:,:,ceil(numslices)-4)); %user select - change the '1' to select which frame to choose bg from
cellfill_gm = cellfill_g>threshval_g;
cellfill.mask = cellfill_lm | cellfill_gm;
GeneralAnalysis.viewMaskOverlay(cellfill.image,cellfill.mask)
%% clean up the cell fill if you want
cellfill.mask = GeneralAnalysis.cleanUpMask_manual_square(dip_image(cellfill.image),cellfill.mask,100);
%% select 'somatic' region
h = dipshow(max(cellfill.image,[],3),'lin')
diptruesize(h,50);
[roi, v] = diproi(h);
s_mask = repmat(roi,[1 1 size(cellfill.image,3)]);
soma_vertices = v;
close(h);
soma.mask = s_mask.*cellfill.mask;
soma.area = sum(soma.mask(:));
joinchannels('rgb',cellfill.mask,soma.mask)
%% select ch1 threshold
tag1_l = GeneralAnalysis.imgLaplaceCutoff(tag1.image); %laplace cuttoff filter
[~ ,threshval_t1,C]= GeneralAnalysis.imgThreshold_fixedUserInput(tag1_l(:,:,ceil(numslices/2))); %user select
tag1.mask = tag1_l>threshval_t1;
GeneralAnalysis.viewMaskOverlay(tag1.image,tag1.mask)
dipshow(tag1.image)
%% clean up ch1 (if you want)
tag1.mask = GeneralAnalysis.cleanUpMask_manual_square(dip_image(tag1.image),tag1.mask,100);
%% select ch2 threshold
tag2_l = GeneralAnalysis.imgLaplaceCutoff(tag2.image); %laplace cuttoff filter
[~ ,threshva2_t1,C]= GeneralAnalysis.imgThreshold_fixedUserInput(tag2_l(:,:,ceil(numslices/2))); %user select
tag2.mask= tag2_l>threshva2_t1;
GeneralAnalysis.viewMaskOverlay(tag2.image,tag2.mask)
%% clean up ch2 (if you want)
tag2.mask = GeneralAnalysis.cleanUpMask_manual_square(dip_image(tag2.image),tag1.mask,100);
%% label tag1 puncta
tag1_lb = label(tag1.mask,1,4);
dipshow(tag1_lb,'labels');
tag1_soma_lb = label(tag1.mask.*soma.mask,1,4); % this is labeled mask inside selected somatic region
dipshow(tag1_soma_lb,'labels');
tag1_cell_lb = label(tag1.mask.*cellfill.mask,1,4); % this is labeled mask inside selected somatic region
dipshow(tag1_cell_lb,'labels');
% label tag2 puncta
tag2_lb = label(tag2.mask,1,4);
dipshow(tag2_lb,'labels'); 
tag2_soma_lb = label(tag2.mask.*soma.mask,1,4); % this is labeled mask inside selected somatic region
dipshow(tag2_soma_lb,'labels'); 
tag2_cell_lb = label(tag2.mask.*cellfill.mask,1,4); % this is labeled mask inside selected somatic region
dipshow(tag2_cell_lb,'labels'); 

%% If you are OK with all the masks then run measurement calculations

ids_withtag2_all = single(unique(tag1_lb.*tag2.mask));
ids_withtag2_all(ids_withtag2_all==0) = [];
ids_withtag2_soma = single(unique(tag1_soma_lb.*tag2.mask));
ids_withtag2_soma(ids_withtag2_soma==0) = [];
ids_withtag2_cell = single(unique(tag1_cell_lb.*tag2.mask));
ids_withtag2_cell(ids_withtag2_cell==0) = [];
tag1.colocWtag2.totalmask = ismember(single(tag1_lb),ids_withtag2_all);
tag1.colocWtag2.cellmask = ismember(single(tag1_cell_lb),ids_withtag2_cell);
tag1.colocWtag2.somamask = ismember(single(tag1_soma_lb),ids_withtag2_soma);

%measure tag1
tag1.total_msr = measure(tag1_lb,tag1.image,{'size','sum'});
tag1.total_numpuncta = size(tag1.total_msr,1);
tag1.soma_msr = measure(tag1_soma_lb,tag1.image,{'size','sum'});
tag1.soma_numpuncta = size(tag1.soma_msr,1);
tag1.cell_msr = measure(tag1_cell_lb,tag1.image,{'size','sum'});
tag1.cell_numpuncta = size(tag1.cell_msr,1);

%measure tag1 colocalized with tag2
tag1.colocWtag2.total_msr = tag1.total_msr(ids_withtag2_all')
tag1.colocWtag2.total_numpuncta = size(ids_withtag2_all',1);
tag1.colocWtag2.soma_msr = tag1.soma_msr(ids_withtag2_soma');
tag1.colocWtag2.soma_numpuncta = size(ids_withtag2_soma',1);
tag1.colocWtag2.cell_msr = tag1.cell_msr(ids_withtag2_cell');
tag1.colocWtag2.cell_numpuncta = size(ids_withtag2_cell',1);

%measure tag2
tag2.total_msr = measure(tag2_lb,tag2.image,{'size','sum'});
tag2.total_numpuncta = size(tag2.total_msr,1);
tag2.soma_msr = measure(tag2_soma_lb,tag2.image,{'size','sum'});
tag2.soma_numpuncta = size(tag2.soma_msr,1);
tag2.cell_msr = measure(tag2_cell_lb,tag2.image,{'size','sum'});
tag2.cell_numpuncta = size(tag2.cell_msr,1);

%write tag1
tag1.table = table(tag1.total_numpuncta,tag1.soma_numpuncta,tag1.cell_numpuncta);
tag1.table.Properties.VariableNames = {'Tag1_NumTotalPuncta','Tag1_NumSomaticPuncta', 'Tag1_NumCellPuncta'};
tag1TOT_resultsmat = table(tag1.total_msr.Size', tag1.total_msr.Sum');
tag1TOT_resultsmat.Properties.VariableNames = {'Tag1_TotalImage_Sizes','Tag1_TotalImage_Sums'};
tag1SOMA_resultsmat = table(tag1.soma_msr.Size', tag1.soma_msr.Sum');
tag1SOMA_resultsmat.Properties.VariableNames = {'Tag1_Soma_Sizes','Tag1_Soma_Sums'};
tag1CELL_resultsmat = table(tag1.cell_msr.Size', tag1.cell_msr.Sum');
tag1CELL_resultsmat.Properties.VariableNames = {'Tag1_Cell_Sizes','Tag1_Cell_Sums'};
writetable(tag1.table,savename,'Sheet',1,'FileType','spreadsheet','Range','A1');
writetable(tag1SOMA_resultsmat,savename,'Sheet','Tag1 Results','FileType','spreadsheet','Range','A1');
writetable(tag1CELL_resultsmat,savename,'Sheet','Tag1 Results','FileType','spreadsheet','Range','C1');
writetable(tag1TOT_resultsmat,savename,'Sheet','Tag1 Results','FileType','spreadsheet','Range','E1');

%write tag1 colocalizing with tag2
tag1.colocWtag2.table = table(tag1.colocWtag2.total_numpuncta,tag1.colocWtag2.soma_numpuncta,tag1.colocWtag2.cell_numpuncta);
tag1.colocWtag2.table.Properties.VariableNames = {'Tag1_colocWtag2_NumTotalPuncta','Tag1_colocWtag2_NumSomaticPuncta', 'Tag1_colocWtag2_NumCellPuncta'};
tag1_colocWtag2_TOT_resultsmat = table(tag1.colocWtag2.total_msr.Size', tag1.colocWtag2.total_msr.Sum');
tag1_colocWtag2_TOT_resultsmat.Properties.VariableNames = {'Tag1_colocWtag2_TotalImage_Sizes','Tag1_colocWtag2_TotalImage_Sums'};
tag1_colocWtag2SOMA_resultsmat = table(tag1.colocWtag2.soma_msr.Size', tag1.colocWtag2.soma_msr.Sum');
tag1_colocWtag2SOMA_resultsmat.Properties.VariableNames = {'Tag1_colocWtag2_Soma_Sizes','Tag1_colocWtag2_Soma_Sums'};
tag1_colocWtag2CELL_resultsmat = table(tag1.colocWtag2.cell_msr.Size', tag1.colocWtag2.cell_msr.Sum');
tag1_colocWtag2CELL_resultsmat.Properties.VariableNames = {'Tag1_colocWtag2_Cell_Sizes','Tag1_colocWtag2_Cell_Sums'};

writetable(tag1.colocWtag2.table,savename,'Sheet',1,'FileType','spreadsheet','Range','A7');
writetable(tag1_colocWtag2SOMA_resultsmat,savename,'Sheet','Tag1_colocWtag Results','FileType','spreadsheet','Range','A1');
writetable(tag1_colocWtag2CELL_resultsmat,savename,'Sheet','Tag1_colocWtag Results','FileType','spreadsheet','Range','C1');
writetable(tag1_colocWtag2_TOT_resultsmat,savename,'Sheet','Tag1_colocWtag Results','FileType','spreadsheet','Range','E1');

%write tag2
tag2.table = table(tag2.total_numpuncta,tag2.soma_numpuncta,tag2.cell_numpuncta);
tag2.table.Properties.VariableNames = {'Tag2_NumTotalPuncta','Tag2_NumSomaticPuncta', 'Tag2_NumCellPuncta'};
tag2TOT_resultsmat = table(tag2.total_msr.Size', tag2.total_msr.Sum');
tag2TOT_resultsmat.Properties.VariableNames = {'Tag2_TotalImage_Sizes','Tag2_TotalImage_Sums'};
tag2SOMA_resultsmat = table(tag2.soma_msr.Size', tag2.soma_msr.Sum');
tag2SOMA_resultsmat.Properties.VariableNames = {'Tag2_Soma_Sizes','Tag2_Soma_Sums'};
tag2CELL_resultsmat = table(tag2.cell_msr.Size', tag2.cell_msr.Sum');
tag2CELL_resultsmat.Properties.VariableNames = {'Tag2_Cell_Sizes','Tag2_Cell_Sums'};
writetable(tag2.table,savename,'Sheet',1,'FileType','spreadsheet','Range','A4');
writetable(tag2SOMA_resultsmat,savename,'Sheet','Tag2 Results','FileType','spreadsheet','Range','A1');
writetable(tag2CELL_resultsmat,savename,'Sheet','Tag2 Results','FileType','spreadsheet','Range','C1');
writetable(tag2TOT_resultsmat,savename,'Sheet','Tag2 Results','FileType','spreadsheet','Range','E1');
% write soma size
areatable = table(soma.area,sum(cellfill.mask(:)));
areatable.Properties.VariableNames = {'Soma_Num_Pixels', 'Cell_Mask_Num_Pixels'};
writetable(areatable,savename,'Sheet',1,'FileType','spreadsheet','Range','D1');
% save the results as .mat file
save(savename,'cellfill','tag1','tag2','soma');
%% --- example commands

dipshow(cellfill.image,'lin')
dipshow(tag1.image,'lin')
dipshow(tag2.image,'lin');
joinchannels('rgb',tag1.mask,tag2.mask,cellfill.mask)
joinchannels('rgb',soma.mask,tag1.colocWtag2.	mask)
GeneralAnalysis.viewMaskOverlay(cellfill.image,tag2.mask)
