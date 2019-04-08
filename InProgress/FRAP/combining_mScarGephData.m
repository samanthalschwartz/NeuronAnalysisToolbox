% load and average GephIB files
files = uipickfiles();
allnormvals = [];
for ii = 1:numel(files)
   currdat = load(files{ii}); 
   allnormvals = [allnormvals, mean(currdat.normvals,2)];
end
gphib = mean(allnormvals,2);
gphibst = std(allnormvals,0,2)/sqrt(4)

cry2files = uipickfiles();
cry2allnormvals = [];
for ii = 1:numel(cry2files)
   currdat_cry2 = load(cry2files{ii}); 
   cry2allnormvals = [cry2allnormvals,  mean(currdat_cry2.normvals,2)];
end
cry2 = mean(cry2allnormvals,2);
cry2st = std(cry2allnormvals,0,2)/sqrt(3);

figure;errorbar(cry2,cry2st);
hold on;
errorbar(gphib,gphibst)

%%
cry2 = load('G:\FromMicroscopeComputer\190303 mScarGeph FRAP_FingR\GephIB\cell1highexp_mScarGeph_GephIBGFP_20190303_70622 PM\Results.mat');
cry2.dataim
for ii = 1:size(cry2.frappa,1)
    rect = rectangle('Position',[(cry2.frappa(ii,1))-ceil(cry2.boxsize/2),(cry2.frappa(ii,2))-ceil(cry2.boxsize/2),cry2.boxsize,cry2.boxsize],...
        'EdgeColor','red',...
        'LineWidth',1);
end
for ii = 1:numel(cry2.ub_rois)
   rect_unbleached = rectangle('Position',[cry2.ub_rois{ii}(1,1),cry2.ub_rois{ii}(1,2),cry2.ub_rois{ii}(2,1),cry2.ub_rois{ii}(2,2)],...
        'EdgeColor','blue',...
        'LineWidth',1); 
end
rect_background = rectangle('Position',[cry2.C1(1,1),cry2.C1(1,2),cry2.C1(2,1),cry2.C1(2,2)],...
                    'EdgeColor','green',...
                    'LineWidth',1);
%%                
cry2 = load('G:\FromMicroscopeComputer\190303 mScarGeph FRAP_FingR\Cry2Olig_GephIB\cell1_mScarGeph_Cry2OligGephIBGFP_20190303_92949 PM\Results.mat');
cry2.dataim
for ii = 1:size(cry2.frappa,1)
    rect = rectangle('Position',[(cry2.frappa(ii,1))-ceil(cry2.boxsize/2),(cry2.frappa(ii,2))-ceil(cry2.boxsize/2),cry2.boxsize,cry2.boxsize],...
        'EdgeColor','red',...
        'LineWidth',1);
end
for ii = 1:numel(cry2.ub_rois)
   rect_unbleached = rectangle('Position',[cry2.ub_rois{ii}(1,1),cry2.ub_rois{ii}(1,2),cry2.ub_rois{ii}(2,1),cry2.ub_rois{ii}(2,2)],...
        'EdgeColor','blue',...
        'LineWidth',1); 
end
rect_background = rectangle('Position',[cry2.C1(1,1),cry2.C1(1,2),cry2.C1(2,1),cry2.C1(2,2)],...
                    'EdgeColor','green',...
                    'LineWidth',1);
