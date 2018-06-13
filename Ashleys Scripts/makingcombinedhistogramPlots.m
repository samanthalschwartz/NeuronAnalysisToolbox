datafilepath = 'F:\Sam\Data\Ashley\180607 Results\histogram distances';
GluA1files= {'HistogramResults-DHFR-GFP-GluA1-051718-slip1',...
    'HistogramResults-GluA1release_REs-060118-slip1',...
    'HistogramResults-GluA1release_REs-060118-slip2'};
%     'HistogramResults-DHFR-GFP-GluA1-051718-slip2',...
    

NL1files= {'HistogramResults-050118-NL1-insertion',...
    'HistogramResults-051318-NL1-insertion'};

distances = [0:30:120];
% collect all the results.....

% loop through and make matrix for each distance
glu11distances = cell(1,numel(distances));
for dd = 1:numel(distances) %append data to matrix corresponding to distance plot dd
    distancematrix = nan(50,4); % make too big and then shrink down later
    distmatid = 1;
    for gg = 1:numel(GluA1files) %loop through day by day to extract for distances dd
        aa = load(fullfile(datafilepath,GluA1files{gg}));
        for rr = 1:size(aa.savedhistograms,1) %loop through each file within the day to add to distancematrix
            temphist = aa.savedhistograms{rr,dd};
            distancematrix(distmatid,:) = temphist.yvals;
            distmatid = distmatid + 1;
        end
    end
    glu11distances{dd} = distancematrix(~isnan(distancematrix(:,1)),:);
end
NL1distances = cell(1,numel(distances));
for dd = 1:numel(distances) %append data to matrix corresponding to distance plot dd
    distancematrix = nan(50,4); % make too big and then shrink down later
    distmatid = 1;
    for gg = 1:numel(NL1files) %loop through day by day to extract for distances dd
        aa = load(fullfile(datafilepath,NL1files{gg}));
        for rr = 1:size(aa.savedhistograms,1) %loop through each file within the day to add to distancematrix
            temphist = aa.savedhistograms{rr,dd};
            distancematrix(distmatid,:) = temphist.yvals;
            distmatid = distmatid + 1;
        end
    end
    NL1distances{dd} = distancematrix(~isnan(distancematrix(:,1)),:);
end
%%
xplotvals = temphist.xvals;
f = figure; hold on;
for tt = 1:numel(distances)
    subplot(numel(distances),1,tt); hold on;
    currGlu1 = glu11distances{tt};
    currNL1 = NL1distances{tt};
    
    currGlu1(currGlu1==0) = nan;
    currGlu1_meanvals = nanmean(currGlu1);
%     currGlu1_stdvals = nanstd(currGlu1);
    currGlu1_stdvals = nanstd(currGlu1)./sqrt(sum(~isnan(currGlu1)));
    errorbar(xplotvals,currGlu1_meanvals,currGlu1_stdvals,'DisplayName','GluA1','LineWidth',1.5); hold on;    
    currNL1(currNL1==0) = nan;
    currNL1_meanvals = nanmean(currNL1);
%     currNL1_stdvals = nanstd(currNL1);
    currNL1_stdvals = nanstd(currNL1)./sqrt(sum(~isnan(currNL1)));
    errorbar(xplotvals,currNL1_meanvals,currNL1_stdvals,'DisplayName','NL1','LineWidth',1.5);
    
    if tt ==1
        l = legend(); set(l,'Position', [0.1451    0.7968    0.1535    0.1050])
    end
    if tt < numel(distances)
        titlestr = (['Events Occuring between '...
            num2str(distances(tt)) ' - ' num2str(distances(tt+1)) ' \mum of Cell Soma']);
    else
        titlestr = (['Events Occuring at greater than '...
            num2str(distances(tt)) ' \mum of Cell Soma']);
    end
    title(titlestr)
%     ylim([0 0.2]);
    xlabel('Time (min)','FontSize',10);
    ylabel('% Total Cargo Events','FontSize',10);

end
set(gcf,'Position',[1          41         912        1040]);
% suptitle(['Cargo Delivery For File: ' datapath])


