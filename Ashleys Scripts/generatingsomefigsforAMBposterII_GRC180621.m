%% plotting glua1 vs nl1

% distance 50-100
nl1_100 = []; % copy in NL1 file
nl1_100(nl1_100 ==0) = NaN;

nl1_100mean = nanmean(nl1_100(2:end,:));
nl1_100counts = sum(~isnan(nl1_100(2:end,:)),1);
nl1_100std = nanstd(nl1_100(2:end,:))./sqrt(nl1_100counts);

% figure; plot(nl1_100(1,:),nl1_100mean)
figure;
errorbar(nl1_100(1,:),nl1_100mean,nl1_100std); hold on;
plot(nl1_100(1,:),nl1_100mean);

glu1_100 = []; % copy in NL1 file
glu1_100(glu1_100 ==0) = NaN;

glu1_100mean = nanmean(glu1_100(2:end,:));
glu1_100counts = sum(~isnan(glu1_100(2:end,:)),1);
glu1_100std = nanstd(glu1_100(2:end,:))./sqrt(glu1_100counts);

% figure; plot(glu1_100(1,:),glu1_100mean)
errorbar(glu1_100(1,:),glu1_100mean,glu1_100std)
plot(glu1_100(1,:),glu1_100mean)

xlabel('Time (min)');
ylabel('Relative Intensity (AU)') 
yL = get(gca,'YLim');
line([0 0],yL,'Color','k');
%%
% distance 50-100
nl1_200 = []; % copy in NL1 file
nl1_200(nl1_200 ==0) = NaN;

nl1_200mean = nanmean(nl1_200(2:end,:));
nl1_200counts = sum(~isnan(nl1_200(2:end,:)),1);
nl1_200std = nanstd(nl1_200(2:end,:))./sqrt(nl1_200counts);

% figure; plot(nl1_200(1,:),nl1_200mean)
figure;
errorbar(nl1_200(1,:),nl1_200mean,nl1_200std); hold on;
plot(nl1_200(1,:),nl1_200mean);

gluA1_200 = []; % copy in NL1 file
gluA1_200(gluA1_200 ==0) = NaN;

glu1_200mean = nanmean(gluA1_200(2:end,:));
glu1_200counts = sum(~isnan(gluA1_200(2:end,:)),1);
glu1_200std = nanstd(gluA1_200(2:end,:))./sqrt(glu1_200counts);

% figure; plot(glu1_200(1,:),glu1_200mean)
errorbar(gluA1_200(1,:),glu1_200mean,glu1_200std)
plot(gluA1_200(1,:),glu1_200mean)


xlabel('Time (min)');
ylabel('Relative Intensity (AU)') 
yL = get(gca,'YLim');
line([0 0],yL,'Color','k');
%%
%%
% distance 0 50
nl1_50 = []; % copy in NL1 file
nl1_50(nl1_50 ==0) = NaN;

nl1_50mean = nanmean(nl1_50(2:end,:));
nl1_50counts = sum(~isnan(nl1_50(2:end,:)),1);
nl1_50std = nanstd(nl1_50(2:end,:))./sqrt(nl1_50counts);

% figure; plot(nl1_50(1,:),nl1_50mean)
figure;
errorbar(nl1_50(1,:),nl1_50mean,nl1_50std); hold on;
plot(nl1_50(1,:),nl1_50mean);

gluA1_50 = []; % copy in NL1 file
gluA1_50(gluA1_50 ==0) = NaN;

glu1_50mean = nanmean(gluA1_50(2:end,:));
glu1_50counts = sum(~isnan(gluA1_50(2:end,:)),1);
glu1_50std = nanstd(gluA1_50(2:end,:))./sqrt(glu1_50counts);

% figure; plot(glu1_200(1,:),glu1_200mean)
errorbar(gluA1_50(1,:),glu1_50mean,glu1_50std)
plot(gluA1_50(1,:),glu1_50mean)


xlabel('Time (min)');
ylabel('Relative Intensity (AU)') 
yL = get(gca,'YLim');
line([0 0],yL,'Color','k');