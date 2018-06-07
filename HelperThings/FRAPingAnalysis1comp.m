function [h1 frapvars] = FRAPingAnalysis1comp(time,blfrm,bleachcontrol,rois,background)
 % this function takes in arrays for all FRAP info, plots and calculates
 % IF,tau etc for 1 component.
% time- [n x 1] in seconds for n data points/timesteps
% blfrm- frames that were used for bleaching can be number or vector;
% bleach control- [n x 1] non-bleached ROI over time to determine overall
%       photobleaching in image.
% ROIs- [n x m] where m is number of bleached ROI(s) for FRAP analysis (number or array for multiple ROIs)
% background- [n x 1] optional off cell ROI to determine background subtraction.

startblframe=min(blfrm);
stopblframe=max(blfrm);
cmp=lines(size(rois,2));
h1=figure;
lgd={};
for jj=1:size(rois,2)
    bl=rois(:,jj)-background;
    pbcor=bleachcontrol-background;
    %Correct for photobleaching
    blcor=bl./(pbcor./pbcor(1));
    %plot(pb./(pbcor./pbcor(1)));
    %Remove pre-bleach to start fitting after bleach
    blcorRm=blcor(stopblframe+1:end);
    timeRm=time(stopblframe+1:end)-time(stopblframe+1);
    x=squeeze(double(timeRm));
    y=squeeze(double(blcorRm));
    preint=mean(blcor(1:startblframe-1));
    postint=blcor(startblframe);
    opts = fitoptions('method','NonlinearLeastSquares','Lower',[0,0,-Inf]);
    
    %First plot and fit for single exponential
    yend=y(1)*2;
    opts.StartPoint=[y(1),yend,5];
    opts.Lower=[y(1),1,0.1];
    ftype =fittype('y0+A1*(1-exp(-x/tau1))','options',opts,'coeff',{'y0','A1','tau1'},'indep','x','depen','y');
    [results,goodness] = fit(x,y,ftype)
    c=coeffvalues(results)
    y0(jj,1)=c(1);
    A1(jj,1)=c(2);
    tau1(jj,1)=c(3);
    I0(jj,1)=preint;
    IF(jj,1) = 1 - ((A1(jj,1))/(preint-postint))
    A2(jj,1)=0;
    tau2(jj,1)=0;
    hold on
    plot(squeeze(x),squeeze(blcorRm),'Marker','o','Color', cmp(jj,:),'LineStyle','none')
    fity=[c(1)+c(2).*(1-exp(-x./c(3)))];
    plot(x,fity,'Color', cmp(jj,:),'LineWidth',2 )
    xlabel('Time (s)','FontSize',14);
    ylabel('Intensity','FontSize',14);
    lgd{(jj-1)*2+1}=['Data from Region ' num2str(jj)];
    lgd{(jj-1)*2+2}=['tau = ', num2str(c(3)) ', IF = ', num2str(IF(jj,1))];
    title('FRAP: Single Exponential Fit','FontSize',14);
    
    
end
legend(lgd);

%  now save output info
frapvars.y0=y0;
frapvars.A1=A1;
frapvars.tau1=tau1;
frapvars.I0=I0;
frapvars.IF=IF;
frapvars.x=x;
frapvars.blcorRm=blcorRm;
frapvars.fity=fity;
frapvars.legend=lgd;
frapvars.results=results;
frapvars.goodness=goodness;


end


%% old to move
% this function takes in a list of text file(s) and a cell array of indices
% designating which columns should be used for - time, bleach control,
% rois, and background