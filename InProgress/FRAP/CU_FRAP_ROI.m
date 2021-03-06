classdef CU_FRAP_ROI < handle
    properties
        Intensity = []; % n x 2 image of intensities, where n is number of frames 
        time = [];
        Background = [];
        Control = [];
        bleachframes = [];
        RecoveryCurve = [];
        FitResults = []; 
   end
   methods 
       function makeRecoveryCurve(obj)
           background = mean(mean(mean(single(obj.Background))));
           rROI = mean(obj.Intensity-background,[],[1 2]); % remove background from ROI
           control = mean(obj.Control - background,[],[1 2]); %remove background from control
           corrROI = rROI./(control./control(1)); % adjust ROI intensity relative to control bleaching
           normROI = corrROI./mean(squeeze(corrROI(:,:,0:(obj.bleachframes-1)))); %normalize intensity to pre-bleaching frame intensity
           obj.RecoveryCurve.x = single(squeeze(obj.time));
           obj.RecoveryCurve.y = single(squeeze(normROI));
           obj.fitRecoveryCurve_1component;
           
       end
       function f = plotRecoveryCurve(obj)
           obj.fitRecoveryCurve_1component;
            f = figure;
           plot(obj.RecoveryCurve.x,obj.RecoveryCurve.y); hold on;
           x = obj.FitResults.timecorrect;
           fity=[obj.FitResults.y0+obj.FitResults.A1*(1-exp(-x./obj.FitResults.tau1))];
           plot(obj.FitResults.timecorrect+obj.time(obj.bleachframes),fity,'--k');
%            ylim([0 1]);
           xlabel('Time (seconds)');
           ylabel('Relative Fluorescence');
       end
       function fitRecoveryCurve_1component(obj)
           if isempty(obj.RecoveryCurve)
               obj.makeRecoveryCurve;
           end
           blcor = obj.RecoveryCurve.y;
           blcorRm=blcor(obj.bleachframes+1:end);
           timeRm=obj.time(obj.bleachframes+1:end)-obj.time(obj.bleachframes+1);
           x=squeeze(double(timeRm));
           y=squeeze(double(blcorRm));
           preint=mean(blcor(1:obj.bleachframes));
           postint=blcor(obj.bleachframes+1);
           opts = fitoptions('method','NonlinearLeastSquares','Lower',[0,0,-Inf]);
           
           %First plot and fit for single exponential
           yend=y(1)*2;
           opts.StartPoint=[y(1),yend,5];
           opts.StartPoint=[yend,5];
           opts.Lower=[1,0.1];
           opts.Lower=[y(1),1,0.1];
%            ftype = fittype('A*(1-exp(-tau*t))','options',opts,'independent','t','depen','y','coefficients',{'A','tau'});
           ftype =fittype('y0-A1*(exp(-x.*tau1))','options',opts,'coeff',{'y0','A1','tau1'},'indep','x','depen','y');
           ftype =fittype('A1*(exp(-x.*tau1))','options',opts,'coeff',{'A1','tau1'},'indep','x','depen','y');
          [results,goodness] = fit(x',y',ftype);
           c=coeffvalues(results)
           y0=c(1);
           A1=c(2);
           tau1=c(3);
           I0=preint;
           IF = 1 - (-A1)/1-(y0+A1);
           obj.FitResults.timecorrect = x;
           obj.FitResults.coeff = c;
           obj.FitResults.y0 = y0;
           obj.FitResults.A1 = A1;
           obj.FitResults.IF = IF;
           obj.FitResults.tau1 = tau1;
           obj.FitResults.I0 = I0;    
           obj.RecoveryCurve.fitx = x+obj.time(obj.bleachframes);
           obj.RecoveryCurve.fity=y0-A1*(exp(-x.*tau1));
           
       end
   end
end