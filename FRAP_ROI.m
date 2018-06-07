classdef FRAP_ROI < handle
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
           normROI = corrROI./mean(corrROI(1:(obj.bleachframes(1)-1))); %normalize intensity to pre-bleaching frame intensity
           obj.RecoveryCurve.x = single(squeeze(obj.time));
           obj.RecoveryCurve.y = single(squeeze(normROI));
       end
       function f = plotRecoveryCurve(obj)
            f = figure;
           plot(obj.RecoveryCurve.x,obj.RecoveryCurve.y);
%            ylim([0 1]);
           xlabel('Time (seconds)');
           ylabel('Relative Fluorescence');
       end
       function fitRecoveryCurve_1component(obj)
           if isempty(obj.RecoveryCurve)
               obj.makeRecoveryCurve;
           end
       end
   end
end