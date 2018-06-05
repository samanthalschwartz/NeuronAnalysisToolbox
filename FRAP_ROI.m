classdef FRAP_ROI < handle
    properties
        Intensity = []; % n x 1 vector of intensities, where n is number of frames 
        time = [];
        Background = [];
        Control = [];
        bleachframes = [];
        RecoveryCurve = [];
        FitResults = []; 
   end
   methods 
       function makeRecoveryCurve(obj)
           rROI = obj.Intensity - obj.Background; % remove background from ROI
           control = obj.Control - obj.Background; %remove background from control
           corrROI = rROI./(control./control(1)); % adjust ROI intensity relative to control bleaching
           normROI = corrROI./mean(corrROI(1:(obj.bleachframes(1)-1))); %normalize intensity to pre-bleaching frame intensity
           obj.RecoveryCurve.x = obj.time;
           obj.RecoveryCurve.y = normROI;
       end
       function fitRecoveryCurve_1component(obj)
           if isempty(obj.RecoveryCurve)
               obj.makeRecoveryCurve;
           end
       end
   end
end