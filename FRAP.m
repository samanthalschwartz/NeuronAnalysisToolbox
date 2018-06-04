classdef FRAP < handle
    properties
       ROI = []; % n x 1 array of intensity within the bleached region time series
       bleachframes = []; % m x 1 array of frames where photobleaching occured (m<n). First frame = 1. Example: Bleaching from frames 10 - 15 would be the array [10:15]; 
       framerate = 1; % frame rate in seconds/frame
       backgroundROI = 0; % n x 1 array representing the intensity of a region outside the cell used for background offset. Can also be a single value b if it is the same over all frames.
       controlROI = 1; % n x 1 array of intensity for a region within the cell but outside the bleached region used to correct for photobleaching over the recovery period
       RecoveryCurve = [];
       FitParams = [];
    end
    
    methods
        function makeRecoveryCurve(obj)
            rROI = obj.ROI - obj.backgroundROI; % remove background from ROI
            control = obj.controlROI - obj.backgroundROI; %remove background from control
            corrROI = rROI./(control./control(1)); % adjust ROI intensity relative to control bleaching
            normROI = corrROI./mean(corrROI(1:(obj.bleachframes(1)-1))); %normalize intensity to pre-bleaching frame intensity
            time = obj.framerate * (1:numel(obj.ROI)); % set array of times
            obj.RecoveryCurve.x = time;
            obj.RecoveryCurve.y = normROI;
        end
        function fitRecoveryCurve_1component(obj)
            if isempty(obj.RecoveryCurve)
                obj.makeRecoveryCurve
            end
        end
    end
    methods (Static)
        function fit1componentexponential(x,y)
            
        end
        function h = plotRecoveryCurve(curves)
            h = figure(); hold on;
            for ii = 1:numel(curves)
                plot(h,curves{ii}.x,curves{ii}.y);
                xlabel = 'Time (seconds)';
                ylabel = 'Normalized Fluorescence Intensity';
            end
        end
        function image = loadtiffSeries(datadir,uniquefilestr,zproject)
            if nargin<3
                zproject = 'sum';
            end
            files = dir(fullfile(datadir,uniquefilestr));
            frame1 = loadtiff(fullfile(datadir,files(1).name));
            image = dip_image(zeros(size(frame1,1),size(frame1,2),numel(files)));
            switch zproject
                case 'sum'
                    for ff = 1:numel(files)
                        frame = loadtiff(fullfile(datadir,files(ff).name));
                        image(:,:,ff-1) = sum(dip_image(frame),[],3);
                    end
                case 'max'
                    for ff = 1:numel(files)
                        frame = loadtiff(fullfile(datadir,files(ff).name));
                        image(:,:,ff-1) = max(dip_image(frame),[],3);
                    end
            end
        end
        function getmeanROIfromCoords(img_in,coords)
            
        end
        function [ch1, ch2] = loadAndorFile(datadir,unique_ch1,unique_ch2)
            if nargin < 1 || isempty(datdir)
                datadir = uigetdir('.','Select the folder containing the .tif series for Channel1 and Channel2');
            end
            if nargin < 2
            unique_ch1 = '*t*_w0000.tif';
            unique_ch2 = '*t*_w0001.tif';
            end
            if nargin == 2 % -- only need to load in 1 channel of data series
                ch1 = FRAP.loadtiffSeries(datadir,unique_ch1,'sum');
                return;
            end
            ch1 = FRAP.loadtiffSeries(datadir,unique_ch1,'sum');
            ch2 = FRAP.loadtiffSeries(datadir,unique_ch2,'sum');
        end
    end
end