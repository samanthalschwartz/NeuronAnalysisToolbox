classdef FRAP < handle
    properties
        %         fp.ROIs = {[328, 186],[209, 241],[211, 302],[332, 335]};
        % fp.ROIs = {[305, 205],[49, 38],[330, 317],[407, 294]};
        filestr = {'*t*_w0000.tif','*t*_w0001.tif','*.tif'};
        image = [];
        framerate = [];
        background = [];
        backgroundROI = [];
        control = [];
        controlROI = [];
        FRAP_ROI = [];
        %        Ch = struct(...
        %            'image',[],...
%            'framerate',1,... % frame rate in seconds/frame
%            'background',0,...
%            'control',[],...
%            'backgroundROI',[],...
%            'controlROI',[],...
%            'FRAP_ROI',{});
       ROIs = []; % list of ROIs to be analyzed
       %cell array of FRAP_ROI objects
       bleachframes = []; % m x 1 array of frames where photobleaching occured (m<n). First frame = 1. Example: Bleaching from frames 10 - 15 would be the array [10:15]; 
       inc = 10; %size of bleach square
    end
        
    methods
        function makeRecoveryCurve(obj)
            
        end
        function fitRecoveryCurve_1component(obj)
            if isempty(obj.RecoveryCurve)
                obj.makeRecoveryCurve
            end
        end
        function load(obj,datadir,uniquestr)
            if nargin<3
                uniquestr = obj.filestr{3};
                if nargin<2
                    datadir = [];
                end
            end
            obj.image = loadtiffSeries(datadir,uniquestr,'sum');
        end
        function getBackgroundROI(obj)
            h = dipshow(sum(obj.image,[],3),'lin');
            diptruesize(h,100);
            [B,C] = dipcrop(h);
            obj.backgroundROI.image = B;
            obj.backgroundROI.coords = C;
            close(h);
        end
        function h = showROIs(obj)
           h = dipshow(obj.image,'log');
            diptruesize(h,100);
            for pp = 1:numel(obj.ROIs)
                test = rectangle('Position',[obj.ROIs{pp}(1)-obj.inc,obj.ROIs{pp}(2)-obj.inc,10,10],...
                    'EdgeColor','w',...
                    'LineWidth',1);
                %                 roi{pp} = data(obj.ROIs{pp}(1)-inc:obj.ROIs{pp}(1)+inc,obj.ROIs{pp}(2)-inc:obj.ROIs{pp}(2)+inc,:);
            end            
        end
        function getControlROI(obj)
            h = dipshow(sum(obj.Ch.image,[],3),'lin');
            diptruesize(h,100);
            for pp = 1:numel(obj.ROIs)
                test = rectangle('Position',[obj.ROIs{pp}(1)-obj.inc,obj.ROIs{pp}(2)-obj.inc,10,10],...
                    'EdgeColor','w',...
                    'LineWidth',1);
                %                 roi{pp} = data(obj.ROIs{pp}(1)-inc:obj.ROIs{pp}(1)+inc,obj.ROIs{pp}(2)-inc:obj.ROIs{pp}(2)+inc,:);
            end
            [B,C] = dipcrop(h);
            obj.controlROI.image = B;
            obj.controlROI.coords = C;
            close(h);
        end
        function makeRecoveryCurves(obj)
          %---- Channel1 
            % make background vector
            % make control vector
          %---- Channel2 
            % make background vector
            % make control vector
            
           %for each ROI load in and make curve 
           Ch_FRAP_ROIs = []; %cell array of FRAP_ROIs for ch1
           ch2_FRAP_ROIs = []; %cell array of FRAP_ROIs for ch1
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
        function out = loadAndorFile(datadir,unique_ch1,unique_ch2)
            if nargin < 1 || isempty(datadir)
                datadir = uigetdir('.','Select the folder containing the .tif series for Channel1 and Channel2');
            end
            if nargin < 2
            unique_ch1 = '*t*_w0000.tif';
            unique_ch2 = '*t*_w0001.tif';
            end
            if nargin == 2 % -- only need to load in 1 channel of data series
                out.ch1 = FRAP.loadtiffSeries(datadir,unique_ch1,'sum');
                return;
            end
            out.ch1 = FRAP.loadtiffSeries(datadir,unique_ch1,'sum');
            out.ch2 = FRAP.loadtiffSeries(datadir,unique_ch2,'sum');
        end
    end
end

% %%
% if isempty(obj.Background)
%                obj.Background = obj.calcROImean(obj.BackgroundROI);
%            end
%             if isempty(obj.Control)
%                obj.Control = obj.calcROImean(obj.ControlROI);
%             end
%            
%             methods (Static)
%         function ROImean = calcROImean(in_ROI)
%             im = dip_image(in_ROI);
%             ROImean = single(squeeze(mean(im,[],[1 2])));            
%         end
%         function ROIsum = calcROIsum(in_ROI)
%             im = dip_image(in_ROI);
%             ROIsum = single(squeeze(sum(im,[],[1 2])));  
%         end 
