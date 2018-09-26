classdef SIM < handle
    properties
        filepath = '';        
        channelordering = [2 1 3]; % array indicating which channel is which;
        %         channelordering(1) = abeta channel in image
        %         channelordering(2) = channel in image corresponding to 'ch1';
        %         channelordering(3) = channel in image corresponding to 'ch2';
        channelorderingstr = {'chABeta','ch1','ch2'};
        planeTOP = [];
        planeBOTTOM = [];
        abeta = struct('image',[],'mask',[],'distance_mask',[],'labeled_mask',[],'msr',[],'sizes',[],'densities',[]);
        ch1 = struct('image',[],'distance_mask',[],'mask',[],'name','',...
            'thisCh',struct('msr',[],'sizes',[],'densities',[]),...
            'abetaCh',struct('msr',[],'sizes',[],'densities',[],'radialdensity_norm',[],'radialdensity_raw',[],...
                        'radialnumberdensity_norm',[],'radialnumberdensity_raw',[],...
                        'cumulative_radialnumberdensity_norm',[],'cumulative_radialnumberdensity_raw',[]));
        ch2 = struct('image',[],'distance_mask',[],'mask',[],'name','',...
            'thisCh',struct('msr',[],'sizes',[],'densities',[]),...
            'abetaCh',struct('msr',[],'sizes',[],'densities',[],'radialdensity_norm',[],'radialdensity_raw',[],...
                        'radialnumberdensity_norm',[],'radialnumberdensity_raw',[],...
                        'cumulative_radialnumberdensity_norm',[],'cumulative_radialnumberdensity_raw',[]));
        measurements = {'size','sum'};
        XYpxsize = 0.0321;
        XYpxsize_units = '�m';
        Zpxsize = 0.2;
        Zpxsize_units = '�m'; 
        metadata = [];
        cellmask = [];
    end
    
    methods
        function loadNDfile(obj,filepath)
            if nargin>1
                obj.filepath = filepath;
            else
                if isempty(obj.filepath)
                    filepath = uipickfiles('Prompt','Pick a .nd File');
                    obj.filepath = filepath;
                end
            end
            if iscell(filepath)
                filepath = filepath{1};
            end
            try
            [image,obj.metadata] = obj.ndFileloader(obj.filepath); % load file
            catch
            disp('File did not load - make sure it is an ND file?');
            end
            if size(image,4)<3  % assign channels
                obj.abeta.rawimage = image(:,:,:,obj.channelordering(1));
                obj.ch1 = image(:,:,:,obj.channelordering(2));
            else
                obj.abeta.image = image(:,:,:,obj.channelordering(1));
                obj.ch1.image = image(:,:,:,obj.channelordering(2));
                obj.ch2.image = image(:,:,:,obj.channelordering(3));
            end
            %-- get some meta data info
            obj.XYpxsize = double(obj.metadata.getPixelsPhysicalSizeX(0).value());           % returns value in default unit
            obj.XYpxsize_units = char(obj.metadata.getPixelsPhysicalSizeX(0).unit().getSymbol());
            obj.Zpxsize = double(obj.metadata.getPixelsPhysicalSizeZ(0).value());           % returns value in default unit
            obj.Zpxsize_units = char(obj.metadata.getPixelsPhysicalSizeZ(0).unit().getSymbol());
        end
        function make_cellmask(obj)
            sumim = obj.ch1.image+obj.ch2.image+obj.abeta.image;
            maxsumim = max(sumim,[],3);
            gmaxsumim = gaussf(maxsumim,10);
            %  threshs = multithresh(single(gmaxsumim),2);
            %  cellmask = gmaxsumim>threshs(1);
            openeninimage_out = opening(gmaxsumim,20,'rectangular');
            threshs = multithresh(single(openeninimage_out),2);
            msk_open = openeninimage_out>threshs(1);
            obj.cellmask = repmat(msk_open,1,1,size(obj.ch1.image,3));
        end
        function make_masks(obj)
            disp('Generating Whole Cell Mask');
            obj.make_cellmask;
            disp('Generating ABeta Mask');
            obj.make_maskchAB;
            disp('Generating CH1 Mask');
            obj.make_maskch1;
            disp('Generating CH2 Mask');
            obj.make_maskch2;
        end
        function make_distancemasks(obj,maxval)
            if nargin<2
                maxval = 40;
            end
            zscale = obj.Zpxsize/obj.XYpxsize;
            disp('Calculating ABeta Distance');
            obj.abeta.distance_mask = bwdistsc1(single(obj.abeta.mask),[1 1 zscale],maxval);
            disp('Calculating Ch1 Distance');
            obj.ch1.distance_mask = bwdistsc1(single(obj.ch1.mask),[1 1 zscale],maxval);
            disp('Calculating Ch2 Distance');
            obj.ch2.distance_mask = bwdistsc1(single(obj.ch2.mask),[1 1 zscale],maxval);
        end
        function calculateDensities_abetaINch1(obj,bins)
            if nargin<2
                bins = 0:30;
            end
            obj.ch1.abetaCh = obj.calculateDensities_abetaINch(obj.ch1,obj.abeta.mask,obj.abeta.image,obj.cellmask,bins);
        end
        
        function calculateDensities_abetaINch2(obj,bins)
            if nargin<2
                bins = 0:30;
            end
            obj.ch2.abetaCh = obj.calculateDensities_abetaINch(obj.ch2,obj.abeta.mask,obj.abeta.image,obj.cellmask,bins);
        end
        %-- helpers for making masks. if you want to change how the masks
        %are made then change which static method these are calling
        function make_maskchAB(obj)
            premask  = SIM.make_maskABeta(obj.abeta.image);
            m = zeros(1,size(obj.cellmask,3));
            image = dip_image(obj.abeta.image);
            for ii = 0:(size(obj.cellmask,3)-1)
                currframe = image(:,:,ii);
                currcellmask = ~obj.cellmask(:,:,ii);
                currcellmask = bdilation(currcellmask,2);
                m(ii+1) = single(sum(currframe,currcellmask,[1 2]))./sum(currcellmask);
                premask(:,:,ii) = premask(:,:,ii).*~(currframe<m(ii+1));
            end
            obj.abeta.mask = premask;
            obj.abeta.labeled_mask = label(obj.abeta.mask,1);
        end
        function make_maskch1(obj)
            if strcmp(obj.channelorderingstr{2},'PSD95i')
                premask  = SIM.make_maskPSD95ib(obj.ch1.image);
            else
                premask  = SIM.make_maskSynapseMarker(obj.ch1.image);
            end
            m = zeros(1,size(obj.cellmask,3));
            image = dip_image(obj.ch1.image);
            for ii = 0:(size(obj.cellmask,3)-1)
                currframe = image(:,:,ii);
                currcellmask = ~obj.cellmask(:,:,ii);
                currcellmask = bdilation(currcellmask,2);
                m(ii+1) = single(sum(currframe,currcellmask,[1 2]))./sum(currcellmask);
                premask(:,:,ii) = premask(:,:,ii).*~(currframe<m(ii+1));
            end
            obj.ch1.mask = premask;
        end
        function make_maskch2(obj)
            premask  = SIM.make_maskSynapseMarker(obj.ch2.image);
            m = zeros(1,size(obj.cellmask,3));
            image = dip_image(obj.ch2.image);
            for ii = 0:(size(obj.cellmask,3)-1)
                currframe = image(:,:,ii);
                currcellmask = ~obj.cellmask(:,:,ii);
                currcellmask = bdilation(currcellmask,2);
                m(ii+1) = single(sum(currframe,currcellmask,[1 2]))./sum(currcellmask);
                premask(:,:,ii) = premask(:,:,ii).*~(currframe<m(ii+1));
            end
            obj.ch2.mask = premask;
        end
        %----
        function measure_allthings(obj)
            disp('Measuring values inside masks');
            obj.measure_AB;
            obj.measure_ch1;
            obj.measure_ch2;
            obj.measure_AB_inch1;
            obj.measure_AB_inch2;
        end
        function measure_AB(obj)
            lb = label(obj.abeta.mask,1);
            msr = measure(lb,obj.abeta.image,obj.measurements);
            obj.abeta.msr = msr;
            obj.abeta.sizes = msr.Size;
            obj.abeta.densities = msr.sum./msr.size;
        end
        function measure_ch1(obj)
            msr = measure(obj.ch1.mask,obj.ch1.image,obj.measurements);
            obj.ch1.thisCh.msr = msr;
            obj.ch1.thisCh.sizes = msr.Size;
            obj.ch1.thisCh.densities = msr.sum./msr.size;
        end
        function measure_ch2(obj)
            msr = measure(obj.ch2.mask,obj.ch2.image,obj.measurements);
            obj.ch2.thisCh.msr = msr;
            obj.ch2.thisCh.sizes = msr.Size;
            obj.ch2.thisCh.densities = msr.sum./msr.size;
        end
        function measure_AB_inch1(obj)
            msr = measure(obj.ch1.mask,obj.abeta.image,obj.measurements);
            obj.ch1.abetaCh.msr = msr;
            obj.ch1.abetaCh.sizes = msr.Size;
            obj.ch1.abetaCh.densities = msr.sum./msr.size;
        end
        function measure_AB_inch2(obj)
            msr = measure(obj.ch2.mask,obj.abeta.image,obj.measurements);
            obj.ch2.abetaCh.msr = msr;
            obj.ch2.abetaCh.sizes = msr.Size;
            obj.ch2.abetaCh.densities = msr.sum./msr.size;
        end
        
        function dothething(obj)
            obj.make_masks();
            obj.make_distancemasks();
            obj.calculateDensities_abetaINch1();
            obj.calculateDensities_abetaINch2();
            obj.measure_allthings();
        end
        
        function save(obj,inputsavedir)
            if nargin<2
                savedir = [obj.filepath(1:end-4) '_SIM'];
            else
                [~,NAME,~] = fileparts(obj.filepath);
                savedir = fullfile(inputsavedir,[NAME '_SIM']);
            end
           save(savedir,'obj'); 
        end
        
    end
    
    methods (Static)
        function mask = make_maskABeta(image,params)
            ab = gaussf(image);
            out = GeneralAnalysis.imgLaplaceCutoff(ab,[2 2 1],[2 2 1]);
%             out = gaussf(out,[2 2 1]);
            thr = multithresh(single(out),2);
            mask_start = out>thr(2);
            sch2 = sum(dip_image(image),[],3);
            sch2_series = repmat(sch2,1,1,size(image,3));
            wshed = GeneralAnalysis.watershed_timeseries(-gaussf(sch2_series,1),1);
            mask = mask_start.*~wshed;
        end
        
        function mask = make_maskSynapseMarker(image,params)
            gch = gaussf(image);
            out = GeneralAnalysis.imgLaplaceCutoff(gch,[2 2 1],[2 2 1]);
            thr = multithresh(single(out),2);
            mask = bdilation(out>thr(2),1);
        end
        
        function mask = make_maskSynapseMarkerLow(image,params)
            gch = gaussf(image);
            out = GeneralAnalysis.imgLaplaceCutoff(gch,[2 2 1],[2 2 1]);
            thr = multithresh(single(out),1);
            mask = bdilation(out>thr(1),1);
        end
        
        function mask = make_maskPSD95ib(image,params)
            gch = gaussf(image);
            out = GeneralAnalysis.imgLaplaceCutoff(gch,[2 2 1],[2 2 1]);
            out = gaussf(out,[2 2 1]);
            thr = multithresh(single(out),2);
            mask = bdilation(out>thr(2),1); 
        end
        
        
        
        function [image,metadata] = ndFileloader(filename)
            if nargin<1
                [name, path] = uigetfile('*.nd2');
                filename = fullfile(path,name);
            end
            hannah = bfopen(filename);
            % use first input to determine #channels, #z-planes
            imsize = size(hannah{1,1}{1,1}); %size of single image frame
            infostr = hannah{1,1}{1,2};
            info = SIM.parseNDtext(infostr);
            image = zeros([imsize,info.totalZ,info.totalCol]);
            for ii = 1:size(hannah{1,1},1)
                infostr = hannah{1,1}{ii,2};
                info = SIM.parseNDtext(infostr);
                image(:,:,info.currZ,info.currCol) = hannah{1,1}{ii,1};
                metadata = hannah{1, 4};
            end
        end
        function info = parseNDtext(infostr)
            
            splstr1 = strsplit(infostr,{';'});
            
            planebool = cellfun(@(x) contains(x,'plane '),splstr1);
            planestr = splstr1{planebool};
            planeparts = strsplit(planestr,{' plane ','/'});
            info.currPlane = str2double(planeparts{2});
            info.totalPlane = str2double(planeparts{3});
            
            zbool = cellfun(@(x) contains(x,'Z=')|contains(x,'Z?='),splstr1);
            zstr = splstr1{zbool};
            zparts = strsplit(zstr,{'=','/'});
            info.currZ = str2double(zparts{2});
            info.totalZ = str2double(zparts{3});
            
            colbool = cellfun(@(x) contains(x,'C=')|contains(x,'C?='),splstr1);
            colstr = splstr1{colbool};
            colparts = strsplit(colstr,{'=','/'});
            info.currCol = str2double(colparts{2});
            info.totalCol = str2double(colparts{3});
        end
        
        function output = calculateDensities_abetaINch(channel,abeta_mask,abeta_image,cellmask,bins)
            nbins = numel(bins);
            %-- alocate empty arrays for values calculated
            density = zeros(nbins,1);
            numb = zeros(nbins,1);
            cumulative_numb = zeros(nbins,1);
%             distmasks = zeros(size(obj.ch2.distance_mask,1),size(obj.ch2.distance_mask,1),nbins);
            %-- label abeta mask to get individuals
            lb = label(abeta_mask,1);
            %-- loop through and calculate some things
            wb = waitbar(0,'Calculating Density Information...');
            for nn = 1:nbins
                cumulative_distmask = (channel.distance_mask<=bins(nn)).*cellmask;
                if nn == 1
                    curr_distmask = (channel.distance_mask<=bins(nn)).*cellmask;
                else
                    curr_distmask = (channel.distance_mask>bins(nn-1) & channel.distance_mask<=bins(nn)).*cellmask;
                end
                density(nn) = sum(abeta_image.*curr_distmask)./sum(curr_distmask);
                inmasklb = single(lb.*curr_distmask);
                numb(nn) = size(unique(inmasklb),1);%./sum(curr_distmask);
                cumulative_numb(nn) = size(unique(lb.*cumulative_distmask));%./sum(curr_distmask); 
                try
                waitbar(nn/nbins,wb);
                catch
                   wb = waitbar(nn/nbins,'Calculating Density Information...');
                end
            end
            close(wb);
            avgdensity = sum(abeta_image.*cellmask)./sum(cellmask);
            incelllb = single(lb.*cellmask);
            avgnumberdensity = size(channel.thisCh.msr,1);
%             avgnumberdensity = size(unique(incelllb),2)./sum(cellmask);
            xs = bins;
            output.radialdensity_norm = [xs',density/avgdensity];
            output.radialdensity_raw = [xs',density];
            output.radialnumberdensity_norm = [xs',numb/avgnumberdensity];
            output.radialnumberdensity_norm_raw = [xs',numb];
            output.cumulative_radialnumberdensity_norm_norm = [xs',cumulative_numb/avgnumberdensity];
            output.cumulative_radialnumberdensity_norm_raw = [xs',cumulative_numb];
        end 
        function results = calculateRadialDensity(mask1_distance,image,cellmask,bins)
            % mask1 is the mask that distances are calculated from (ie: the synapse marker)
            % mask1_distance mask - is distance transform
            % image is the raw image that the intensity will be calculated from (ie: abeta image)
            % cellmask is the whole cell mask that the area will be taken and normalized from
            nbins = numel(bins);
            %-- alocate empty arrays for values calculated
            intensity = zeros(nbins,1);
            volume = zeros(nbins,1);
            %-- loop through and calculate some things
            wb = waitbar(0,'Calculating Intensity Density Information...');
            for nn = 1:nbins
                if nn == 1
                    curr_distmask = (mask1_distance<=bins(nn)).*cellmask;
                else
                    curr_distmask = (mask1_distance>bins(nn-1) & mask1_distance<=bins(nn)).*cellmask;
                end
                intensity(nn) = sum(image.*curr_distmask);
                volume(nn) = sum(curr_distmask); 
                try
                waitbar(nn/nbins,wb);
                catch
                   wb = waitbar(nn/nbins,'Calculating Intensity Density Information...');
                end
            end
            close(wb);
            results.d = bins;
            results.radialintensity = intensity;
            results.cumulativeradialintensity = cumsum(intensity);
            results.volume = volume;
            results.cumulativevolume = cumsum(volume);
            results.totalintensity = sum(image.*cellmask);
            results.totalvolume = sum(cellmask);
            results.radial_density = (results.radialintensity./results.volume)./(results.totalintensity/results.totalvolume);
            results.cumulative_radial_density = (results.cumulativeradialintensity./results.cumulativevolume)./(results.totalintensity/results.totalvolume);
        end
        function results = calculateNumberDensity(mask1_distance,lb,cellmask,bins)
            % mask1_distance mask - is distance transform of mask (ie: synapse marker mask)
            % image is the raw image that the intensity will be calculated from (ie: abeta image)
            % cellmask is the whole cell mask that the area will be taken and normalized from
            % radial_numberdensity is: (#abeta puncta/area)/#gephyrin
            
            nbins = numel(bins);
            %-- alocate empty arrays for values calculated
            numbers = zeros(nbins,1);
            volume = zeros(nbins,1);
            %-- loop through and calculate some things
            wb = waitbar(0,'Calculating Number Density Information...');
            mask10s = mask1_distance==0;
            mask1_lb = label(mask10s,1);
            nummask1 = size(unique(mask1_lb));
            for nn = 1:nbins
                if nn == 1
                    curr_distmask = (mask1_distance<=bins(nn)).*cellmask;
                else
                    curr_distmask = (mask1_distance>bins(nn-1) & mask1_distance<=bins(nn)).*cellmask;
                end
                inmasklb = single(lb.*curr_distmask);
                numbers(nn) = size(unique(inmasklb),1);
                volume(nn) = sum(curr_distmask); 
                try
                waitbar(nn/nbins,wb);
                catch
                   wb = waitbar(nn/nbins,'Calculating Number Density Information...');
                end
            end
            close(wb);
            results.d = bins;
            results.radialnumber = numbers;
            results.cumulativeradialnumber = cumsum(numbers);
            results.volume = volume;
            results.cumulativevolume = cumsum(volume);
            results.totalnumber = size(unique(lb.*cellmask));
            results.totalvolume = sum(cellmask);
            results.radial_numberdensity = results.radialnumber;
            results.cumulative_number_density = results.cumulativeradialnumber;
            results.nummask1 = nummask1;
        end 
        function results = calculateOverlap(mask1_distance,mask2,cellmask,bins)
            % mask1_distance mask - is distance transform of mask (ie: synapse marker mask)
            % mask2 is the mask to find the area overlap from (ie: abeta mask
            % cellmask is the whole cell mask that the area will be taken and normalized from
            % radialareaoverlap is: area of abeta within distance from gephyrin/total area of abeta
            
            nbins = numel(bins);
            %-- alocate empty arrays for values calculated
            inmaskoverlap = zeros(nbins,1);
            %-- loop through and calculate some things
            wb = waitbar(0,'Calculating Number Density Information...');
            for nn = 1:nbins
                if nn == 1
                    curr_distmask = (mask1_distance<=bins(nn)).*cellmask;
                else
                    curr_distmask = (mask1_distance>bins(nn-1) & mask1_distance<=bins(nn)).*cellmask;
                end
                inmaskoverlap(nn) = mask2.*curr_distmask;
                try
                waitbar(nn/nbins,wb);
                catch
                   wb = waitbar(nn/nbins,'Calculating Number Density Information...');
                end
            end
            close(wb);
            results.d = bins;
            results.radialareaoverlap = sum(inmaskoverlap);
            results.cumulativeareaoverlap = cumsum(results.radialareaoverlap);
            results.totalarea = sum(mask2);
            results.radialareaoverlap =  results.radialareaoverlap./results.totalarea;
            results.cumulative_radialareaoverlap = results.cumulativeareaoverlap./results.totalarea;
        end  
    end
end