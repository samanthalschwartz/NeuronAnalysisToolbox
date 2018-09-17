classdef SIM < handle
    properties
        filepath = '';
        channelorderingstr = {'a-beta ch','ch1','ch2'};
        channelordering = [2 1 3]; % array indicating which channel is which; 
            %         channelordering(1) = abeta channel in image
            %         channelordering(2) = channel in image corresponding to 'ch1';
            %         channelordering(3) = channel in image corresponding to 'ch2';
        abeta = struct('image',[],'distance_mask',[],'mask',[],'msr',[],'sizes',[],'densities',[]);
        ch1 = struct('image',[],'distance_mask',[],'mask',[],'name','',...
            'thisCh',struct('msr',[],'sizes',[],'densities',[]),...
            'abetaCh',struct('msr',[],'sizes',[],'densities',[]));
        ch2 = struct('image',[],'distance_mask',[],'mask',[],'name','',...
            'thisCh',struct('msr',[],'sizes',[],'densities',[]),...
            'abetaCh',struct('msr',[],'sizes',[],'densities',[]));
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
            disp('Generating ABeta Mask');
            obj.make_maskchAB;
            disp('Generating CH1 Mask');
            obj.make_maskch1;
            disp('Generating CH2 Mask');
            obj.make_maskch2;
            disp('Generating Whole Cell Mask');
            obj.make_cellmask;
        end
        function make_distancemasks(obj,maxval)
            if nargin<2
                maxval = 30;
            end
            zscale = obj.Zpxsize/obj.XYpxsize;
            disp('Calculating ABeta Distance');
            obj.abeta.distance_mask = bwdistsc1(single(obj.abeta.mask),[1 1 zscale],maxval);
            disp('Calculating Ch1 Distance');
            obj.ch1.distance_mask = bwdistsc1(single(obj.ch1.mask),[1 1 zscale],maxval);
            disp('Calculating Ch2 Distance');
            obj.ch2.distance_mask = bwdistsc1(single(obj.ch2.mask),[1 1 zscale],maxval);
        end
        
        function calculateDensity_abetaINch1(obj)
            nbins = 30;
            density = zeros(nbins,1);
            for nn = 1:nbins
                if nn == 1
                    curr_distmask = (obj.ch1.distance_mask<=nn).*obj.cellmask;
                else
                    curr_distmask = (obj.ch1.distance_mask>(nn-1) & obj.ch1.distance_mask<=nn).*obj.cellmask;
                end
                density(nn) = sum(obj.abeta.image.*curr_distmask)./sum(curr_distmask);
            end
            avgdensity = sum(obj.abeta.image.*obj.cellmask)./sum(obj.cellmask);
            xs = 1:nbins;
            figure; plot(xs,density/avgdensity)
        end
        function calculateDensity_abetaINch2(obj)
            nbins = 30; 
            density = zeros(nbins,1);
            for nn = 1:nbins
                curr_distmask = [];
                if nn == 1
                    curr_distmask = (obj.ch2.distance_mask<=nn).*obj.cellmask;
                else
                    curr_distmask = (obj.ch2.distance_mask>(nn-1) & obj.ch2.distance_mask<=nn).*obj.cellmask;
                end
                density(nn) = sum(obj.abeta.image.*curr_distmask)./sum(curr_distmask);
            end
            avgdensity = sum(obj.abeta.image.*obj.cellmask)./sum(obj.cellmask);
            xs = 1:nbins;
            figure; plot(xs,density/avgdensity)
        end
        
        
        
        %-- helpers for making masks. if you want to change how the masks
        %are made then change which static method these are calling
        function make_maskchAB(obj)
             obj.abeta.mask = SIM.make_maskABeta(obj.abeta.image);
        end
        function make_maskch1(obj)
            obj.ch1.mask = SIM.make_maskSynapseMarker(obj.ch1.image);
        end
        function make_maskch2(obj)
            obj.ch2.mask = SIM.make_maskSynapseMarker(obj.ch2.image);
        end
        %----
        function measure_AB(obj)
        end
        function measure_ch1(obj)
        end
        function measure_ch2(obj)
        end
        function measure_AB_inch1(obj)
        end
        function measure_AB_inch2(obj)
        end
    end
    
    methods (Static)
        function mask = make_maskABeta(image,params)
            ab = gaussf(image);
            out = GeneralAnalysis.imgLaplaceCutoff(ab,[2 2 1],[2 2 1]);
            thr = multithresh(single(out),2);
            mask_start = out>thr(1);
            sch2 = sum(dip_image(image),[],3);
            sch2_series = repmat(sch2,1,1,size(image,3));
            wshed = GeneralAnalysis.watershed_timeseries(-gaussf(sch2_series,1),1);
            mask = mask_start.*~wshed;
        end
        
        function mask = make_maskSynapseMarker(image,params)
            gch = gaussf(image);
            out = GeneralAnalysis.imgLaplaceCutoff(gch,[2 2 1],[2 2 1]);
            thr = multithresh(single(out),2);
            mask = bdilation(out>thr(1),1);
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
            
            zbool = cellfun(@(x) contains(x,'Z='),splstr1);
            zstr = splstr1{zbool};
            zparts = strsplit(zstr,{'=','/'});
            info.currZ = str2double(zparts{2});
            info.totalZ = str2double(zparts{3});
            
            colbool = cellfun(@(x) contains(x,'C='),splstr1);
            colstr = splstr1{colbool};
            colparts = strsplit(colstr,{'=','/'});
            info.currCol = str2double(colparts{2});
            info.totalCol = str2double(colparts{3});
        end
    end
end