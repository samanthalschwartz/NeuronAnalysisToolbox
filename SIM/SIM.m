classdef SIM < handle
    properties
        filepath = '';   
        savepath = [];
        channelordering = [2 1 3]; % array indicating which channel is which;
        %         channelordering(1) = abeta channel in image
        %         channelordering(2) = channel in image corresponding to 'ch1';
        %         channelordering(3) = channel in image corresponding to 'ch2';
        channelorderingstr = {'chABeta','ch1','ch2'};
        planeTOP = [];
        planeBOTTOM = [];% in dipimage format: 0 is first plane
        sizecutoff = 10; 
        abeta = struct('rawimage',[],'image',[],'mask',[],'distance_mask',[],'labeled_mask',[],'msr',[],'sizes',[],'densities',[],'COM_image',[]);
        ch1 = struct('rawimage',[],'image',[],'distance_mask',[],'mask',[],'mask_highsense',[],'name','',...
            'thisCh',struct('msr',[],'sizes',[],'densities',[]),'results',[],...
            'abetaCh',struct('msr',[],'sizes',[],'densities',[],'radialdensity_norm',[],'radialdensity_raw',[],...
                        'radialnumberdensity_norm',[],'radialnumberdensity_raw',[],...
                        'cumulative_radialnumberdensity_norm',[],'cumulative_radialnumberdnnnpppppensity_raw',[]));
        ch2 = struct('rawimage',[],'image',[],'distance_mask',[],'mask',[],'mask_highsense',[],'name','',...
            'thisCh',struct('msr',[],'sizes',[],'densities',[]),'results',[],...
            'abetaCh',struct('msr',[],'sizes',[],'densities',[],'radialdensity_norm',[],'radialdensity_raw',[],...
                        'radialnumberdensity_norm',[],'radialnumberdensity_raw',[],...
                        'cumulative_radialnumberdensity_norm',[],'cumulative_radialnumberdensity_raw',[]));
        measurements = {'size','sum','Gravity'};
        results = [];
        results_moreselective = [];
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
            if iscell(obj.filepath)
                filepath = obj.filepath{1};
            end
            try
            [image,obj.metadata] = obj.ndFileloader(obj.filepath); % load file
            catch
            disp('File did not load - make sure it is an ND file?');
            end
            if size(image,4)<3  % assign channels
                obj.abeta.rawimage = image(:,:,:,obj.channelordering(1));
                obj.ch1.rawimage = image(:,:,:,obj.channelordering(2));
            else
                obj.abeta.rawimage = image(:,:,:,obj.channelordering(1));
                obj.ch1.rawimage = image(:,:,:,obj.channelordering(2));
                obj.ch2.rawimage = image(:,:,:,obj.channelordering(3));
            end
            %-- get some meta data info
            obj.XYpxsize = double(obj.metadata.getPixelsPhysicalSizeX(0).value());           % returns value in default unit
            obj.XYpxsize_units = char(obj.metadata.getPixelsPhysicalSizeX(0).unit().getSymbol());
            obj.Zpxsize = double(obj.metadata.getPixelsPhysicalSizeZ(0).value());           % returns value in default unit
            obj.Zpxsize_units = char(obj.metadata.getPixelsPhysicalSizeZ(0).unit().getSymbol());
        end
        function resetCh2(obj)
           obj.ch2 = struct('rawimage',[],'image',[],'distance_mask',[],'mask',[],'mask_highsense',[],'name','',...
            'thisCh',struct('msr',[],'sizes',[],'densities',[]),'results',[],...
            'abetaCh',struct('msr',[],'sizes',[],'densities',[],'radialdensity_norm',[],'radialdensity_raw',[],...
                        'radialnumberdensity_norm',[],'radialnumberdensity_raw',[],...
                        'cumulative_radialnumberdensity_norm',[],'cumulative_radialnumberdensity_raw',[])); 
        end
        function resetCh1(obj)
            obj.ch1 = struct('rawimage',[],'image',[],'distance_mask',[],'mask',[],'mask_highsense',[],'name','',...
                'thisCh',struct('msr',[],'sizes',[],'densities',[]),'results',[],...
                'abetaCh',struct('msr',[],'sizes',[],'densities',[],'radialdensity_norm',[],'radialdensity_raw',[],...
                'radialnumberdensity_norm',[],'radialnumberdensity_raw',[],...
                'cumulative_radialnumberdensity_norm',[],'cumulative_radialnumberdnnnpppppensity_raw',[]));
        end
        function loadtiff(obj)
            if nargin>1
                obj.filepath = filepath;
            else
                if isempty(obj.filepath)
                    filepath = uipickfiles('Prompt','Pick a .nd File');
                    obj.filepath = filepath;
                end
            end
            if iscell(obj.filepath)
                filepath = obj.filepath{1};
            end
            try
            image = loadtiff(obj.filepath);
            catch
            disp('File did not load - make sure it is an ND file?');
            end
            if size(image,4)<3  % assign channels
                obj.abeta.rawimage = image(:,:,:,obj.channelordering(1));
                obj.ch1.rawimage = image(:,:,:,obj.channelordering(2));
            else
                obj.abeta.rawimage = image(:,:,:,obj.channelordering(1));
                obj.ch1.rawimage = image(:,:,:,obj.channelordering(2));
                obj.ch2.rawimage = image(:,:,:,obj.channelordering(3));
            end
            %-- get some meta data info
%             obj.XYpxsize = double(obj.metadata.getPixelsPhysicalSizeX(0).value());           % returns value in default unit
%             obj.XYpxsize_units = char(obj.metadata.getPixelsPhysicalSizeX(0).unit().getSymbol());
%             obj.Zpxsize = double(obj.metadata.getPixelsPhysicalSizeZ(0).value());           % returns value in default unit
%             obj.Zpxsize_units = char(obj.metadata.getPixelsPhysicalSizeZ(0).unit().getSymbol());            
        end
        function setimage(obj) 
            % sets obj.ch1.image etc to only include the selected frames (if obj.planeBOTTOM/obj.planeTOP are set)
            % obj.rawimage is unaltered
            % automatically orders top and bottom in case they were set wrong
            if ~isempty(obj.planeTOP) || ~isempty(obj.planeBOTTOM) % check that these values are set
                old_top  = str2double(obj.planeTOP);
                old_bottom = str2double(obj.planeBOTTOM);
                if ~isnan(old_bottom) % checks if bottom is a number
                    bottom = min(old_top,old_bottom);
                else
                    bottom = 0;
                end
                if ~isnan(old_top) % checks if bottom is a number
                    top = max(old_top,old_bottom);
                else
                    top = size(obj.ch1.image,3);
                end
                obj.ch1.image = obj.ch1.rawimage(:,:,bottom+1:top+1);
                obj.ch2.image = obj.ch2.rawimage(:,:,bottom+1:top+1);
                obj.abeta.image = obj.abeta.rawimage(:,:,bottom+1:top+1);
            else
                obj.ch1.image = obj.ch1.rawimage;
                obj.ch2.image = obj.ch2.rawimage;
                obj.abeta.image = obj.abeta.rawimage;
                
            end
        end
        function make_laplacemasks(obj)
           obj.ch1.mask = obj.laplacemask(obj.ch1.image);
           obj.ch2.mask = obj.laplacemask(obj.ch2.image);
           obj.abeta.mask = obj.laplacemask(obj.abeta.image);
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
        function simulationAbeta(obj,distance)
            if nargin<2
                distance = 20;
            end
            if isempty(obj.abeta.COM_image)
                obj.abetaCOM;
            end
            % channnel 1
            testim = 0.*obj.ch1.distance_mask;
            for zz = 1:size(obj.abeta.image,3)
                %
                goodmask = obj.ch1.distance_mask(:,:,zz)<distance;
                % find how many abeta COM in this frame
                currmask = obj.abeta.COM_image(:,:,zz-1).*goodmask;
                numabeta = sum(currmask(:));
                sz_x = size(goodmask,1);
                sz_y = size(goodmask,2);
                % simulation
                ids = find(goodmask);
                simframemask = goodmask.*0;
                randids = ids(1+round((numel(ids)-1).*rand(1,numabeta)));
                simframemask(randids) = 1;
                testim(:,:,zz) = simframemask;
            end
            obj.ch1.abetaSIM.COMimage = testim;
            obj.ch1.abetaSIM.image = bdilation(logical(testim),3);
            obj.ch1.abetaSIM.dist = distance;
            clear testim;
            % channnel 2
            testim = 0.*obj.ch2.distance_mask;
            for zz = 1:size(obj.abeta.image,3)
                %
                goodmask = obj.ch2.distance_mask(:,:,zz)<distance;
                % find how many abeta COM in this frame
                currmask = obj.abeta.COM_image(:,:,zz-1).*goodmask;
                numabeta = sum(currmask(:));
                % simulation
                ids = find(goodmask);
                simframemask = goodmask.*0;
                randids = ids(1+round((numel(ids)-1).*rand(1,numabeta)));
                simframemask(randids) = 1;
                testim(:,:,zz) = simframemask;
            end
            obj.ch2.abetaSIM.COMimage = testim;
            obj.ch2.abetaSIM.image = bdilation(logical(testim),3);
            obj.ch2.abetaSIM.dist = distance; 
        end
        
        function abetaDensityAlongPrePost(obj,ROIstring)
            switch ROIstring
                case 'selectedROIs'
                    regions = obj.results.selectedROIs;
                case 'moreselectiveROIs'
                    regions = obj.results.moreselectiveROIs;
            end
            totab = sum(obj.abeta.COM_image);
            plotvals = zeros(2,totab);
            prepostdis_list = zeros(1,numel(regions));
            numabeta_aroundsynapse = zeros(1,numel(regions));
            wb = waitbar(0);
            cnt = 0;
            prepostdis_cnt = 0;
            for aa = 1:numel(regions)
                prepostdis_cnt = prepostdis_cnt+1;
                vertices = regions{aa};
                if size(vertices,1)<3
                    disp('Not enough vertices to make a region');
                    continue;
                end
                zscale = obj.Zpxsize/obj.XYpxsize;
                coords_pre = getCOMfromVerts(vertices,obj.ch1.mask);
                coords_post = getCOMfromVerts(vertices,obj.ch2.mask);
                if isempty(coords_pre) || isempty(coords_post)
                    continue;
                end
                
                PrPo_vec = coords_post - coords_pre;
                PrPo_vec(3) = PrPo_vec(3)*zscale;
                Lmag = sqrt(sum(PrPo_vec).^2);
                prepostdis_list(prepostdis_cnt) = Lmag;
                
                lowbound = max(min(coords_pre,coords_post)-[20;20;3],[0; 0; 0]);
                maxbound = min(max(coords_pre,coords_post)+[20;20;3],(size(obj.abeta.image)-1)');
                
                abeta_im = obj.abeta.COM_image;
                labab = label(logical(abeta_im),1);
                cropped_labab = labab(floor(lowbound(1)):ceil(maxbound(1)),floor(lowbound(2)):ceil(maxbound(2)),floor(lowbound(3)):ceil(maxbound(3)));
                ids = unique(single(cropped_labab(:))); ids = ids(ids~=0);
                numabeta_aroundsynapse(prepostdis_cnt) = numel(ids);
                for ii =1 :numel(ids)
                    cnt= cnt+1;
                    c = findcoord(labab==ids(ii));
                    p1 = c';
                    
                    %     pre_img = 0*obj.ch1.mask;
                    %     pre_img( round(coords_pre(1)),round(coords_pre(2)),round(coords_pre(3)) ) = 1;
                    %     pre_imgsub = pre_img(floor(lowbound(1)):ceil(maxbound(1)),floor(lowbound(2)):ceil(maxbound(2)),floor(lowbound(3)):ceil(maxbound(3)))
                    %
                    %     post_img = 0*obj.ch1.mask;
                    %     post_img( round(coords_post(1)),round(coords_post(2)),round(coords_post(3)) ) = 1;
                    %     post_imgsub = post_img(floor(lowbound(1)):ceil(maxbound(1)),floor(lowbound(2)):ceil(maxbound(2)),floor(lowbound(3)):ceil(maxbound(3)))
                    L1 = coords_pre;
                    L2 = coords_post;
                    c1 = (L1+L2)/2;
                    L_vec = L2 - c1;
                    L_vec(3) = L_vec(3)*zscale;
                    pvec = p1 - c1;
                    pvec(3) = pvec(3)*zscale;
                    [Zdist, rdist] = vectordotproduct(L_vec,pvec);
                    plotvals(:,cnt) = [rdist;Zdist];
                end
                %     joinchannels('rgb',abetaCOM,pre_imgsub,post_imgsub)
                waitbar(aa/numel(regions),wb);
            end
            close(wb);
            switch ROIstring
                case 'selectedROIs'
                    obj.results.prepostdis_list = prepostdis_list;
                    obj.results.numabeta_aroundsynapse = numabeta_aroundsynapse;
                    obj.results.plotvals = plotvals;
                case 'moreselectiveROIs'
                    obj.results_moreselective.selectedROIs = obj.results.moreselectiveROIs;
                    obj.results_moreselective.prepostdis_list = prepostdis_list;
                    obj.results_moreselective.numabeta_aroundsynapse = numabeta_aroundsynapse;
                    obj.results_moreselective.plotvals = plotvals;
            end
            
        end
        function boolval = pre_postOverlap(obj,vertices,overlapcutoff)
            assert(isequal(size(obj.ch1.mask),size(obj.ch2.mask)));
            if nargin<3
                overlapcutoff = 20;
            end
            out = poly2mask(vertices(:,1)',vertices(:,2)',size(obj.ch1.mask,1),size(obj.ch1.mask,2));
            lbch1 = label(obj.ch1.mask,1);
            lbch2 = label(obj.ch2.mask,1);
            if size(vertices,1)<3
                disp('Not enough vertices to make a region');
                return;
            end
            if size(vertices,2) == 2
                mask_ch1 = lbch1.*out;
                mask_ch2 = lbch2.*out;
            else
                mask_ch1 = lbch1(:,:,vertices(1,3)).*out;
                mask_ch2 = lbch2(:,:,vertices(1,3)).*out;
            end
            overlap = mask_ch1.*mask_ch2;
            if sum(overlap(:))>overlapcutoff
                boolval = 0;
            else
                boolval = 1;
            end
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
                % only take abeta that is above average background
                m(ii+1) = single(sum(currframe,currcellmask,[1 2]))./sum(currcellmask);
                premask(:,:,ii) = premask(:,:,ii).*~(currframe<m(ii+1));
            end
            obj.abeta.mask = premask;
            obj.abeta.labeled_mask = label(obj.abeta.mask,1);
        end
        function make_maskchX(obj)
           inputim = obj.abeta.image;
            premask_1  = SIM.make_maskSynapseMarker(inputim);
            obj.abeta.mask = removeWeirdSIMartifact(obj,inputim,premask_1);
            % second set of masks
            premask_2  = SIM.make_maskSynapseMarkerI(inputim);
            obj.abeta.mask_highsense = removeWeirdSIMartifact(obj,inputim,premask_2);  
        end
        function make_maskch1(obj)
            inputim = obj.ch1.image;
            if strcmp(obj.channelorderingstr{2},'PSD95i')
                premask_1  = SIM.make_maskPSD95ib(obj.ch1.image);
            else
                 premask_1  = SIM.make_maskSynapseMarker(inputim);
            end
            obj.ch1.mask = removeWeirdSIMartifact(obj,inputim,premask_1);
            % second set of masks
            premask_2  = SIM.make_maskSynapseMarkerI(inputim);
            obj.ch1.mask_highsense = removeWeirdSIMartifact(obj,inputim,premask_2);
        end
        function make_maskch2(obj)
            inputim = obj.ch2.image;
            % first set of masks
            premask_1  = SIM.make_maskSynapseMarker(inputim);
            obj.ch2.mask = removeWeirdSIMartifact(obj,inputim,premask_1);
            % second set of masks
            premask_2  = SIM.make_maskSynapseMarkerI(inputim);
            obj.ch2.mask_highsense = removeWeirdSIMartifact(obj,inputim,premask_2);
        end
        function mask = removeWeirdSIMartifact(obj,image,premask)
            m = zeros(1,size(obj.cellmask,3));
            if ~isa(image,'dip_image')
            image = dip_image(image);
            end
            for ii = 0:(size(obj.cellmask,3)-1)
                currframe = image(:,:,ii);
                currcellmask = ~obj.cellmask(:,:,ii);
                currcellmask = bdilation(currcellmask,2);
                m(ii+1) = single(sum(currframe,currcellmask,[1 2]))./sum(currcellmask);
                premask(:,:,ii) = premask(:,:,ii).*~(currframe<m(ii+1));
            end
            mask = premask;
        end
        function viewpatches(obj)
            joinchannels('rgb',obj.abeta.mask,obj.ch1.mask,obj.ch2.mask)
            hold on;
            for ii =1:numel(obj.results_moreselective.selectedROIs)
                patch(obj.results_moreselective.selectedROIs{ii}(:,1)',obj.results_moreselective.selectedROIs{ii}(:,2),'white','EdgeColor','white','FaceColor','none','LineWidth',2)
            end
        end
        function selectPrePostROI(obj,alternate)  
            ch1 = obj.ch1.mask.*(obj.abeta.distance_mask<20);
            ch2 = obj.ch2.mask.*(obj.abeta.distance_mask<20);
            ov = joinchannels('rgb',ch1.*0,ch1,ch2);
            g = dipfig('ov');
            dipshow(ov);
            diptruesize(g,100);
            selectedROIs = [];
            patches = [];
            while(ishandle(g))
                try
                    w = waitforbuttonpress;
                    a = gcf;
                    if strcmp(a.CurrentCharacter,'t')
                        a = gcf;
                        try
                            vertices = dipdrawpolygon(a);
                            gcfinfo = get(g,'UserData');
                            currtime = gcfinfo.curslice;
                            vertices = cat(2,vertices,repmat(currtime,size(vertices,1),1));
                            selectedROIs = cat(1,selectedROIs,{vertices});
                            patch('Vertices',vertices,'EdgeColor',[1 0 0],'Faces',1:size(vertices,1),'FaceAlpha',0);
                            a.CurrentCharacter = 'f';
                        catch
                            break;
                        end
                       
                    end
                    %                     if strcmp(a.CurrentCharacter,'x')
%                         try
%                             patches
%                         catch
%                             break
%                         end
%                     end
                catch
                end
            end
            if nargin==2
                obj.results.moreselectiveROIs = selectedROIs;
            else
                obj.results.selectedROIs = selectedROIs;
            end
        end
        function abim = abetaCOM(obj)
            obj.measure_AB;
            alllabels = obj.abeta.msr.ID;
            abim = dip_image(0*obj.abeta.image);
            vals = round(obj.abeta.msr.Gravity);
            for nn = alllabels
                xval = vals(1,nn);
                yval = vals(2,nn);
                zval = vals(3,nn);
                abim(xval,yval,zval) = 1;
            end
            obj.abeta.COM_image = abim;
            %             for nn = alllabels
            %                 xval = round(obj.abeta.msr.Gravity(1,nn));
            %                 yval = round(obj.abeta.msr.Gravity(2,nn));
            %                 zval = round(obj.abeta.msr.Gravity(3,nn));
            %                 zmin = max(0,zval-1);
            %                 zmax = min(size(obj.abeta.mask,3)-1,zval+1);
            %                 abim(xval-1:xval+1,yval-1:yval+1,zmin:zmax) = 1;
            %             end
        end
        
         function calculateNumberDensityCOM(obj,reset,simbool)
            if nargin>1
                resetbool = reset;
                if nargin>2
                else
                    simbool = 0;
                end
            else
                resetbool = 0;
                simbool = 0;
            end
            
            % determine how many objects from Ch1 and Ch2 there are
            ch1_lb = label(obj.ch1.mask,1);
            ch1_numobj = max(ch1_lb(:));
            ch2_lb = label(obj.ch2.mask,1);
            ch2_numobj = max(ch2_lb(:));
%             ch1_maxdist = max(obj.ch1.distance_mask(:));
%             ch2_maxdist = max(obj.ch2.distance_mask(:));
            ch1_maxdist = obj.ch1.abetaSIM.dist;
            ch2_maxdist = obj.ch2.abetaSIM.dist;
            % make com mask from lb
            if ~isfield(obj.abeta,'COM_image') || isempty(obj.abeta.COM_image) || resetbool
                obj.abetaCOM();
            end
            if ~simbool
               abetCOMim_ch1 = double(obj.abeta.COM_image); 
               abetCOMim_ch2 = double(obj.abeta.COM_image);
            else
                if ~isfield(obj.ch1,'abetaSIM')
                    obj.simulationAbeta(20);
                end
               abetCOMim_ch1 = obj.ch1.abetaSIM.COMimage; 
               abetCOMim_ch2 = obj.ch2.abetaSIM.COMimage; 
            end
                
            bins = 0:(max(ch1_maxdist,ch2_maxdist));
            nbins = size(bins,2);
            ch1_numabeta = zeros(nbins,1);
            ch2_numabeta = zeros(nbins,1);
            ch1_volume= zeros(nbins,1);
            ch2_volume = zeros(nbins,1);
            wb = waitbar(0,'Calculating Number Density Information...');
            for nn = 1:nbins
                if nn == 1
                    ch1_curr_distmask = (obj.ch1.distance_mask<=bins(nn));
                    ch2_curr_distmask = (obj.ch2.distance_mask<=bins(nn));
                else
                    ch1_curr_distmask = (obj.ch1.distance_mask>bins(nn-1) & obj.ch1.distance_mask<=bins(nn));
                    ch2_curr_distmask = (obj.ch2.distance_mask>bins(nn-1) & obj.ch2.distance_mask<=bins(nn));
                end
                ch1_inmasklb = single(abetCOMim_ch1.*ch1_curr_distmask);
                ch1_numabeta(nn) = sum(ch1_inmasklb(:));
                ch1_volume(nn) = sum(ch1_curr_distmask(:));
                
                ch2_inmasklb = single(abetCOMim_ch2.*ch2_curr_distmask);
                ch2_numabeta(nn) = sum(ch2_inmasklb(:));  
                ch2_volume(nn) = sum(ch2_curr_distmask(:));
                try
                    waitbar(nn/nbins,wb);
                catch
                    wb = waitbar(nn/nbins,'Calculating Number Density Information...');
                end
            end
            close(wb);
            if ~simbool
                obj.ch1.results.bins = bins;
                obj.ch1.results.numabeta = ch1_numabeta;
                obj.ch1.results.numobj = ch1_numobj;
                obj.ch1.results.volume = ch1_volume;
                
                obj.ch2.results.bins = bins;
                obj.ch2.results.numabeta = ch2_numabeta;
                obj.ch2.results.numobj = ch2_numobj;
                obj.ch2.results.volume = ch2_volume;
            else
                obj.ch1.abetaSIM.results.bins = bins;
                obj.ch1.abetaSIM.results.numabeta = ch1_numabeta;
                obj.ch1.abetaSIM.results.numobj = ch1_numobj;
                obj.ch1.abetaSIM.results.volume = ch1_volume;
                
                obj.ch2.abetaSIM.results.bins = bins;
                obj.ch2.abetaSIM.results.numabeta = ch2_numabeta;
                obj.ch2.abetaSIM.results.numobj = ch2_numobj;
                obj.ch2.abetaSIM.results.volume = ch2_volume;
            end
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
                if isempty(obj.savepath)
                    obj.savepath = [obj.filepath(1:end-4) '_SIM.mat'];
                end
            else
                if isfolder(inputsavedir)
                    [~,NAME,~] = fileparts(obj.filepath);
                    obj.savepath = fullfile(inputsavedir,[NAME '_SIM.mat']);
                else
                    obj.savepath = inputsavedir;
                end
            end
            save(obj.savepath,'obj','-v7.3','-nocompression');
        end
        
        function [tempvals,msr] = mindistCh1tofirstAbeta(obj,redodistancebool)
            if isempty(obj.abeta.distance_mask) || nargin>1
                obj.make_distancemasks;
            end
            img_dist2ch1 = obj.abeta.distance_mask .* obj.ch1.mask;
%             img_dist2ch1(img_dist2ch1==0) = Inf;
            msr = measure(obj.ch1.mask,img_dist2ch1,{'MinVal'});      
            tempvals = msr.MinVal;
            tempvals(tempvals >= 40) = [];
            obj.results_closest_ch1_2ab = tempvals;
        end
        
        function [tempvals,msr] = mindistCh2tofirstAbeta(obj,redodistancebool)
            if isempty(obj.abeta.distance_mask) || nargin>1
                obj.make_distancemasks;
            end
            img_dist2ch2 = obj.abeta.distance_mask .* obj.ch2.mask;
%             img_dist2ch1(img_dist2ch1==0) = Inf;
            msr = measure(obj.ch2.mask,img_dist2ch2,{'MinVal'});      
            tempvals = msr.MinVal;
            tempvals(tempvals >= 40) = [];
            obj.results_closest_ch2_2ab = tempvals;
        end
        
        function calculate_minDist2Abeta(obj)
            [ch1_edges, ch1_N] = obj.static_calculate_minDist2Abeta(obj.ch1.mask,obj.abeta.distance_mask);
            [ch2_edges, ch2_N] = obj.static_calculate_minDist2Abeta(obj.ch2.mask,obj.abeta.distance_mask);
            obj.ch1.results.mindist2abeta.N = ch1_N;
            obj.ch1.results.mindist2abeta.edges = ch1_edges;
            obj.ch2.results.mindist2abeta.N = ch2_N;
            obj.ch2.results.mindist2abeta.edges = ch2_edges;
        end
        
        function calculate_SimulatedminDist2Abeta(obj)
            [ch1_edges, ch1_N] = obj.static_calculate_minDist2Abeta(obj.ch1.mask,obj.ch1.abetaSIM.distance_mask);
            [ch2_edges, ch2_N] = obj.static_calculate_minDist2Abeta(obj.ch2.mask,obj.ch2.abetaSIM.distance_mask);
            obj.ch1.abetaSIM.results.mindist2abeta.N = ch1_N;
            obj.ch1.abetaSIM.results.mindist2abeta.edges = ch1_edges;
            obj.ch2.abetaSIM.results.mindist2abeta.N = ch2_N;
            obj.ch2.abetaSIM.results.mindist2abeta.edges = ch2_edges;
        end
        
    end
    
    methods(Static)
        function mask = make_maskABeta(image)
            ab = gaussf(image);
            %             out = GeneralAnalysis.imgLaplaceCutoff(ab,[1 1 1],[1 1 1]);
            out = GeneralAnalysis.imgLaplaceCutoff(ab);
            %            out = gaussf(out,[2 2 1]);
            thr = multithresh(single(out),2);
            mask_start = out>thr(1);
            sch2 = sum(dip_image(image),[],3);
            sch2_series = repmat(sch2,1,1,size(image,3));
            wshed = GeneralAnalysis.watershed_timeseries(-gaussf(sch2_series,1),1);
            pre_mask = mask_start.*~wshed;
            lb = label(pre_mask,1,100,10^10);
            mask = lb>0;
        end
        
        function mask = make_maskSynapseMarker(image,params)
            gch = gaussf(image);
            out = GeneralAnalysis.imgLaplaceCutoff(gch,[2 2 1],[2 2 1]);
            thr = multithresh(single(out),2);
            mask = out>thr(2);
%             mask = bdilation(out>thr(2),1);
        end
        function mask = make_maskSynapseMarkerI(image,params)
            gch = gaussf(image);
            out = GeneralAnalysis.imgLaplaceCutoff(gch);
            thr = multithresh(single(out),2);
            mask = out>thr(1);
%             mask = bdilation(out>thr(2),1);
        end
        function mask = make_maskSynapseMarkerLow(image,params)
            gch = gaussf(image);
            out = GeneralAnalysis.imgLaplaceCutoff(gch,[2 2 1],[2 2 1]);
            thr = multithresh(single(out),1);
            mask = out>thr(1);
%             mask = bdilation(out>thr(1),1);
        end
        
        function mask = make_maskPSD95ib(image,params)
            gch = gaussf(image);
            out = GeneralAnalysis.imgLaplaceCutoff(gch,[2 2 1],[2 2 1]);
            out = gaussf(out,[2 2 1]);
            thr = multithresh(single(out),2);
            mask = bdilation(out>thr(2),1); 
        end
        
        function mask = laplacemask(image,params)
            if nargin<2 || isempty(params)
                params.minsize = 50;
                params.maxsize = 30000;
            end
            lp = laplace(image,2);
            img_laplcutoff = -lp;
            img_laplcutoff(img_laplcutoff<0) = 0;
            threshval = multithresh(single(img_laplcutoff),2);
            firstmask = img_laplcutoff>threshval(2);
            lb = label(firstmask,1,params.minsize,params.maxsize);
            mask = lb>0;
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
%         function results = calculateNumberDensity(mask1_distance,lb,cellmask,bins)
%             % mask1_distance mask - is distance transform of mask (ie: synapse marker mask)
%             % image is the raw image that the intensity will be calculated from (ie: abeta image)
%             % cellmask is the whole cell mask that the area will be taken and normalized from
%             % radial_numberdensity is: (#abeta puncta/area)/#gephyrin
%             
%             nbins = numel(bins);
%             %-- alocate empty arrays for values calculated
%             numbers = zeros(nbins,1);
%             volume = zeros(nbins,1);
%             %-- loop through and calculate some things
%             wb = waitbar(0,'Calculating Number Density Information...');
%             mask10s = mask1_distance==0;
%             mask1_lb = label(mask10s,1);
%             nummask1 = size(unique(mask1_lb));
%             for nn = 1:nbins
%                 if nn == 1
%                     curr_distmask = (mask1_distance<=bins(nn)).*cellmask;
%                 else
%                     curr_distmask = (mask1_distance>bins(nn-1) & mask1_distance<=bins(nn)).*cellmask;
%                 end
%                 inmasklb = single(lb.*curr_distmask);
%                 numbers(nn) = size(unique(inmasklb),1);
%                 volume(nn) = sum(curr_distmask); 
%                 try
%                 waitbar(nn/nbins,wb);
%                 catch
%                    wb = waitbar(nn/nbins,'Calculating Number Density Information...');
%                 end
%             end
%             close(wb);
%             results.d = bins;
%             results.radialnumber = numbers;
%             results.cumulativeradialnumber = cumsum(numbers);
%             results.volume = volume;
%             results.cumulativevolume = cumsum(volume);
%             results.totalnumber = size(unique(lb.*cellmask));
%             results.totalvolume = sum(cellmask);
%             results.radial_numberdensity = results.radialnumber;
%             results.cumulative_number_density = results.cumulativeradialnumber;
%             results.nummask1 = nummask1;
%         end 
             

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
        
        function [edges, N] = static_calculate_minDist2Abeta(ch_mask,abeta_mask)
            ch_lab = label(ch_mask,1);
            ch_msr = measure(ch_lab,abeta_mask,{'MinVal'}); 
            [N,edges] = histcounts(ch_msr.MinVal,'BinWidth',1,'Normalization','cdf');
        end
        function distance_mask = static_make_distancemasks(mask,zscale,maxval)
            if nargin<3
                maxval = 40;
                if nargin<2
                    zscale = 1;
                end
            end
            disp('Calculating Distance Mask');
            distance_mask = bwdistsc1(single(mask),[1 1 zscale],maxval);
        end
        function pts = BresenhamPoints(p1,p2)  
            d = p2-p1;
            N = max(abs(d));
            s = d/N;
            pts = zeros(N,3);
            pts(1,:) = p1;         
            for ii=1:N
                pts(ii+1,:) = pts(ii,:)+s;
            end
            pts = round(pts);
        end
    end
end