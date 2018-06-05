% datadir = '/nfs/olah/home/mjo/LidkeLab/Data/Syk Data/For Lifetimes';
% good filename examples:
% cell1_100ng#0001-2014-12-20-14-59-25 or KOWTcell5_preonly#0002-2016-7-8-11-47-1
% ---- some example calls:
% datdir = 'Z:\Sam\Syk\For Lifetimes';
% CellTypes = {'Peter161011KO','Peter161011Parental','Peter161011KO10%LP','Peter161011Bilayers'};
% SykTypes = {'WT','WTdas','YE','YF'};
% testdat = LifetimeAnalysis(datdir,CellTypes,SykTypes);
% allfilts = createFilters(testdat);
% missingIgE = testdat.linkIgEfilenames(filter);
classdef LifetimeAnalysis < handle
    properties (Constant=true, Hidden=true)
%         DefaultCellTypes = {'KO','Parental'};
%         DefaultSykTypes = {'WT','YE','YF'};
        SykRPTPattern = '*SYK.rpt';
        IGERPTPattern = '*IGE.rpt';
        
    end

    properties
        CellTypes;
        SykTypes;
        dataDir; % The directory where the data is
        obs;          % date
                      % dateStr
                      % dateTime
                      % cellType
                      % sykType
                      % conc
                      % preDate
                      % preTracks 
                      % preRPT
                      % seqTracks [#1, #2, #3, ...]                      
                      % seqRPT [#1, #2, ...]
    end

    properties 
        NTotalObs; %Total number of observed cells
        frameT; %Frame time in seconds
    end

    methods 
        function obj = LifetimeAnalysis(dataDir,cellTypes, sykTypes)
            if nargin>0
                obj.dataDir = dataDir;
                obj.CellTypes = cellTypes;
                obj.SykTypes = sykTypes;
                obj.reloadTracks();
            end
        end

        function gmm = makeGMM(obj, filters, seqs, offset, max_intensity)
            if nargin<4
                offset = 2;
                max_intensity=[];
            end
            if nargin<5
                max_intensity=[];
            end
            groupTs = obj.getTracks(filters, seqs);
            if ~isempty(max_intensity)
                groupTs = cellmap(@(Ts) Ts(obj.trackFilterIntensity(Ts, max_intensity)), groupTs);
            end
            Ds = cellmap(@(Ts) obj.trackDurationFrames(Ts), groupTs); % removed offset subtraction because GMM now incorporates this offset properly!
            Ds = cellmap(@(d) d(d>=1), Ds);
            gmm = GeometricMixtureModelMLE(Ds);
            gmm.offset = offset;
            gmm.setData(gmm.D);    
        end

        function Is = getallIntensities_noCPT(obj, filter, seqs)
            Is = cellmap(@(Ts) cellcat(obj.trackIntensities(Ts)), makecell(obj.getAllTracks_noCPT(filter,seqs)));
        end
        
        function Is = getIntensities(obj, filter, seqs)
            Is = cellmap(@(Ts) cellcat(obj.trackIntensities(Ts)), makecell(obj.getTracks(filter,seqs)));
        end

        function Is = getMeanIntensities(obj, filter, seqs)
            Is = cellmap(@(Ts) obj.trackMeanIntensities(Ts), makecell(obj.getTracks(filter,seqs)));
        end

        function Ds = getDurationFrames(obj, filter, seqs, Ithreshold)
            if nargin==4
                Ds = cellmap(@(Ts) obj.trackDurationFrames(Ts(obj.trackFilterIntensity(Ts,Ithreshold))), obj.getTracks(filter,seqs));
            else
                Ds = cellmap(@(Ts) obj.trackDurationFrames(Ts), obj.getTracks(filter,seqs));
            end
        end

        function [Ds, cellnames] = getDurationFramesByCell(obj, filter, seqs, max_intensity)
            assert(islogical(filter)); % Filter should be just one list of cells 
            Ncells = sum(filter);
            Ds = cell(1,Ncells);
            cellnames = cell(1,Ncells);
            filteredobs = obj.obs(filter);
            for cc = 1:Ncells
                cellnames{cc} = [filteredobs(cc).cellName '-' filteredobs(cc).v];
            end
            cellidxs = find(filter);
            for n=1:Ncells
                cell_filter = false(1,obj.NTotalObs);
                cell_filter(cellidxs(n)) = true;
                Ts = obj.getTracks(cell_filter,seqs);
                if iscell(Ts) && numel(Ts) == 1
                Ts = Ts{1};% not sure that this was right but changed to call Ts{1} so input not a 1x1 cell
                end
                if isempty(Ts)
                    Ds{n} = [];
                else
                    if ~isempty(max_intensity)
                        threshTsId = obj.trackFilterIntensity(Ts, max_intensity);
%                         Ds{n} = obj.trackDurationFrames(Ts(threshTsId)) - offset;
                         Ds{n} = obj.trackDurationFrames(Ts(threshTsId));
                    else
                        Ds{n} = obj.trackDurationFrames(Ts);
                    end
                end
            end
        end
        
        function savename = plotCPD70ValperCell(obj, allfilts, allseqs,max_intensity,offset,savedir)
            if nargin<6
                savedir = pwd;
            end
            fffs  = fieldnames(allfilts);
            savename = cell(1,numel(fffs));
            for ff = 1:numel(fffs)
                filt = getfield(allfilts,fffs{ff});
                filtername = strrep(fffs{ff},'_','-');
                Ds = {};
                for ii = allseqs
                    seqs = allseqs(ii+1);
                    [Ds(ii+1,:), cellnames] = obj.getDurationFramesByCell(filt, seqs, max_intensity);
                end
                % now loop through and plot?
                output = zeros(size(Ds));
                for dd = 1:size(Ds,1)
                    for rr = 1:size(Ds,2)
                        ts = Ds{dd,rr};
                        tlens = ts(ts>offset);
                        if ~isempty(tlens)
                            [F,X] = ecdf(tlens);
                            tVal = interp1(F,X,0.7);
                            output(dd,rr) = tVal;
                        else
                            output(dd,rr) = NaN;
                        end
                    end
                end
                f = figure;plot(output,'XData',0:8,'LineWidth',2); ylim([0 30]);
                ylabel('70% CDF Tracklength (seconds)');
                xlabel('Movie #');
                title(['Tracklength Summary for: ' filtername]);
                legend(cellnames,'Location','northwest')
                savename{ff} = ['timeCourse_' fffs{ff} '_photonThresh=' num2str(max_intensity) 'offset=' num2str(offset)];
                print(fullfile(savedir,savename{ff}),'-dpng','-r0');
                close(f);
            end
        end
        
        function f = plotDurationFramesHist(obj, filter, seqs, names, frameTime, Ithreshold)
            f = figure('Position',[100, 100, 800,600]);
%             whitebg(f);
%             axH = axes();
            hold('on');
            if nargin>5
                 Ds = obj.getDurationFrames(filter, seqs, Ithreshold);
            else
                Ds = obj.getDurationFrames(filter, seqs);
            end
            cs = lines(10);
            sykTypes = unique({obj.obs(:).sykType});
            concs = [1,10,100,1000];
            if nargin>4
                bins = logspace(log10(frameTime),log10(1000*frameTime),30);
                bins = linspace(3*frameTime,1000*frameTime,50);
            else
                frameTime =1;
                bins = logspace(0,3,30);
                bins = linspace(3*frameTime,1000*frameTime,50);
            end
            xlim([3,1000].*frameTime);
            for n=flip(1:numel(Ds))
%                 temphist = histogram(Ds{n}*frameTime, bins(2:end),'Normalization','PDF');
                temphist = histogram(Ds{n}*frameTime,'Normalization','PDF');
%                 C = C ./sum(C);
%                 [counts, bins] = histcounts(Ds{n}.*frameTime,bins,'Normalization','probability');
%                 bin_centers = (bins(2:end) + bins(1:end-1))/2;
%                 H = bar(bins(2:end-1),C); 
                temphist.DisplayName = names{n};
%                 H.LineWidth = 1.5;
%                 if strfind(names{n},'pre')
%                     H.LineStyle='--';
%                 end
%                 O = obj.obs(find(filter{n},1,'first'));
%                 H.Color = cs{find(O.conc==concs ,1,'first')};
%                 H.Color = cs(find(strcmp(sykTypes,O.sykType) ,1,'first'),:);        
            end
            legend('location','best');
            ylabel('Normalized Probability');
            if nargin>4                
                xlabel('Duration (Secs)');
            else
                xlabel('Duration (frames)');
            end
            axH.XScale='log';
            axH.YScale='linear';
        end
        function f = plotDurationFramesHist_nonNorm(obj, filter, seqs, names, frameTime, Ithreshold)
            f = figure('Position',[100, 100, 800,600]);
            %             whitebg(f);
            axH = axes();
            hold('on');
            if nargin>5
                Ds = obj.getDurationFrames(filter, seqs, Ithreshold);
            else
                Ds = obj.getDurationFrames(filter, seqs);
            end
            cs = lines(10);
            sykTypes = unique({obj.obs(:).sykType});
            concs = [1,10,100,1000];
            if nargin>4
                bins = logspace(log10(frameTime),log10(1000*frameTime),30);
%                 bins = linspace(3*frameTime,1000*frameTime,50);
            else
                frameTime =1;
                bins = logspace(0,3,30);
%                 bins = linspace(3*frameTime,1000*frameTime,50);
            end
            xlim([3,1000].*frameTime);
            for n=1:numel(Ds)
                [H] = histogram(Ds{n}*frameTime,bins,'FaceColor','none','DisplayName', names{n}); %--- added in. uncomment to use equal bins
%                 H.DisplayName = names{n};
                C = H.Values./sum(H.Values);
                f = plot(C,bins(2:end));
%                 f.DisplayName = names{n};
%                 [C,~] = histc(Ds{n}*frameTime, bins(2:end));
                %                 [counts, bins] = histcounts(Ds{n}.*frameTime,bins,'Normalization','probability');
                %                 bin_centers = (bins(2:end) + bins(1:end-1))/2;
%                 H = plot(bins(2:end),C.Values);
                
%                 H.LineWidth = 1.5;
                %                 if strfind(names{n},'pre')
                %                     H.LineStyle='--';
                %                 end
                O = obj.obs(find(filter{n},1,'first'));
                %                 H.Color = cs{find(O.conc==concs ,1,'first')};
                %                 H.Color = cs(find(strcmp(sykTypes,O.sykType) ,1,'first'),:);
            end
            legend('location','best');
            ylabel('# Tracks');
            if nargin>4
                xlabel('Duration (Secs)');
            else
                xlabel('Duration (frames)');
            end
            axH.XScale='log';
            axH.YScale='linear';
        end
        function f = plotDurationFramesCDF(obj, filter, seqs, names, frameTime, Ithreshold,offset)
            
            f = figure('Position',[100, 100, 800,600]);
%             whitebg(f);
            axH = axes();
            hold('on');
            if nargin>5
                Ds = obj.getDurationFrames(filter, seqs, Ithreshold);
            else
                Ds = obj.getDurationFrames(filter, seqs);
            end
            if nargin<7
                offset = 0;
            end
            
%             cs = lines(20);
            cs = distinguishable_colors(20);
            sykTypes = unique({obj.obs(:).sykType});
            if nargin<=4
                frameTime =1;
            end
            xlim([0,1000].*frameTime);
            
            for n=1:numel(Ds)
                Ds{n}(Ds{n} < offset) = [];
                [x,f] = ecdf(Ds{n}.*frameTime);
                numtracks = numel(Ds{n});
                H = plot(f(2:end),x(2:end));
                H.DisplayName = [names{n} ' # Tracks = ' num2str(numtracks)];
                H.LineWidth = 1.5;
                H.Color = cs(n,:);
                if strfind(names{n},'pre')
                    H.LineStyle='--';
                end
                % % uncomment/comment to make all syk types the same color
%                 O = obj.obs(find(filter{n},1,'first'));
%                 H.Color = cs(find(strcmp(sykTypes,O.sykType) ,1,'first'),:);
            end
            ylabel('Cumulative Probability');
            if nargin>4                
                xlabel('Duration (Secs)');
            else
                xlabel('Duration (frames)');
            end
            legend('location','best');
            axH.XScale='log';
            axH.YScale='linear';
        end
        
        function fig = plotDurationFramesCDF_andFit(obj, filter, seqs, names, fitinfo, frameTime, Ithreshold)
            % fitinfo is a cell array - size of filter
            % fitinfo{1}.selected_model and fitinfo{1}.theta_mle
            fig = figure('Position',[100, 100, 800,600]);
%             whitebg(f);
            axH = axes();
            hold('on');
            offset = 2;
            if nargin>5
                pre = obj.getDurationFrames(filter, seqs, Ithreshold);
                Ds = cellmap(@(x) x - offset,pre);
            else
                pre = obj.getDurationFrames(filter, seqs);
                Ds = cellmap(@(x) x - offset,pre);
            end
            cs = lines(10);
            sykTypes = unique({obj.obs(:).sykType});
            if nargin<=4
                frameTime =1;
            end
            xlim([0,1000].*frameTime);
            for n=1:numel(Ds)
                [x,f] = ecdf(Ds{n}.*frameTime);
                H = plot(f(2:end),x(2:end));
                H.DisplayName = names{n};
                H.LineWidth = 1.5;
                if strfind(names{n},'pre')
                    H.LineStyle='--';
                end
                O = obj.obs(find(filter{n},1,'first'));
                H.Color = cs(find(strcmp(sykTypes,O.sykType) ,1,'first'),:);
                
                % now check the gmm for this
                sim_gmm = GeometricMixtureModelMLE();
                sim_gmm.K = 1;
                sim_gmm.Ntot = 50000;
                switch fitinfo{n}.selected_model
                    case '2comp_ind_alpha'
%                         sim_gmm.simulate_2comp_ind_alpha(numel(Ds{n}), fitinfo{n}.theta_mle);
                        sim_gmm.simulate_2comp_ind_alpha(sim_gmm.Ntot, fitinfo{n}.theta_mle);
                    case '1comp'
%                         sim_gmm.simulate_1comp(numel(Ds{n}), fitinfo{n}.theta_mle);
                        sim_gmm.simulate_1comp(sim_gmm.Ntot, fitinfo{n}.theta_mle)
                    case '2comp_ind'
%                         sim_gmm.simulate_2comp_ind(numel(Ds{n}), fitinfo{n}.theta_mle);
                        sim_gmm.simulate_2comp_ind(sim_gmm.Ntot, fitinfo{n}.theta_mle)
                end
                simdata = (sim_gmm.D{1}(sim_gmm.D{1}>offset)-offset)*frameTime; %hard coded 2 here - matches with tracking data cutoff
                [y,g] = ecdf(simdata);
                I = plot(g(2:end),y(2:end));
                I.DisplayName = [names{n} ' Sim'];
                I.LineWidth = 1.5;
                I.Color = [0 0 0];
                I.LineStyle = '--';
                
                
            end
            ylabel('Cumulative Probability');
            if nargin>4                
                xlabel('Duration (Secs)');
            else
                xlabel('Duration (frames)');
            end
            legend('location','best');
            axH.XScale='log';
            axH.YScale='linear';
     %----- NEED TO FIX THIS!!!% ----------       
%             %-- plot residuals (assuming only 1 conditions
%             % x is actual data, y is sim data
%             if numel(Ds)==1
%                 f = figure('Position',[100, 100, 800,600]);
%                 datsize = size(y,1);
%                 R = plot(g(2:end),x(2:datsize)-y(2:end));
%                 ylabel('Absolute Difference (Data - Sim)');
%             if nargin>4                
%                 xlabel('Duration (Secs)');
%             else
%                 xlabel('Duration (frames)');
%             end
%             title('Residual Plot');
%             end
        end
        
        function f = plotDurationFramesByCellCDF(obj, filter, seqs, plotPre, maxIntensity)
            
            if nargin<3
                plotPre = false;
            end
            if nargin<4
                maxIntensity = Inf;
            end
            f = figure('Position',[100, 100, 800,600]);
            whitebg(f);
            axH = axes();
            concs = unique([obj.obs(:).conc]);
            sykTypes = unique({obj.obs(:).sykType});
            cs = {[1,0,0],[0,1,0],[0,0.2,1],[1,0,1],[0,1,1],[1,0.5,0],[0.5,0,1]};
            hold('on');
            cellIdxs = find(filter);
            
            Ds = obj.getDurationFramesByCell(filter, seqs,maxIntensity);
            colmap = hsv(numel(Ds));
            for n=1:numel(Ds)
                if isempty(Ds{n}); continue; end
                [x,f] = ecdf(Ds{n});
                O = obj.obs(cellIdxs(n));  
                %
                H = plot(f,x);
                H.DisplayName = sprintf('%s-%s',O.cellName, O.dateStr);
                H.LineWidth = 1;
                H.Color = colmap(n,:);
%                 H.Color = cs{find(O.conc==concs ,1,'first')};
%                 H.Color = cs{find(strcmp(sykTypes,O.sykType) ,1,'first')};
            end
            if plotPre
                Ds = obj.getDurationFramesByCell(filter, 0,maxIntensity);
                for n=1:numel(Ds)
                    if isempty(Ds{n}); continue; end
                    [x,f] = ecdf(Ds{n});
                    H = plot(f,x);
                    O = obj.obs(cellIdxs(n));
                    H.DisplayName = sprintf('%s-%ing-pre%',O.cellName, O.conc);
                    H.LineWidth = 1;
                    H.LineStyle = '--';
%                     H.Color = cs{find(O.conc==concs ,1,'first')};
                    H.Color = cs{find(strcmp(O.sykType,sykTypes) ,1,'first')};
                end
            end
            ylabel('Cumulative Probability');
            xlabel('Duration (frames)');
            legend('location','best');
            axH.XScale='log';
            axH.YScale='linear';
        end

        function f = plotIntensitiesHist(obj, filter, seqs, names)
            f = figure('Position',[100, 100, 800,600]);
            axH = axes();
            Is = obj.getIntensities(filter, seqs);
            BaseRPT.plot_logDistributions(axH, Is, names, true);
            xlabel('Intensity (photons)');
            axH.XScale='log';
        end

        function f = plotMeanIntensitiesHist(obj, filter, seqs, names)
            f = figure('Position',[100, 100, 800,600]);
%             whitebg(f);
            axH = axes();
            Is = obj.getMeanIntensities(filter, seqs);
            BaseRPT.plot_distributions(axH, Is, names, true);
            xlabel('Intensity (photons)');
            axH.XScale='linear';
        end
 
        function f = plotIntensityLengthScatter(obj, filter, seqs, names)
            filter = makecell(filter);
            seqs = makecell(seqs);
            names = makecell(names);
            N = numel(filter);
            Ts = cell(N,1);
            Is = cell(N,1);
            Ls = cell(N,1);
            for n = 1:N
                Ts{n} = cellcat(obj.getTracks(filter{n}, seqs{n}));
                Is{n} = obj.trackMaxIntensities(Ts{n});
                Ls{n} = obj.trackDurationFrames(Ts{n}).*obj.frameT;
            end
            mxL = max(cellfun(@max,Ls));
            mxI = max(cellfun(@max,Is));

%             mxL = 200.*obj.frameT;
%             mxI = 600;
            f = figure('Position',[100, 100, 800,600]);
            whitebg(f);
            for n = 1:N
                subplot(ceil(N/2), 2, n);
%                 plot( Is, Ls, 'LineStyle','none','DisplayName', names{n},...
%                               'Marker','o','MarkerSize',3,...
%                               'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,1]);
%                 Is_edge = logspace(log10(10),log10(600),60);
%                 Ls_edge = logspace(log10(3),log10(200),60);
%                 [C,G] = hist3([Is{n}',Ls{n}'], 'Edges',{Ls_edge,Is_edge});
%                 surface(G{1}, G{2}, C,'EdgeColor','none','DisplayName', names{n});
                 [C,G] = hist3([Is{n}',Ls{n}'],[30, 30]);
%                 surface(G{1}, G{2}, C,'EdgeColor','none','DisplayName', names{n});
                contour(G{1}, G{2}, C,'DisplayName', names{n},'Fill','on');
                ylim([min(G{2}), max(G{2})]);
                xlim([min(G{1}), max(G{1})]);
                ylabel('Track Max Intensity (photons)');
                xlabel('Track Duration (s)');
                legend('location','best');
            end
        end

        function Ts = getTracks(obj, filter, seqs)
            % seqs = [0 1 2 3] to get pre, #0001 #0002 and #0003 tracks.
            if ~iscell(filter)
                filter = {filter};
            end
            assert(all(cellfun(@numel,filter) == numel(obj.obs)));  % filter is the wrong size - is it for the correct database?
            function Ts = selectTracks(Obs, seqIdx) % local function to help out
                if seqIdx==0
                    if isfield(Obs,'changePointIDpreTracks')
                        Ts = Obs.preTracks(logical(Obs.changePointIDpreTracks));
                    else
                        Ts = Obs.preTracks;
                    end
                elseif seqIdx <= numel(Obs.seqTracks)
                    if isfield(Obs,'changePointIDseqTracks')
                        Ts = Obs.seqTracks{seqIdx}(logical(Obs.changePointIDseqTracks{seqIdx}));
                    else
                        Ts = Obs.seqTracks{seqIdx};
                    end
                else
                    Ts = [];
                end
            end
            Ts = cellmap(@(F,S) ... % for each filter and seqs list
                cellcat(cellmap(@(O) ... % for each observations
                cellcat(cellmap(@(i) selectTracks(O,i),S)), ... % Collect all tracks from those sequences
                obj.obs(F))), makecell(filter),makecell(seqs));
            
            %             if numel(Ts)==1
            %                 Ts = Ts{1};
            %             end
        end
        function Ts = getAllTracks_noCPT(obj, filter, seqs)
            % this function gets all the tracks before the change point
            % thresholding. Mainly used for supplemental figure. Same as
            % the way the function used to be before changePointIDpreTracks
            % was added to getTracks function
            % seqs = [0 1 2 3] to get pre, #0001 #0002 and #0003 tracks.
            if ~iscell(filter)
                filter = {filter};
            end
            assert(all(cellfun(@numel,filter) == numel(obj.obs)));  % filter is the wrong size - is it for the correct database?
            function Ts = selectTracks(Obs, seqIdx) % local function to help out
                if seqIdx==0
                        Ts = Obs.preTracks;
                elseif seqIdx <= numel(Obs.seqTracks)
                    
                        Ts = Obs.seqTracks{seqIdx};
                else
                    Ts = [];
                end
            end
            Ts = cellmap(@(F,S) ... % for each filter and seqs list
                cellcat(cellmap(@(O) ... % for each observations
                cellcat(cellmap(@(i) selectTracks(O,i),S)), ... % Collect all tracks from those sequences
                obj.obs(F))), makecell(filter),makecell(seqs));
            
            %             if numel(Ts)==1
            %                 Ts = Ts{1};
            %             end
        end
         function make_normIgECh(obj,filter,frameInfo)
            % function generates normalized IgE images for each observation sequence
            % associated with the files in the LifetimeAnalysis object. 
            % Optional 'filter' input to only generate for subset of
            % observations
            % Defaults are that it uses the #2 sequence first frames to
            % normalize the data
            if nargin > 1
                obsIds = find(filter);
            elseif nargin == 1
                obsIds = 1:numel(obj.obs);
            end
         
            for ii = obsIds
                % -- get starting fluorescence to normalize data from
                nseqIdx = 1; % use the IgE image from the second movie on file to normalize the image
                % load in rpt #2
                if nargin<3
                    ige_start = LifetimeAnalysis.getIgEChannel(obj.obs(ii), nseqIdx);
                elseif nargin == 3
                ige_start = LifetimeAnalysis.getIgEChannel(obj.obs(ii), nseqIdx,frameInfo);
                end
                startval = mean(sum(ige_start(:,:,10:10:30),[],[1 2])); % 
                % -- now use this startval to normalize all the IgE data
                igenames = cell(1,numel(obj.obs(ii).seqTracks)); % cell array of ige file names to add to database
                for seqIdx = 1:numel(obj.obs(ii).seqTracks)
                    if nargin<3
                        ige_dcp = LifetimeAnalysis.getIgEChannel(obj.obs(ii), seqIdx);
                    elseif nargin == 3
                    ige_dcp = LifetimeAnalysis.getIgEChannel(obj.obs(ii), seqIdx,frameInfo);
                    end

                    % only take in frames where laser was on
                    ige=[]; corrfactor_arr = [];
                    for tt=10:10:size(ige_dcp,3) % loop through and make new IgE movie from only DutyCyclePulse Frames
                        corrfactor = startval/double(sum(ige_dcp(:,:,tt-1),[],[1 2]));
                        corrfactor_arr = [corrfactor_arr; repmat(corrfactor,[10 1])];
                        ige=cat(3,ige,double(repmat(ige_dcp(:,:,tt-1).*corrfactor,[1 1 10])));
                        %                     ige_nc=cat(3,ige_nc,repmat(ige_dcp(:,:,tt-1),[1 1 10]));
                    end
                    if ~isempty(ige)
                    % now properly save image with proper image name
                    igenames{seqIdx} = LifetimeAnalysis.getIgEfileName(obj.obs(ii), seqIdx);
                    
                    save(igenames{seqIdx},'ige','corrfactor_arr');
                    end
                end
                f = LifetimeAnalysis.plotIgE_correctionFactorArray(obj.obs(ii));
                close(f);
              obj.obs(ii).seqIgEImage = igenames;
                % then add property to obj observation that puts in the filenames
            end
        end
        
        function missingIgE = linkIgEfilenames(obj,filter)
            % pass in the object and the filter of observation to link to
            % the appropriate IgE files (should be in IgE subfolder)
            % also returns the observation index of any files that should
            % be there but are missing
            if nargin == 2
                obsIds = find(filter);
            elseif nargin == 1
                obsIds = 1:numel(obj.obs);
            end
            missingIgE = [];
            for ii = obsIds
                IgEfilenames = {};
                for seqIdx = 1:numel(obj.obs(ii).seqTracks)
                    IgEfilenames{seqIdx} = LifetimeAnalysis.getIgEfileName(obj.obs(ii), seqIdx);
                end
                if any(cell2mat(cellfun(@(x) ~exist([x '.mat'],'file'), IgEfilenames,'UniformOutput',false)))
                  missingIgE = [missingIgE; ii];
                  fprintf('No IgE Channel Image Information for observation # %g\n',ii)
                end
                obj.obs(ii).seqIgEImage = IgEfilenames;
            end
        end
        
        function IntensityChangePoint(obj,logBayesThreshold,filter)
            if nargin == 3
                obsIds = find(filter);
            elseif nargin == 2
                obsIds = 1:numel(obj.obs);
            end
            objFlag = 0; %flag for if IntensityCPA object has been created or not.
            for ii = obsIds
                %--- create change point booleans for seqTracks
                chngptIDseq = cell(1,numel(obj.obs(ii).seqTracks));%first initialize a cell array for ID booleans
                for seqIdx = 1:numel(obj.obs(ii).seqTracks)
                    if isempty(obj.obs(ii).seqTracks{seqIdx})
                        chngptIDseq{seqIdx} = [];
                        continue;
                    end
                    ds = cellfun(@(x) size(x,1), obj.obs(ii).seqTracks{seqIdx}); %only care about number of localizations - not actual time
                    ts = cellmap(@(x) ceil(x(:,4)), obj.obs(ii).seqTracks{seqIdx}); %intensity trace
                    chngptIDseq{seqIdx} = true(size(ts)); %first initialize all IDs to be one.
                    for tt = 1:size(ts,2)
                        if ds(tt) > 3
                            if ~objFlag
                                iCPA=IntensityCPA(ts{tt}, logBayesThreshold);
                                objFlag = 1;
                            else
                                iCPA.setData(ts{tt},logBayesThreshold);
                            end
                            if iCPA.changePoints>0
                                chngptIDseq{seqIdx}(tt) = false;
                            end
                        else
                            continue
                        end
                    end
                end
                %--- create change point booleans for preTracks
                % assumes only 1 pre file for each observation 
                if ~isempty(obj.obs(ii).preTracks)
                    ds = cellfun(@(x) size(x,1), obj.obs(ii).preTracks); %only care about number of localizations - not actual time
                    ts = cellmap(@(x) ceil(x(:,4)), obj.obs(ii).preTracks); %intensity trace
                    chngptIDpre = true(size(ts)); %first initialize all IDs to be one.
                    for tt = 1:size(ts,2)
                        if ds(tt) > 3
                            if ~objFlag
                                iCPA=IntensityCPA(ts{tt}, logBayesThreshold);
                                objFlag = 1;
                            else
                                iCPA.setData(ts{tt},logBayesThreshold);
                            end
                            if iCPA.changePoints>0
                                chngptIDpre(tt) = false;
                            end
                        else
                            continue
                        end
                    end   
                else
                    chngptIDpre = [];
                end
                %--- now set these values in the database
                obj.obs(ii).changePointIDseqTracks = chngptIDseq;
                obj.obs(ii).changePointIDpreTracks = chngptIDpre;
                display(['Observation ' num2str(ii) '/' num2str(numel(obj.obs))]);
            end
        end
        
        function make_Masksfrom_normIgE(obj,filter,input_seqIdx)
            obsIds = find(filter);
            for ii = obsIds
                if nargin < 3
                    input_seqIdx=1:numel(obj.obs(ii).seqTracks);
                end
                for seqIdx = input_seqIdx %= 1:numel(obj.obs(ii).seqTracks)
                    if isempty(obj.obs(ii).seqIgEImage{seqIdx})
                        continue
                    end
                    s = load(obj.obs(ii).seqIgEImage{seqIdx});
                    display(['Masking file: ' obj.obs(ii).seqIgEImage{seqIdx}]);
                    [IgEMask, CellMask] = LifetimeAnalysis.masking_makeMasks(s.ige);
                    filename = LifetimeAnalysis.getIgEfileName(obj.obs(ii), seqIdx);
                    save(filename,'CellMask', 'IgEMask', '-append');
                end
                f = LifetimeAnalysis.plotIgE_seriesMeanImIntensity(obj.obs(ii),input_seqIdx);
                close(f);
            end
        end    
        
        function make_LabelandMeasurefrom_normIgE(obj,filter,input_seqIdx)
            obsIds = find(filter);
            for ii = obsIds
                if nargin<3
                    input_seqIdx = 2:numel(obj.obs(ii).seqTracks);
                end
                for seqIdx = input_seqIdx
                    if isempty(obj.obs(ii).seqIgEImage{seqIdx})
                        display(['Skipping file: ' obj.obs(ii).seqIgEImage{seqIdx} ' because it does not exist']);
                        continue;
                    end
                     display(['Labeling and Measuring file: ' obj.obs(ii).seqIgEImage{seqIdx}]);
                    [msr,lM] = LifetimeAnalysis.measureIgEChannel(obj.obs(ii),seqIdx,{'size','mass'},[],1);
                    filename = LifetimeAnalysis.getIgEfileName(obj.obs(ii), seqIdx);
                    save(filename,'msr', 'lM', '-append');
                    f = LifetimeAnalysis.plotIgE_measuredSizevsMass(obj.obs(ii),seqIdx);
                    close(f);
                    % call measure, gives back labele and measurement -
                    % then save in file
                end
            end
        end
        
        function make_AllIgEInfo(obj,filter,seqIdx,frameInfo)
            %  first do a check to make sure all observations from filter have
            if nargin<4
                obj.make_normIgECh(filter);
            else
                obj.make_normIgECh(filter,frameInfo);
            end
            if nargin < 3
                obj.make_Masksfrom_normIgE(filter);
                obj.make_LabelandMeasurefrom_normIgE(filter);
            else
                obj.make_Masksfrom_normIgE(filter,seqIdx);
                obj.make_LabelandMeasurefrom_normIgE(filter,seqIdx);
            end
        end
        
        function make_onlyLabelingIgEInfo(obj,filter,seqIdx)
            obj.make_LabelandMeasurefrom_normIgE(filter,seqIdx);
        end
         
         function make_onlyMaskingLabelIgEInfo(obj,filter,seqIdx)
                 obj.make_Masksfrom_normIgE(filter,seqIdx);
                 obj.make_LabelandMeasurefrom_normIgE(filter,seqIdx);
         end
        
       function filter = makeFilter(obj, cellTypes, sykTypes, concs)
            filter = true(1,obj.NTotalObs);
            if nargin>=2 && ~isempty(cellTypes) 
                cellTypes = makecell(cellTypes);
                cellType_filter = false(1,obj.NTotalObs);
                for n=1:numel(cellTypes)
                    cellType_filter = cellType_filter | cellfun(@(T) strcmpi(cellTypes{n},T), {obj.obs(:).cellType});
                end
                filter = filter & cellType_filter;
            end
            if nargin>=3 && ~isempty(sykTypes)
                sykTypes = makecell(sykTypes);
                sykType_filter = false(1,obj.NTotalObs);
                for n=1:numel(sykTypes)
                    sykType_filter = sykType_filter | cellfun(@(T) strcmpi(sykTypes{n},T), {obj.obs(:).sykType});
                end
                filter = filter & sykType_filter;
            end
            if nargin>=4 && ~isempty(concs)
                conc_filter = false(1,obj.NTotalObs);
                cs = [obj.obs(:).conc];
                for n=1:numel(concs)
                    conc_filter = conc_filter | (cs == concs(n));
                end
                filter = filter & conc_filter;
            end
        end
        
        function val = get.NTotalObs(obj)
            val = numel(obj.obs);
        end
    end

    methods (Static=true)
        function tfilt = trackFilterIntensity(Ts, max_intensity)
            tfilt = cellfun(@(T) max(T(:,4))<=max_intensity, Ts); 
        end
        function ds = trackDurationFrames(Ts)
            ds = cellfun(@(T) T(end,end)-T(1,end)+1, Ts); 
        end
        function Is = trackIntensities(Ts)
            Is = cellmap(@(T) T(:,4)', Ts); 
        end
        function Is = trackMeanIntensities(Ts)
            Is = cellfun(@(T) mean(T(:,4)), Ts); 
        end
        function Is = trackMaxIntensities(Ts)
            Is = cellfun(@(T) max(T(:,4)), Ts); 
        end

        function [outmask] = masking_image_threshingDCC(in_im,thresh_id)  
            im = smooth(in_im,[2 2 0]);
            outmask = false(size(logical(im)));
            for ii =  1:size(im,3)
                laplace_smooth_im = dcc(im(:,:,ii-1));
                thr = multithresh(single(-laplace_smooth_im),thresh_id);
            om = -laplace_smooth_im > thr(thresh_id);
            outmask(:,:,ii) = logical(om);   
            end
        end
        
        function out_mask=masking_image_threshingNiblack(imagetothreshold,gaussianfilter,localthresholdradius) %this needs to be improved - not actually a local operations
            % loosely related to M.Meddens et al 2012 Microscopy and
            % Microanalysis for Podosome detection
            
            % gaussian filter
            gaussimage=gaussf(imagetothreshold,[gaussianfilter gaussianfilter 0]);
            out_mask = false(size(im));
            % laplacian of gaussian filter
            for ii = 1:size(imagetothreshold,3)
            laplaceimage=gaussimage(:,:,ii-1)-laplace(gaussimage(:,:,ii-1));
            meanimgval = unif(laplaceimage,[localthresholdradius localthresholdradius 0],'elliptic');
            varofimg  = varif(laplaceimage,[localthresholdradius localthresholdradius 0],'elliptic');
            LocalThreshold = meanimgval + (0.5 * (varofimg^0.5));
            out_mask(:,:,ii) = round(laplaceimage > LocalThreshold);
            end
%             laplaceimage=gaussimage-laplace(gaussimage,[1 1 0]);
%             % calculate mean and variance of image
%             meanimgval = unif(1*laplaceimage,[localthresholdradius localthresholdradius 0],'elliptic');
%             varofimg  = varif(1*laplaceimage,[localthresholdradius localthresholdradius 0],'elliptic');
%             % use mean and variance to calculate local threshold: eq from Meddens
%             % Microscopy and Microanalysis Paper (podosome identification)
%             LocalThreshold = meanimgval + (0.5 * (varofimg^0.5));
%             % Apply local threshold
%             out_mask = round(laplaceimage > LocalThreshold);
        end
        
        function L = masking_frame_watershed(in_im,conn)
            % calls the watershed function on the negative of an input matrix (image single frame). 
            % First smooths the image then call
            % this function is to find regions in between clusters so can
            % break apart linked clustes in masking
            im = smooth(in_im);
            if nargin>1
                L = watershed(-im,conn);
            else
                L = watershed(-im);
            end
        end
        
        function L = masking_movie_watershed(im,conn)
            % calls the frame_watershed function using slice_op to calculate watershed of a multiframe movie 3D matrix.
            % frame_watershed operates on the negative of an input matrix (image single frame). 
            % First smooths the image then call
            % this function is to find regions in between clusters so can
            % break apart linked clustes in masking
            if nargin>1
                L = slice_op('LifetimeAnalysis.masking_frame_watershed',im,conn);
            else
                L = slice_op('LifetimeAnalysis.masking_frame_watershed',im);
            end
        end
        
       function [IgEMask, CellMask] = masking_makeMasks(in_im, option, val,watershedconn,localthresholdradius)
           % function to threshold an image and get back the mask and a mask of the cell area
           % outputs: 
           %    IgEMask - binary mask of intensity threshold (logical not dipimage) same size as in_im
           %    CellMask - binary mask of cell output (logical not dipimage) same size as in_im
           % inputs:
           %    in_im: input image to threshold - doesn't need to be dipimage
           %    option: 'LG' for laplacian of the gaussian thresholding
           %        'NiB' for Niblack filtering. Default is LG.
           %    val: depends on option choise: LG - number of thresholds in multithresh call. NiB- gaussianfilter
           %        size. Default is 1 for both.
           %    localthresholdradisu: input for NiB option. Default is 10
           %    watershedconn: input for NiB option because need to break up so many connections (could also
           %        incorporate for LG option? Default is 1.
            if nargin<2
               option = 'LG';
           end
           if nargin<4
               watershedconn = 1;
           end
           thr = multithresh(single(gaussf(in_im,[2 2 0])),5);
           aa = gaussf(in_im,[2 2 0])>thr(1);
           CellMask = logical(closing(aa,7));
           L = LifetimeAnalysis.masking_movie_watershed(in_im,watershedconn); 
           switch option
               case 'LG'
                   if nargin<3
                       numthreshs =1;
                   else
                       numthreshs = val;
                   end
                   out_im = LifetimeAnalysis.masking_image_threshingDCC(in_im,numthreshs);
                   out_im(L) = 0;
                   IgEMask = logical(out_im.*CellMask);
               case 'NiB'    
                   if nargin<5
                       localthresholdradius = 10;
                   end
                   if nargin<3
                       gaussianfilter = 1;
                   end
                   if nargin>=3
                       gaussianfilter = val;
                   end
                   %--- second threshold out IgE areas
                   %
                   out_im=LifetimeAnalysis.masking_image_threshingNiblack(in_im,gaussianfilter,localthresholdradius);
                   out_im(L) = 0;
                   IgEMask = logical(out_im.*CellMask);
           end
       end
       
       function [lM] = labelIgEChannel(obs, seqIdx)
           % calls 'getIgEChannelMask' and then labels mask using
           % 'label_3DImage' to ensure that labels in between frames aren't
           % connected.
           [im] = LifetimeAnalysis.getIgEMask(obs, seqIdx);
           lM = LifetimeAnalysis.label_3DImage(dip_image(im));
       end
       
       function lM = label_3DImage(object_in)
           % This function calls the label function on a mask time series. It insures that masks between
           % frames are not connected (in time). Slice-op call of the measure function
           % doesn't work because mask numbering/naming is repeated per frame. Should
           % just be able to call normal measure function on this labelled lM image.
           % this function creates a measurement
           lM= slice_op('label',object_in,1);
           lbs =single(unique(lM(:,:,0)));
           addon = max(lbs);
           for fr = 1:size(lM,3)-1
               tempframe = lM(:,:,fr);
               tempframe(find(tempframe)) = tempframe(find(tempframe)) + addon;
               lM(:,:,fr) = tempframe;
               lbs =single(unique(lM(:,:,fr)));
               addon = max(lbs);
           end
       end
       
       function [msr,lM] = measureIgEChannel(obs,seqIdx,measurementIDs,objectIDs,...
               connectivity,minSize,maxSize)
           % this function calls 'labelIgEChannel' and also gets the
           % normalized IgE channel using 'getNormedIgEChannel' to measure
           % the time series mask labels. Other inputs are normal inputs
           % for measure function.
           % set defaults:
           if nargin < 7
               maxSize = 0;
           end
           if nargin < 6
               minSize = 0;
           end
           if nargin < 5
               connectivity = 1;
           end
           if nargin < 4
               objectIDs = [];
           end
           if nargin < 3
               measurementIDs = 'size';
           end 
           fprintf('Labeling IgE Channel\n');
           [lM] = LifetimeAnalysis.labelIgEChannel(obs, seqIdx);
           fprintf('Extracting Normed IgE Channel\n');
           normIgEim = dip_image(LifetimeAnalysis.getNormedIgEim(obs, seqIdx));
           fprintf('Measuring IgE Channel\n');
           msr = measure(lM,normIgEim,measurementIDs,objectIDs,...
               connectivity,minSize,maxSize);
       end
       
       function f = plotIgE_measuredSizevsMass(obs,seqIdx)
           filename = LifetimeAnalysis.getIgEfileName(obs, seqIdx);
           s = load(filename);
           f = figure; 
           plot(s.msr.Size,s.msr.Mass,'.');
           xlabel('Aggregate Size (pixels)'); ylabel('Aggregate Mass (AU)'); title('Masking');
           [PATHSTR,NAME,EXT] = fileparts(LifetimeAnalysis.getIgEfileName(obs, seqIdx));
            figdir = fullfile(PATHSTR,'MeasuredSizevsMass');
            if ~exist(figdir)
                mkdir(figdir);
            end
            saveas(f,fullfile(figdir,[NAME '_MeasuredSizevsMass']),'png');
            saveas(f,fullfile(figdir,[NAME '_MeasuredSizevsMass']),'fig');
       end
        function f = plotIgE_seriesMeanImIntensity(obs,seqIdx)
            meanIntensity = [];
            for sId = seqIdx
                if isempty(LifetimeAnalysis.getIgEfileName(obs, sId))
                    continue;
                end
             s = load(LifetimeAnalysis.getIgEfileName(obs, sId));
             mI = squeeze(double(mean(dip_image(s.ige),dip_image(s.CellMask),[1 2])));
%              mI = squeeze(mean(mean(s.ige.*IgEMask,1),2));
             meanIntensity = cat(1,meanIntensity,mI);
            end
            f = figure; plot(meanIntensity);
            xlabel('Frames'); ylabel('Mean Intensity Inside Cell Mask');  
            [PATHSTR,NAME,EXT] = fileparts(LifetimeAnalysis.getIgEfileName(obs, 1));
            figdir = fullfile(PATHSTR,'MeanIntensities');
            if ~exist(figdir)
                mkdir(figdir);
            end
            saveas(f,fullfile(figdir,[NAME '_MeanInCellIntensity']),'png');
            saveas(f,fullfile(figdir,[NAME '_MeanInCellIntensity']),'fig');
        end
        
        function [f, N, bns] = plotinputHistogram(Ds,names,frameTime,clmap)
            if nargin < 3
                frameTime = 1;
                clmap = lines;
            end
             if nargin < 4
                clmap = lines;
             end
            f = figure; hold('on'); fa = gca;
            N = cell(numel(Ds),1);
            for rr = 1:numel(Ds)
                if rr == 1
                [N{rr},edges] = histcounts(Ds{rr}.*frameTime,'Normalization','probability');
                if numel(edges)<80
                   [N{rr},edges] = histcounts(Ds{rr}.*frameTime,80,'Normalization','probability'); 
                end
                binwidth = diff(edges(1:2));
                bns = edges(2:end)-binwidth;
                else
                    [N{rr},edges] = histcounts(Ds{rr}.*frameTime,edges,'Normalization','probability');
                end
                hst = plot(bns,N{rr});
%                 HH.DisplayName = 'Data';
                hst.LineWidth = 1.5;
                hst.LineStyle = '-';
                hst.Color = clmap(rr,:);
                hst.DisplayName = names{rr};
            end
            ylabel('Probability Density');
            if nargin>=3
                xlabel('Duration (Secs)');
            elseif nargin<3
                xlabel('Duration (frames)');
            end
            legend('location','best');
            fa.XScale='log';
            fa.YScale='log';
            fa.Title.String = sprintf('PDF of Track Lengths'); 
        end
        
        function [f,N,x] = plotinputCDF(Ds,names,frameTime,clmap)
            if nargin < 3
                frameTime = 1;
                clmap = lines;
            end
            if nargin < 4
                clmap = lines;
            end
            f = figure; hold('on'); fa = gca;
            x = cell(numel(Ds),1);
            N = cell(numel(Ds),1);
            for rr = 1:numel(Ds)
                [x{rr},N{rr}] = ecdf(Ds{rr}.*frameTime);
                cline = plot(N{rr}(2:end),x{rr}(2:end));
                cline.LineWidth = 1.5;
                cline.LineStyle = '-';
                cline.Color = clmap(rr,:);
                cline.DisplayName = names{rr};
            end
            ylabel('Probability Density');
            if nargin>=3
                xlabel('Duration (Secs)');
            elseif nargin<3
                xlabel('Duration (frames)');
            end
            legend('location','best');
            fa.XScale='log';
            fa.Title.String = sprintf('CDF of Track Lengths'); 
        end
        
        function f = plotIgE_correctionFactorArray(obs)
            totcorarr = [];
            for seqIdx = 1:numel(obs.seqTracks)   
                igename = LifetimeAnalysis.getIgEfileName(obs, seqIdx);
                if ~isempty(igename)
                    s = load(LifetimeAnalysis.getIgEfileName(obs, seqIdx));
                else
                    s.corrfactor_arr = zeros(40,1); % file doesn't actually exist so just put some zeros for plotting
                end
            totcorarr = [totcorarr;s.corrfactor_arr];
            end
            f = figure; plot(totcorarr);
            xlabel('Frame'); ylabel('Multiplicative Correction Factor');
            [PATHSTR,NAME,EXT] = fileparts(LifetimeAnalysis.getIgEfileName(obs, 1));
            figdir = fullfile(PATHSTR,'CorrectionFactor');
            if ~exist(figdir)
                mkdir(figdir);
            end
            saveas(f,fullfile(figdir,[NAME '_corrfactorArray']),'png');
            saveas(f,fullfile(figdir,[NAME '_corrfactorArray']),'fig');
        end
              
        function im = getNormedIgEim(obs, seqIdx) % this maybe obsolete
            % important: returns image in standard image format
            % returns this as a single array. call dip_image(im) to
            % get dipimage. Also need to index in the other way
            IgEfilename = LifetimeAnalysis.getIgEfileName(obs, seqIdx);
            s = load(IgEfilename);
            im = s.ige;
        end
        
        function msk = getIgECellMask(obs,seqIdx)
            IgEfilename = LifetimeAnalysis.getIgEfileName(obs, seqIdx);
            s = load(IgEfilename);
            msk = s.CellMask;
        end
        
        function [msk] = getIgEMask(obs,seqIdx)
           IgEfilename = LifetimeAnalysis.getIgEfileName(obs, seqIdx);
           msk = [];
           if ~isempty(IgEfilename)
            s = load(IgEfilename);
            msk = s.IgEMask; 
           end
        end
        
        function [lM] = getLabeledIgEim(obs,seqIdx)
            IgEfilename = LifetimeAnalysis.getIgEfileName(obs, seqIdx);
            s = load(IgEfilename);
            lM = s.lM;
        end
        
        function [msr] = getMeasureIgEim(obs,seqIdx)
            IgEfilename = LifetimeAnalysis.getIgEfileName(obs, seqIdx);
            s = load(IgEfilename);
            msr = s.msr; 
        end
        
        function ige = getIgEChannel(obs, seqIdx, frameInfo)
            % important: returns image in standard image format
            % returns this as a single array. call dip_image(im) to
            % get dipimage. Also need to index in the other way
            % frameInfo: optional variable that signifies the x, y range to use for the IgE Channel. Must match
            % what was used for the ChannelRegistration object!! frameInfo(1,:) is x range. frameInfo(2,:) is y
            % range. if not included then just assume that you are using 1/2 the ccd.
            rptFileName = LifetimeAnalysis.getRPTSequenceFileName(obs,seqIdx);
            if isempty(rptFileName)
                display(['No file for ' obs.v ' file at series ' num2str(seqIdx)]);
                ige = [];
                return;
            end
            rr = RPT(rptFileName);
            % now get SPData from this
            spd = SPData(rr.data.saveFilePath);
            filename= spd.saveFileBaseName;
            datf = load(fullfile(spd.rawDataPath));
            if nargin>2 && size(datf.sequence,1) == 512
                ige = (datf.sequence(frameInfo(1,:),frameInfo(2,:),:)-spd.CCDBackground)/spd.CCDGain;
            else
            ige = (datf.sequence(0:floor(spd.frameSize(1)/2)-1,:,:)-spd.CCDBackground)/spd.CCDGain;
            end
        end
        
         function IgEfilename = getIgEfileName(obs, seqIdx)
            rptFileName = LifetimeAnalysis.getRPTSequenceFileName(obs,seqIdx);
            if isempty(rptFileName)
                IgEfilename = [];
                return;
            end
            rr = RPT(rptFileName);
            % now get SPData from this
            spd = SPData(rr.data.saveFilePath);
            filename= [spd.saveFileBaseName '_IgECh'];
            igedir = fullfile(spd.workingDir,'IgE');
            if ~exist(igedir)
                mkdir(igedir)
            end
            IgEfilename = fullfile(igedir,filename); 
        end
        
        
         function [inMask_tracks,outMask_tracks] = compareTrackswithMask(obs,seqIdx,Ithresh,inmaskThreshFrames)
            % [inMask_tracks,outMask_tracks] = compareTracksMask(obs,seqIdx)
            % works for a single observation of database (obs). pass in
            % sequence id: example [2 3]
            inMask_tracks = [];
            outMask_tracks = [];
            for ii =1:numel(obs)
                for ss = 1:numel(seqIdx)
                    if isempty(obs(ii).seqTracks{seqIdx(ss)})
                        continue;
                    end
                ts = LifetimeAnalysis.shiftTracksPixels(obs(ii),seqIdx(ss));
%                 if isempty(ts)
%                      continue;
%                 end
                igemask = LifetimeAnalysis.getIgEMask(obs(ii),seqIdx(ss));
                inframe = cellmap(@(T) T(:,2)>=0 & T(:,2)<= size(igemask,2) & T(:,3)>=0 & T(:,3)<= size(igemask,1), ts);
                outofframe = cellfun(@(f) any(f==0), inframe); % boolean. true = entire track is out of frame.
                inframe = inframe(~outofframe); % remove tracks that are entirely out of frame
                ts = ts(~outofframe); % remove tracks that are entirely out of frame
                if nargin==3
                    tfilt = LifetimeAnalysis.trackFilterIntensity(ts, Ithresh);
                    ts = ts(tfilt);
                    inframe=inframe(tfilt);
                end
                inmask = cellmap(@(T,f) igemask(sub2ind(size(igemask),ceil(T(:,3)),ceil(T(:,2)),T(:,end))),ts,inframe);
                if nargin>3
                    inMask_tracksId = cellfun(@(x) sum(x)/size(x,1)>inmaskThreshFrames,inmask);
                else
                    inMask_tracksId = cellfun(@any, inmask);
                end
                inMask_tracks = [inMask_tracks, ts(inMask_tracksId)];
                outMask_tracks = [outMask_tracks, ts(~inMask_tracksId)];
                end
            end
         end
         
        function [f,fa] = plotCDF_compareTrackswithMask(obs,seqIdx)
            [inMask_tracks,outMask_tracks] = LifetimeAnalysis.compareTrackswithMask(obs,seqIdx);
            inmask_tracklength = cellfun(@ (T) T(end,end)-T(1,end) , inMask_tracks);
            outmask_tracklength = cellfun(@ (T) T(end,end)-T(1,end) , outMask_tracks);
            f = figure; fa = axes('Parent',f); hold on;
            cdfplot(inmask_tracklength);
            cdfplot(outmask_tracklength);
            legend({'Tracks in Mask','Tracks out of Mask'});
            xlabel('Track Length (Frames)'); ylabel('CDF');
            
        end
        
        function [tracklengths, pixelvalues] = compareTrackwithMaxPixelIntensity(obs,seqIdx,Ithresh)
            %[tracklengths, pixelvalues] = compareTrackwithMaxPixelIntensity(obs,seqIdx)
            % only include tracks that are within the IgE image after
            % shifting. 
            ts = LifetimeAnalysis.shiftTracksPixels(obs,seqIdx);
            igemask = LifetimeAnalysis.getNormedIgEChannel(obs,seqIdx);
            inframe = cellmap(@(T) T(:,2)>=0 & T(:,2)<= size(igemask,2) & T(:,3)>=0 & T(:,3)<= size(igemask,1), ts);
            outofframe = cellfun(@(f) ~any(f), inframe); % boolean. true = entire track is out of frame. 
            inframe = inframe(~outofframe); % remove tracks that are entirely out of frame
            ts = ts(~outofframe); % remove tracks that are entirely out of frame
            if nargin==3
                tfilt = LifetimeAnalysis.trackFilterIntensity(ts, Ithresh);
                ts = ts(tfilt);
                 inframe=inframe(tfilt);
            end
            igemaskpixel = cellmap(@(T,f) igemask(sub2ind(size(igemask),ceil(T(f,3)),ceil(T(f,2)),T(f,end))),ts,inframe);
            pixelvalues = cellfun(@max, igemaskpixel);
            tracklengths = cellfun(@ (T) T(end,end)-T(1,end) , ts);
        end
        
        function [trackduration_frames,maskAreas,maskMass] = compareTrackwithMeasuredIgECh(obs,seqIdx,Ithresh)
            % need to make sure all the pre thresholding of tracks works
            % here
            % calculates the associated mask size and mask mass per frame
            % for each track. if the track is not in a mask at that frame
            % then the value defaults to 0
%             Ts = shiftTracksPixels(obj,Ts)%--------------------------------------------!!!!!!!!!!!!!!!!!!!!!
            ts = LifetimeAnalysis.shiftTracksPixels(obs,seqIdx);
            [lM] = LifetimeAnalysis.getLabeledIgEim(obs,seqIdx);
            [msr] = LifetimeAnalysis.getMeasureIgEim(obs,seqIdx);
            igemask = double(lM);
            inframe = cellmap(@(T) T(:,3)>=0 & T(:,3)<= size(igemask,1) & T(:,2)>=0 & T(:,2)<= size(igemask,2), ts);
            outofframe = cellfun(@(f) ~any(f), inframe); % boolean. true = entire track is out of frame.
            inframe = inframe(~outofframe); % remove tracks that are entirely out of frame
            ts = ts(~outofframe); % remove tracks that are entirely out of frame
            if nargin==3
                tfilt = LifetimeAnalysis.trackFilterIntensity(ts, Ithresh);
                ts = ts(tfilt);
                inframe=inframe(tfilt);
            end
            
            trackduration_frames = cellfun(@ (T) T(end,end)-T(1,end)+1 , ts);
            maskMass = cellmap(@(x) zeros(size(x,1),1),ts);
            maskAreas = cellmap(@(x) zeros(size(x,1),1),ts);
            
            for tt = 1:numel(ts)
                track = ts{tt};
                tframes = size(track, 1);
                for locid = 1:tframes
                    if inframe{tt}(locid)
                        maskclusterId = igemask(ceil(track(locid,3)),ceil(track(locid,2)),track(locid,end));
                        if maskclusterId
                            maskMass{tt}(locid) = msr(maskclusterId).Mass;
                            maskAreas{tt}(locid) = msr(maskclusterId).Size;
                        end
                    end
                end
                
            end
        end
        
        function [f,all_tracklengths,all_maskmasses,tot_Tracks]  = plotTrackLengthvsMaskMass_averageOverFrames(obs,seqIdx,Ithresh,frac_inmaskframes)
            all_meanmasses = [];
            all_tracklengths = [];
            tot_Tracks = 0;
            all_maskmasses = [];
            for oo = 1:numel(obs)
                for ss = seqIdx
                    if ~isempty(obs(oo).seqTracks{ss})
                    [tlenarray,~,maskSums] = LifetimeAnalysis.compareTrackwithMeasuredIgECh(obs(oo),ss,Ithresh);
                    okid = cell2mat(cellmap(@(x) nnz(x)./numel(x)>=frac_inmaskframes, maskSums));
%                     okid = cell2mat(cellmap(@(x) sum(x>0)>=min_inmaskframes, maskSums));
                    maskSums_ok = maskSums(okid);
                    tracklengths = tlenarray(okid);
                    maskmass = cell2mat(cellmap(@(x) mean(x(x>0)),maskSums_ok)); %-- does this work?
                    all_meanmasses = [all_meanmasses, maskmass];
                    all_tracklengths = [all_tracklengths, tracklengths];
                    for tt = 1:numel(maskSums_ok) % add to array all every mask
                        currMS = maskSums_ok{tt};
                        all_maskmasses = [all_maskmasses; currMS(currMS>0)];
                    end
                    tot_Tracks = tot_Tracks + numel(tlenarray);
                    end
                end
            end
            f = figure; plot(all_tracklengths,all_meanmasses,'.');   
            f.UserData = all_maskmasses;
            xlabel('Tracklengths (frames)');ylabel('Mean of Associated Mask Mass over time');
            title([mat2str(numel(all_tracklengths)) ' tracks Inside Mask out of ' mat2str(tot_Tracks) ' total tracks']);
        end
        
        function [f,all_tracklengths,all_masses,tot_Tracks]  = plotTrackLengthvsMaskMass_MaxOverFrames(obs,seqIdx,Ithresh,frac_inmaskframes)
            all_masses = [];
            all_tracklengths = [];
            tot_Tracks = 0;
            for oo = 1:numel(obs)
                for ss = seqIdx
                    [tlenarray,~,maskSums] = LifetimeAnalysis.compareTrackwithMeasuredIgECh(obs(oo),ss,Ithresh);
                    okid = cell2mat(cellmap(@(x) nnz(x)./numel(x)>=frac_inmaskframes, maskSums));
                    maskmass = cell2mat(cellmap(@(x) max(x),maskSums(okid)));
                    tracklengths = tlenarray(okid);
                    all_masses = [all_masses, maskmass];
                    all_tracklengths = [all_tracklengths, tracklengths];
                    tot_Tracks = tot_Tracks + numel(tlenarray);
                end
            end
            f = figure; plot(all_tracklengths,all_masses,'.');   
            xlabel('Tracklengths (frames)');ylabel('Max of Associated Mask Mass over time');
            title([mat2str(numel(all_tracklengths)) ' tracks Inside Mask out of ' mat2str(tot_Tracks) ' total tracks']);
        end
        
        function h = sim_Randomlypair_MaskMasswithtrackLengths(figfilename)
            % loads in the figfilename. extracts the XData and YData.
            % YData is the means of the mask masses.
            % for each XData - choose a random number to index into YData
            % and plot
            uiopen(figfilename,1);
            h = gcf;
            a1 = get(gcf, 'Children');
            a2 = get(a1, 'Children');
            TrackLengths = a2.XData;
            MaskMasses = sort(a2.YData);
            rmin = 1;
            rmax = numel(MaskMasses);
            randMaskMasses = zeros(1,numel(TrackLengths));
            for tt = 1:numel(TrackLengths) % loop through each track
                MMid = randi([rmin rmax]);
                randMaskMasses(tt) = MaskMasses(MMid);
            end
            hold on; plot(TrackLengths,randMaskMasses,'.r');
            xlabel('Tracklengths (frames)');ylabel('Random Mean of Associated Mask Mass over time');
            title('Randomly Selected (from original data) Mask Means');
            set(gca,'YScale','log');
            set(gca,'Ylim',[10e3 10e7])
            legend({'Real Data','Randomized Data'});
        end
        
        function newsim_obs = sim_TracklengthsFromMaskMass(obs,seqIdx,Ithresh)
            newsim_obs = [];
            for oo = 1:numel(obs)
                for ss = 1:numel(seqIdx)
                    Ts = LifetimeAnalysis.shiftTracksPixels(obs(oo),seqIdx(ss));
                    D = diffusionConstMLE_ensemble(Ts);
                    % get starting 
                    
                    
                    % for this movie - get average Syk diffusion
                    % get IgE Mask
                    % get starting Syk value
                end
            end
        end
        
        function newsim_obs = sim_oldTracklengthsFromMaskMass(obs,seqIdx,Ithresh)
           % obs is structure array of observations
           % Take trajectories and have their tracklengths be a function of the mask intensity they are associated with. 
           % No mobility so x,y positions and all other values in tracks matrix stay the same. 
           % Returns a modified version of the obs variable with the new tracks - so that this can be directly run
           % by the function plotTrackLengthvsMaskMass_averageOverFrames
           for oo = 1:numel(obs)
               for ss = 1:numel(seqIdx)
                   % first only take the tracks that fit the input thresholding params
                   ts = LifetimeAnalysis.shiftTracksPixels(obs(oo),seqIdx(ss));
                   [lM] = LifetimeAnalysis.getLabeledIgEim(obs(oo),seqIdx(ss));
                   [msr] = LifetimeAnalysis.getMeasureIgEim(obs(oo),seqIdx(ss));
                   igemask = double(lM);
                   inframe = cellmap(@(T) T(:,3)>=0 & T(:,3)<= size(igemask,1) & T(:,2)>=0 & T(:,2)<= size(igemask,2), ts);
                   outofframe = cellfun(@(f) ~any(f), inframe); % boolean. true = entire track is out of frame.
                   inframe = inframe(~outofframe); % remove tracks that are entirely out of frame
                   ts = ts(~outofframe); % remove tracks that are entirely out of frame
                   if nargin==3
                       tfilt = LifetimeAnalysis.trackFilterIntensity(ts, Ithresh);
                       ts = ts(tfilt);
                       inframe=inframe(tfilt);
                   end
                   % loop through and get mask masses
                   maskMass = cellmap(@(x) zeros(size(x,1),1),ts);
                   maskAreas = cellmap(@(x) zeros(size(x,1),1),ts);
                   for tt = 1:numel(ts)
                       track = ts{tt};
                       for locid = 1:tframes
                           if inframe{tt}(locid)
                               maskclusterId = igemask(ceil(track(locid,3)),ceil(track(locid,2)),track(locid,end));
                               if maskclusterId
                                   maskMass{tt}(locid) = msr(maskclusterId).Mass;
                                   maskAreas{tt}(locid) = msr(maskclusterId).Size;
                               end
                           end
                       end
                       
                   end
                   % now create trajectoy from the mean of the mask mass (for each track) and save into new obs
                   meanMass = cellmap(@(M) mean(M(M>0)), maskMass);
                   maxframe = max(cell2mat(cellmap(@(T) max(T(:,end)), ts)));
                   for tt = 1:numel(ts)
                       old_track = ts{tt};
                       tlen = geornd(meanMass);
                       tinit = old_track(1,:);
                       new_track = repmat(tinit,tlen,1);
                       frames = tinit(1,end): tinit(1,end) + tlen;
                       new_track(:,end) = frames';
                       obs(oo).seqTracks{ss}{tt} = new_track; 
                   end
               end
           end
           newsim_obs = obs;
        end
        
        function h = sim_RandomlyGen_maskMassfromEmpCDFwithtrackLengths(figfilename)
            uiopen(figfilename,1);
            h = gcf;
            a1 = get(gcf, 'Children');
            a2 = get(a1, 'Children');
            TrackLengths = a2.XData;
            MaskMasses = sort(h.UserData);
            % gen cdf and then draw from it
            % code adapted from emprand function on Matlab FileExchange
            num_mm = numel(MaskMasses);
            x = sort(MaskMasses);
            p = 1:num_mm;
            mm_cdf = p./num_mm;
            rn =  rand(numel(TrackLengths),1);
            % Interpolate ur from empirical cdf and extraplolate for out of range
            % values. 
            MMdrawnfromEmpCDF = interp1(mm_cdf,x,rn,[],'extrap');
            hold on; plot(TrackLengths,MMdrawnfromEmpCDF,'.r');
            xlabel('Tracklengths (frames)');ylabel('Random drawn Mask Mass from Empirical CDF');
            title('Randomly Drawn from CDF (from original data)');
            set(gca,'YScale','log');
            set(gca,'Ylim',[10e3 10e7])
            legend({'Real Data','Randomized Data'});
        end
        
        function [f,fa] = plotCDF_compareTrackwithMaxPixelIntensity(obs,seqIdx,Ithresh)
            [tracklengths, pixelvalues] = LifetimeAnalysis.compareTrackwithMaxPixelIntensity(obs,seqIdx,Ithresh);
            f = figure; fa = axes('Parent',f); hold on;
            plot(tracklengths,pixelvalues,'.b');
            xlabel('Track Length (Frames)'); ylabel('Max Normalized Pixel Value Traveled Through');
        end
        
        function [tracklengths, pixelvalues] = compareTrackwithMeanPixelIntensity(obs,seqIdx,Ithresh)
            ts = LifetimeAnalysis.shiftTracksPixels(obs,seqIdx);
            igemask = LifetimeAnalysis.getNormedIgEChannel(obs,seqIdx);
            inframe = cellmap(@(T) T(:,2)>=0 & T(:,2)<= size(igemask,2) & T(:,3)>=0 & T(:,3)<= size(igemask,1), ts);
            outofframe = cellfun(@(f) ~any(f), inframe); % boolean. true = entire track is out of frame. 
            inframe = inframe(~outofframe); % remove tracks that are entirely out of frame
            ts = ts(~outofframe); % remove tracks that are entirely out of frame
            if nargin==3
                tfilt = LifetimeAnalysis.trackFilterIntensity(ts, Ithresh);
                ts = ts(tfilt);
                inframe=inframe(tfilt);
            end
            igemaskpixel = cellmap(@(T,f) igemask(sub2ind(size(igemask),ceil(T(f,3)),ceil(T(f,2)),T(f,end))),ts,inframe);
            pixelvalues = cellfun(@mean, igemaskpixel);
            tracklengths = cellfun(@ (T) T(end,end)-T(1,end) , ts);
        end
        
        function [f,fa] = plotCDF_compareTrackwithMeanPixelIntensity(obs,seqIdx,Ithresh)
            [tracklengths, pixelvalues] = LifetimeAnalysis.compareTrackwithMaxPixelIntensity(obs,seqIdx,Ithresh);
            f = figure; fa = axes('Parent',f); hold on;
            plot(tracklengths,pixelvalues,'.b');
            xlabel('Track Length (Frames)'); ylabel('Mean of Normalized Pixel Values Traveled Through');
        end
        
        function ra = getRegAnalysisObj(obs)
            % returns registration analysis obj that matches most closely in time with the file. must be located in
            % a folder called 'Registrations' two directories above the original data. If the difference between
            % the time the data was taken and when the registration was acquired, an empty object is returned
            
            % find the first corresponding rpt file
            if ~isempty(obs.preRPT)
            rpt = RPT(obs.preRPT);
            else
                rpt = RPT(obs.seqRPT{1});
                display(['Note: no pre file for observation, using ' obs.seqRPT{1} ' to match ChReg file']);
            end
            % find the date/time stamp for this file
            baseName = rpt.data.saveFileBaseName;
            obs_datnum = datenum(LifetimeAnalysis.fileName2dateTime(baseName));
            if isempty(obs_datnum)
                display(['Could not find time stamp of file ' rpt.saveFileBaseName ' to match up to ChReg file']);
                ra = [];
                return;
            end
            % compare to all channel registrations in the directory - this chreg location is hard coded  
            crfold = fullfile(rpt.data.workingDir,'..','..','Registrations');
            crfiles = dir(fullfile(crfold,'*.reganalysis'));
            chreg_datnums = arrayfun(@(x) datenum(LifetimeAnalysis.fileName2dateTime(x.name(1:end-12))), crfiles,'UniformOutput',false);
            [mindiff,crId] = min(abs(cell2mat(chreg_datnums)-obs_datnum));    
            if mindiff>1
                display(['No good associated Channel Registration File in time for ' rpt.saveFileBaseName ', so skipping this file']);
                ra = [];
                return;
            else
            ra = RegistrationAnalysis(fullfile(crfold,crfiles(crId).name));
            end
        end
        
        function dateTime = fileName2dateTime(baseName)
           parts=strsplit(baseName,'-');
            if numel(parts)==6
                dateTime = datetime([cellfun(@str2double,parts(2:6)), 0]);
            elseif numel(parts)>6
                dateTime = datetime(cellfun(@str2double,parts(end-5:end)));
            else
                dateTime = [];
                return
            end 
        end
        
        
        function tracks = shiftTracksPixels(obs,seqIdx)
            % tracks = shiftTracksPixels(obs,seqIdx) 
             rpt = RPT(LifetimeAnalysis.getRPTSequenceFileName(obs,seqIdx));
             ra = LifetimeAnalysis.getRegAnalysisObj(obs);
             if isempty(ra)
                 tracks = [];
                 return;
             end
             tracks = ra.shiftRPTTracksPixels(rpt);  
        end
        
        function tm = makeTrackMovieOverlay_NormIgECh(obs,seqIdx)
            % tm = makeTrackMovieOverlay(obs,seqIdx)
            % takes in single observation and seqId and returns the
            % TrackMovie object to call tm.viewSequence on
            im_raw = LifetimeAnalysis.getNormedIgEim(obs, seqIdx);
            tracks = LifetimeAnalysis.shiftTracksPixels(obs,seqIdx);   
            frameBounds = [min(cellfun(@(T) min(T(:,end)),tracks)), max(cellfun(@(T) max(T(:,end)),tracks))];
            xBounds = [min(cellfun(@(T) min(T(:,2)),tracks)), max(cellfun(@(T) max(T(:,2)),tracks))];
            im = im_raw(:,:,frameBounds(1):frameBounds(2));
            tm = TrackMovie(tracks,im);
            tm.setTrackColorMethod('TrackLength',logSpaceColorMap(@jet));        
        end

        function tm = makeTrackMovieOverlay_LHS(obs,seqIdx)
            % tm = makeTrackMovieOverlay(obs,seqIdx)
            % takes in single observation and seqId and returns the
            % TrackMovie object to call tm.viewSequence on
            im_raw = LifetimeAnalysis.getNormedIgEChannel(obs, seqIdx);
            tracks = LifetimeAnalysis.shiftTracksPixels(obs,seqIdx);
            frameBounds = [min(cellfun(@(T) min(T(:,end)),tracks)), max(cellfun(@(T) max(T(:,end)),tracks))];
            im = im_raw(:,:,frameBounds(1):frameBounds(2));
            tm = TrackMovie(tracks,im);
            tm.setTrackColorMethod('TrackLength',logSpaceColorMap(@jet));
        end
        

        
        function ts = getTracksSequence(obs,seqIdx)
            % returns the track for a single observation and seqIdx 
            % 0 = pre file
            % otherwise use post files
            if seqIdx == 0
                ts = obs.preTracks;
            else
                ts = obs.seqTracks{seqIdx};
            end
        end
        
        function rptFileName = getRPTSequenceFileName(obs,seqIdx)
            % returns the rpt pathname for a single observation and seqIdx
            % 0 = pre file
            % otherwise use post files
            if seqIdx == 0
                rptFileName = obs.preRPT;
            else
                rptFileName = obs.seqRPT{seqIdx};
            end
        end
        
        function gmm = makeGMMfromTracks(groupTs, offset, max_intensity) %static function
            if ~isempty(max_intensity)
                groupTs = cellmap(@(Ts) Ts(LifetimeAnalysis.trackFilterIntensity(Ts, max_intensity)), groupTs);
            end
            Ds = cellmap(@(Ts) LifetimeAnalysis.trackDurationFrames(Ts) - offset, groupTs);
            Ds = cellmap(@(d) d(d>=1), Ds);
            gmm = GeometricMixtureModelMLE(Ds);
            gmm.offset = offset;
        end
    end
    methods (Access=protected)
        function checkFilters(obj, filters)
            filters=makecell(filters);
            assert(all(cellfun(@(f) numel(f)==obj.NTotalObs, filters)));
%             if numel(filters)>1
%                 for n=1:obj.Nobs
%                     assert(sum(cellfun(@(f) f(n), filters))<=1); % Check at most one filter is on
%                 end
%             end
        end

        function reloadTracks(obj)
            for ctype_idx = 1:numel(obj.CellTypes)
                cellType = obj.CellTypes{ctype_idx};
                for syktype_idx = 1:numel(obj.SykTypes)
                    sykType = obj.SykTypes{syktype_idx};
                    dirpath = fullfile(obj.dataDir,cellType,sykType,'RPT');
                    if ~exist(dirpath,'dir')
                        fprintf('Unable to find RPT directory for CellType:%s SykType:%s\n  ---> Bad Path: "%s"\n',cellType,sykType,dirpath);
                    end
                    rpt_files = Pickle.listExistingFileNames(dirpath, obj.SykRPTPattern);
                    obsIdx = obj.processRPTFiles(cellType, sykType, rpt_files);
                    ige_rpt_files = Pickle.listExistingFileNames(dirpath, obj.IGERPTPattern);
                    obj.processRPTFiles(cellType, sykType, ige_rpt_files, obsIdx);
                end
            end
        end

        function obsIdx = processRPTFiles(obj, cellType, sykType, rpt_files, obsIdx)
            % obseDateIdx a struct mapping dateStr to indexes ino the global observations struct array
            % note: right now only 1 pre file for a dataset is created. 
                    % if this is updated need to change IntensityChangePoint
            if nargin<5
                obsIdx=[];
            end
            Nfiles = numel(rpt_files);
            % sort the files so the 'pre' files go last and we already know about thier observation sequence once we get to them
            PRERegExp = '[Pp][Rr][Ee]';
            PREOnlyRegExp = [PRERegExp '[oO][nN][lL][yY]'];
            pre_file_mask = cellfun(@(m) ~isempty(m), regexp(rpt_files,PRERegExp));
            rpt_files = [rpt_files(~pre_file_mask),rpt_files(pre_file_mask)];
            for n=1:Nfiles
                rpt = RPT(rpt_files{n});
                if isempty(obj.frameT)
                    obj.frameT = rpt.data.frameT;
                end
                fprintf('[%i/%i] %s/%s Processing file: %s\n',n,Nfiles,cellType,sykType,rpt_files{n});
                [~, basename,~] = fileparts(rpt_files{n});
                parts = strsplit(basename,'_');
                if numel(parts)==2
                    channelName = parts{2};
                    baseParts = strsplit(parts{1},'-');
                    cellName = baseParts{1};
                    baseParts = baseParts(2:end);
                else
                    channelName = parts{3};
                    cellName = parts{1};
                    baseParts = strsplit(parts{2},'-');
                end
                if numel(baseParts)~=7
                    error('LifetimeAnalysis:processRPTFiles','Got incorrect file format: "%s"',basename);
                end
                dStr = strjoin(baseParts(2:7),'-');
                dTime = datetime(cellfun(@str2double,baseParts(2:7)));                
                if ~isempty(regexp(baseParts{1},PREOnlyRegExp,'once'))
                    % This is a .pre file that shouldn't have matching activation observations
                    S = strsplit(baseParts{1},'#');                        
                    conc = 0; %delete -- sscanf(S{1},'%ing'); 
                    seqIdx = sscanf(S{2},'%i');
                    if isempty(seqIdx)
                        seqIdx = 1;
                    end
                    if ~isempty(obsIdx)
                        obsDateIdx = find( [obj.obs(obsIdx).dateTime] == dTime, 1, 'first');
                    end
                    if ~isempty(obsIdx) && ~isempty(obsDateIdx)
                        % This one already exits
                        switch channelName
                            case 'SYK' 
                                obj.obs(obsIdx(obsDateIdx)).seqTracks{seqIdx} = rpt.getTracks();
                                obj.obs(obsIdx(obsDateIdx)).seqRPT{seqIdx} = rpt_files{n};
                            case 'IGE'
                                obj.obs(obsIdx(obsDateIdx)).seqIGETracks{seqIdx} = rpt.getTracks();
                                obj.obs(obsIdx(obsDateIdx)).seqIGERPT{seqIdx} = rpt_files{n};
                        end
                    else
                        O.v  = dStr;
                        O.dateTime = dTime;
                        O.preDate = dTime;
                        O.sykType = sykType;
                        O.cellType = cellType;
                        O.conc = conc;
                        O.cellName = cellName;
                        O.seqTracks=[];
                        O.preRPT=[];
                        O.preTracks=[];
                        O.preIGETracks=[];
                        O.preIGERPT=[];
%                         O.seqTracks = [];
                        O.seqTracks{seqIdx} = rpt.getTracks();
                        O.seqIGETracks = [];
                        O.seqRPT{seqIdx} = rpt_files{n};
                        O.seqIGERPT = [];
                        obj.obs = [obj.obs, O];
                        obsIdx(numel(obsIdx)+1) = obj.NTotalObs; %#ok<AGROW>
                    end
                    
                elseif ~isempty(regexp(baseParts{1},PRERegExp,'once'))
                    % This is a .pre file, find its nearest matching observation
                    [dTs, sortIdx] = sort([obj.obs(obsIdx).dateTime]);
                    nextDTi = find( dTime<dTs , 1, 'first');
                    if isempty(nextDTi);
                         error('LifetimeAnalysis:processRPTFiles','Could not associate pre file: "%s" date:%s',basename,dStr);
                    end
                    O = obj.obs(obsIdx(sortIdx(nextDTi)));
                    assert(strcmp(O.cellName, cellName)); % Check the names are the same
                    obj.obs(obsIdx(sortIdx(nextDTi))).preDate = dTime;
                    switch channelName
                        case 'SYK'
                            obj.obs(obsIdx(sortIdx(nextDTi))).preTracks = rpt.getTracks();
                            obj.obs(obsIdx(sortIdx(nextDTi))).preRPT = rpt_files{n};
                        case 'IGE'
                            obj.obs(obsIdx(sortIdx(nextDTi))).preIGETracks = rpt.getTracks();
                            obj.obs(obsIdx(sortIdx(nextDTi))).preIGERPT = rpt_files{n};
                    end
                    
                else
                    % This is an observation sequence file, find it by date or create it if necessary
                    S = strsplit(baseParts{1},'#');                        
                    conc = sscanf(S{1},'%ing');
                    seqIdx = sscanf(S{2},'%u'); %make sure to read in as base 10 here (%i is base 8 by default (with '00' before the number)
                    if isempty(seqIdx)
                        error('LifetimeAnalysis:processRPTFiles','Unable to read sequence file "%s"',basename);
                    end
                    if ~isempty(obsIdx)
                        obsDateIdx = find( [obj.obs(obsIdx).dateTime] == dTime, 1, 'first');
                    end
                    if ~isempty(obsIdx) && ~isempty(obsDateIdx)
                        % This one already exits
                        switch channelName
                            case 'SYK' 
                                obj.obs(obsIdx(obsDateIdx)).seqTracks{seqIdx} = rpt.getTracks();
                                obj.obs(obsIdx(obsDateIdx)).seqRPT{seqIdx} = rpt_files{n};
                            case 'IGE'
                                obj.obs(obsIdx(obsDateIdx)).seqIGETracks{seqIdx} = rpt.getTracks();
                                obj.obs(obsIdx(obsDateIdx)).seqIGERPT{seqIdx} = rpt_files{n};
                        end
                    else
                        O.v  = dStr;
                        O.dateTime = dTime;
                        O.sykType = sykType;
                        O.cellType = cellType;
                        O.conc = conc;
                        O.cellName = cellName;
                        O.seqTracks=[];
                        O.preDate =[];
                        O.preRPT=[];
                        O.preTracks=[];
                        O.preIGETracks=[];
                        O.preIGERPT=[];
%                         O.seqTracks = [];
                        O.seqTracks{seqIdx} = rpt.getTracks();
                        O.seqIGETracks = [];
                        O.seqRPT{seqIdx} = rpt_files{n};
                        O.seqIGERPT = [];
                        obj.obs = [obj.obs, O];
                        obsIdx(numel(obsIdx)+1) = obj.NTotalObs; %#ok<AGROW>
                    end
                end
                
            end
        end

    end
end




