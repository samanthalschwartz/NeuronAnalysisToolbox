load('Z:\Sam\MJK_zapERtrap_for_sam\AMB_previous\GluA1release_REs\060118\merges\slip2\1_merge_AshleyFile.mat')
% load('Z:\Sam\MJK_zapERtrap_for_sam\AMB_previous\050118 NL1 insertion\cell4_AshleyFile.mat')
[distmask,h] = aa.makeDistanceMask;
%%
dmask =(distmask>(100/aa.pxsize)&distmask<(200/aa.pxsize))&distmask~=Inf;
% dmask(isnan(dmask)) = Inf;
displayim = distmask*dmask;
jetblack = jet(255);
jetblack(end,:) = [0 0 0];
jetblack(1,:) = [0 0 0];
h = dipshow(displayim,jetblack);
maxdist = max(distmask(distmask~=Inf));
dipmapping(h,[0 max(distmask(distmask~=Inf))]);

%%
% load('Z:\Sam\MJK_zapERtrap_for_sam\AMB_previous\050118 NL1 insertion\cell4_AshleyFile\MinPaths.mat');
% filename = 'Z:\Sam\MJK_zapERtrap_for_sam\AMB_previous\050118 NL1 insertion\cell4_AshleyFile\MinPaths_movie';

load('Z:\Sam\MJK_zapERtrap_for_sam\AMB_previous\GluA1release_REs\060118\merges\slip2\1_merge_AshleyFile\MinPaths.mat')
filename = 'Z:\Sam\MJK_zapERtrap_for_sam\AMB_previous\GluA1release_REs\060118\merges\slip2\1_merge_AshleyFile\MinPaths_movie';

h = dipshow(distMovie);
options.framerate = 10; %default
options.slices = [1:71]; %default. also could do ex. [0:20]
writeDipImageMovie(h,filename,options)


%% TFR global: this is for normalizing to intensity throughout entire cell at max times
% filename = {'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_global_allMs.mat',...
%     'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_TfR_Global_allMs.mat',...
%     'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_GluA1_Global_allMs.mat',...
%     'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_localSoma_allMs'};
filename = {'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_TfR_Global_allMs.mat'};
titlenames = {'Mean Intensity <50 \mum from Soma','Mean Intensity 50-100 \mum from Soma',...
    'Mean Intensity 100-200 \mum from Soma'};
for ff = 1:numel(filename)
    load(filename{ff})
    cols = hsv(numel(allMs));
    f = figure; fa = gca; hold on;
    [FILEPATH,name,EXT]  = fileparts(filename{ff});
    imtitlename = strrep(name,'_','-');
    for mm = [1:2,4:numel(allMs)]
        values = allMs{mm}.intensity;
        values(values<0) = 0;
        sorted = single(sort(values,2));
        totalcellintensity = sum(mean(sorted(:,end-1:end),2));
        for ss = 1:3
            subplot(3,1,ss); hold on;
            results = values(ss,:)./totalcellintensity;
            totalframes = size(results,2);
            baseframes = 1:6;
            baseframetime = 1;
            postframes = 7:totalframes;
            postframetime = 2;
            baselinevals = flip(-((1:6)-1)*1);
            postvals = [1:size(postframes,2)]*postframetime;
            timevals = [baselinevals postvals];
            plot(timevals,results','Color',cols(mm,:)); hold on;
            yL = get(gca,'YLim');
            line([timevals(baseframes(end)) timevals(baseframes(end))],yL,'Color','k');
            title(titlenames{ss});
            xlabel('Time (min)')
        end
    end
    hout=suptitle(['Mean Delivery Density for Cargo: ' imtitlename ]);
end
%% NL1 global: this is for normalizing to intensity throughout entire cell at max times
% filename = {'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_global_allMs.mat',...
%     'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_TfR_Global_allMs.mat',...
%     'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_GluA1_Global_allMs.mat',...
%     'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_localSoma_allMs'};
filename = {'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_global_allMs.mat'};
titlenames = {'Mean Intensity <50 \mum from Soma','Mean Intensity 50-100 \mum from Soma',...
    'Mean Intensity 100-200 \mum from Soma'};

baseframes{1:9} = 1:11;
baseframetime{1:9} = 1;
postframes{1:9} = 12:totalframes;
postframetime{1:9} = 2;
baseframes{10:20} = 1:6;
for ff = 1:numel(filename)
    load(filename{ff})
    cols = hsv(numel(allMs));
    f = figure; fa = gca; hold on;
    [FILEPATH,name,EXT]  = fileparts(filename{ff});
    imtitlename = strrep(name,'_','-');
    for mm = [1:2,4:numel(allMs)]
        values = allMs{mm}.intensity;
        values(values<0) = 0;
        sorted = single(sort(values,2));
        totalcellintensity = sum(mean(sorted(:,end-2:end),2));
        for ss = 1:3
            subplot(3,1,ss); hold on;
            results = values(ss,:)./totalcellintensity;
            totalframes = size(results,2);
            baseframes = 1:6;
            baseframetime = 1;
            postframes = 7:totalframes;
            postframetime = 2;
            
            baselinevals = flip(-((baseframes)-1)*1);
            postvals = [1:size(postframes,2)]*postframetime;
            timevals = [baselinevals postvals];
            plot(timevals,results','Color',cols(mm,:)); hold on;
            yL = get(gca,'YLim');
            line([timevals(baseframes(end)) timevals(baseframes(end))],yL,'Color','k');
            title(titlenames{ss});
            xlabel('Time (min)')
        end
    end
    hout=suptitle(['Mean Delivery Density for Cargo: ' imtitlename ]);
end
%% GluA1 global: this is for normalizing to intensity throughout entire cell at max times
% filename = {'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_global_allMs.mat',...
%     'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_TfR_Global_allMs.mat',...
%     'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_GluA1_Global_allMs.mat',...
%     'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_localSoma_allMs'};
filename = {'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_GluA1_Global_allMs.mat'};
titlenames = {'Mean Intensity <50 \mum from Soma','Mean Intensity 50-100 \mum from Soma',...
    'Mean Intensity 100-200 \mum from Soma'};
for ff = 1:numel(filename)
    load(filename{ff})
    cols = hsv(numel(allMs));
    f = figure; fa = gca; hold on;
    [FILEPATH,name,EXT]  = fileparts(filename{ff});
    imtitlename = strrep(name,'_','-');
    for mm = [1:2,4:numel(allMs)]
        values = allMs{mm}.intensity;
        values(values<0) = 0;
        sorted = single(sort(values,2));
        totalcellintensity = sum(mean(sorted(:,end-2:end),2));
        for ss = 1:3
            subplot(3,1,ss); hold on;
            results = values(ss,:)./totalcellintensity;
            totalframes = size(results,2);
            baseframes = 1:6;
            baseframetime = 1;
            postframes = 7:totalframes;
            postframetime = 2;
            
            baselinevals = flip(-((1:6)-1)*1);
            postvals = [1:size(postframes,2)]*postframetime;
            timevals = [baselinevals postvals];
            plot(timevals,results','Color',cols(mm,:)); hold on;
            yL = get(gca,'YLim');
            line([timevals(baseframes(end)) timevals(baseframes(end))],yL,'Color','k');
            title(titlenames{ss});
            xlabel('Time (min)')
        end
    end
    hout=suptitle(['Mean Delivery Density for Cargo: ' imtitlename ]);
end

%% Make and save excel files for each of the datasets to manipulate
filename = {'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_localSoma_allMs'};
% filename = {'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_global_allMs.mat',...
%     'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_TfR_Global_allMs.mat',...
%     'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_GluA1_Global_allMs.mat',...
%     'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_localSoma_allMs'};
savedir = 'Z:\Sam\MJK_zapERtrap_for_sam\Presentations\excel results files';
savestrs = {'MI_50','MI_100','MI_200'};
ff =1;
load(filename{ff})
[FILEPATH,name,EXT]  = fileparts(filename{ff});

for ss = 1:3
    results = [];
    savename = fullfile(savedir,[name '_' savestrs{ss}]);
    for mm = [1:numel(allMs)]
        values = allMs{mm}.intensity;
        values(values<0) = 0;
        sorted = single(sort(values,2));
        totalcellintensity = sum(mean(sorted(:,end-2:end),2));
        
        if mm>1 & size(results,2) < size(values(ss,:),2)
            results = cat(2,results,zeros(size(results)));
        end
        results(mm,1:size(values(ss,:),2))= values(ss,:)./totalcellintensity;
    end
    results(results == 0) = NaN;
    xlswrite(savename,results);
end
