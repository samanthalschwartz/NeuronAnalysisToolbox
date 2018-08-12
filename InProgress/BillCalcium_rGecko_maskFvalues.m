% script to get out the intensity over baseline for each event max
% intensity
folderbase = 'E:\Bill\Geph_IB-Cry2Olig-GFP plus blue light\calcium imaging\0 Mg QCTs and dendritic spikes\8.6.18_ibOLIG_QCTs';
excelfileinfo_Fnorm = 'E:\Bill\Geph_IB-Cry2Olig-GFP plus blue light\calcium imaging\0 Mg QCTs and dendritic spikes\MaxIntensities_REDO.xlsx';

foldernames = {'8.6.18_ibOLIG_QCTs_well1','8.6.18_ibOLIG_QCTs_well2','8.6.18_ibOLIG_QCTs_well3','8.6.18_ibOLIG_QCTs_well4',...
    '8.6.18_ibOLIG_QCTs_well9'};
possiblefilenames = {'baseline_movie.tif','postZ_movie.tif'};
savenames = {'base','postZ'};
for ff = 1:numel(foldernames)
    for pf = 1:numel(possiblefilenames)
        filename1 = fullfile(folderbase,foldernames{ff},possiblefilenames{pf});
        fil = dir(filename1);
        filename = fullfile(folderbase,foldernames{ff},fil(1).name);
        ainfo = load([filename(1:end-4) '_analysisinfo.mat']);
        ca = GeneralAnalysis.loadtiff_1ch(filename);
        ca = ca(:,:,0:end-1);
        Fval = zeros(max(ainfo.lbl),1);
        Fnorm = zeros(max(ainfo.lbl),1);
        % loop through each label
        wb = waitbar(0,'waiting....');
        for ll = 1:max(ainfo.lbl)
            % now calculate sum intensity inside mask
            currmask = ainfo.lbl==ll;
            sumvaltrace = sum(ca,currmask,[1 2]);
            [Y,I] = max(single(sumvaltrace),[],3);
            maxmask = currmask(:,:,I-1);
            imagetrace = sum(ca,repmat(maxmask,1,1,size(ca,3)),[1 2]);
            baseline = median(imagetrace);
            Fval(ll) = Y;
            Fnorm(ll) = Y/baseline;
            waitbar(ll/max(ainfo.lbl),wb);
        end
        close(wb);
        sheetname = [foldernames{ff} '-' savenames{pf}];
        xlswrite(excelfileinfo_Fnorm,[Fval,Fnorm],sheetname);
    end
    disp(['Analyzing file ' foldernames{ff} '...']);
end

