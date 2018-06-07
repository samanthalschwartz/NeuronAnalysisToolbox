% make cargominframe file for all AshleyFiles
% alldatadir = {'Z:\Ashley\For Sam\051318 NL1 insertion',...
%     'Z:\Ashley\For Sam\051318 GluA1 insertion',...
%     'Z:\Ashley\For Sam\050118 NL1 insertion'};

alldatadir = {'E:\Zapalog project\GluA1release_REs\060118\merges\slip1'};
alldatadir = {'E:\Zapalog project\DHFR-GFP-GluA1\051718\merges\slip1'};
for dd = 1:numel(alldatadir)
    datadir = alldatadir{dd};
    savedir = fullfile(datadir,'minFrameImages');
    files = dir(fullfile(datadir,'*.mat'));
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    for ff = 1:numel(files)
        thisfile = load(fullfile(datadir,files(ff).name));
        savename = files(ff).name(1:end-4);
        h = thisfile.aa.plot_cargo_minFrame();
        saveas(h,fullfile(savedir,savename),'fig');
        saveas(h,fullfile(savedir,savename),'png');
        close(h);
    end
end