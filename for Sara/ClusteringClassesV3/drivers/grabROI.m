DataDirs = {...
%'F:\sr test data\Results3 - Copy',...
'.', ...
};

RT = ROITools();
for i = 1 : numel(DataDirs)
    DataDir = DataDirs{i};
    Files = dir(fullfile(DataDir, '*.mat'));
    FileName = fullfile(DataDir, Files.name);
    [RoI, Sigma_Reg] = RT.getROI(FileName);
    SaveFile = regexprep(FileName, '\.mat$', '_ROI.mat')
    save(SaveFile, 'RoI', 'Sigma_Reg');
end
