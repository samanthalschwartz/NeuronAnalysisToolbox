function [X1, Y1, X2, Y2] = getData(data_dir)

   file1 = dir(fullfile(data_dir, 'label1', '*_ResultsStruct.mat'));
   file2 = dir(fullfile(data_dir, 'label2', '*_ResultsStruct.mat'));
   SMR1 = load(fullfile(data_dir, 'label1', file1.name));
   SMR2 = load(fullfile(data_dir, 'label2', file2.name));

   X1 = SMR1.SMR.X;
   Y1 = SMR1.SMR.Y;
   X2 = SMR2.SMR.X;
   Y2 = SMR2.SMR.Y;

end
