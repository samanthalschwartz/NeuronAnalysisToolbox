function createROIs(data_dir, results_dir, desc, pixel2nm)
% Create ROIs from the image data in the data_dir.

   [X1, Y1, X2, Y2] = getData(data_dir);
   XY1 = double([X1, Y1] .* pixel2nm);
   XY2 = double([X2, Y2] .* pixel2nm);

   RT = ROITools();
   RT.color = ['g', 'm'];
   RT.ROI_sizes = [3000, 3000];
   txt = regexprep(desc, '_', '\\_');
   [n_ROIs, RoI, ~] = RT.getROI({XY1, XY2}, txt);

   saveas(gcf, fullfile(results_dir, sprintf('%s_ROIs.fig', desc)));
   print(fullfile(results_dir, sprintf('%s_ROIs.png', desc)), '-dpng');
   save(fullfile(results_dir, 'ROIs.mat'), 'n_ROIs', 'RoI');

end
