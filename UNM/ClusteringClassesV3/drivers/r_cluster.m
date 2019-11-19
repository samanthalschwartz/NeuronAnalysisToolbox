E_range = 50 : 10 : 60;
pixel2nm = 1.448;

algorithm = 'Hierarchical';
minPts = 3;

n_range = length(E_range);

indir  = uigetdir('.', 'Input Directory');
outdir = uigetdir('.', 'Output Directory');
D = dir(fullfile(indir, '*.txt'));
n_files = length(D);
if n_files == 0
   error('No coordinate files found in\n%s\n', indir);
end

for i = 1 : n_files
   filename = fullfile(indir, D(i).name);
   [inpath, base_file, ext] = fileparts(filename);
   fprintf('Opening %s ...\n', filename);
   [x, y] = textread(filename, '%*f %f %f %*f', 'headerlines', 1);
   XY = [x, y] .* pixel2nm;

   n_C = zeros(n_range, 1);

   c = Clustering();
   j = 0;
   for E = E_range
      [nC, C, centers, ptsI] = c.cluster(algorithm, XY, E, minPts);
      fprintf('E = %d: clusters = %d', E, nC);

      results = c.clusterStats(XY, C, centers);
      j = j + 1;
      n_C(j) = results.nC;

      clusterFig = c.plotClusters(XY, C, centers, ptsI,                  ...
                                  sprintf('%s: eps = %d', base_file, E), ...
                                  0, 'O1');
      %showm(clusterFig);
      outfile = fullfile(outdir, sprintf('%s_cluster_E=%d', base_file, E));
      saveas(clusterFig, outfile, 'fig');
      print('-dpdf', clusterFig, outfile);
   end
end
