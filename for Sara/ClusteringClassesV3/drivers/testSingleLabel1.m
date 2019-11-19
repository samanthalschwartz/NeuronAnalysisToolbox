clear all
close all

Saving = false;
%Saving = true;

%Plotting = true;
Plotting = false;

pixel2nm = 16000/150;
ShrinkFactor = 0.5;

SD = SimulateDomains();
SD.Printing = false;   % Print statistics

SRc = SRcluster();
SRc.ProduceLegend = false;
SRc.Timing = false;
SRc.E = SD.Sigma_loc;
%SRc.Method = 'keith_HC';

%n_studies = 20;
n_studies = 1;
n_trials = 100;

study = 4;

switch study
case 1
   M = 100;
case 2
   M = 100;
case 3
   M = 100;
case 4
   M = 10;
end

mean_new = zeros(M, n_studies);
stdv_new = zeros(M, n_studies);

if Saving
   fid = fopen('../studies/summary.txt', 'w');
end

for m = 1 : M

   switch study
   case 1
      SD.N_observ = 10;      % observations per molecule
   case 2
      SD.N_observ = 10*m;
   case 3
      SD.N_observ = 10;
   case 4
      SD.N_observ = 10;
   end

   switch study
   case 1
      SD.Sigma_loc = 20;     % localization error in each dimension (nm)
   case 2
      SD.Sigma_loc = 20;
   case 3
      SD.Sigma_loc = 20*m;
   case 4
      SD.Sigma_loc = 20;
   end

   switch study
   case 1
      Sigma_Reg = [10, 10] * (m - 1);
   case 2
      Sigma_Reg = [10, 10];
   case 3
      Sigma_Reg = [10, 10];
   case 4
      Sigma_Reg = [10, 10];
   end

   old = zeros(n_studies, n_trials);
   new = zeros(n_studies, n_trials);
   count = zeros(n_studies, n_trials, 3);

   n_points = 0;
   x = -500;
   dy = 0;
   if Plotting
      h = figure();
      hold on
   end
   for i = 0 : n_studies - 1
      x = x + 1000;
      dy = dy + 10;
      for j = 0 : n_trials - 1 
         y1 = -dy/2 + 500 + 1000*j;
         %y2 = y1 + dy;
         [pts1, sigma1] = SD.generateObservations([x, y1]);
         if study == 4
            y1 = y1 .* m;
            pts1 = pts1 .* m;
         end
         %[pts2, sigma2] = SD.generateObservations([x, y2]);
         %xy_old = [pts1; pts2];
         %std_old = [sigma1; sigma2];
         %center = [x, y1; x, y2];
         xy_old = [pts1];
         std_old = [sigma1];
         center = [x, y1];
         n_pts = size(xy_old, 1);
         n_points = n_points + n_pts;

         old(i + 1, j + 1) = size(xy_old, 1);
         %std_old = SD.Sigma_loc .* ones(size(xy_old));

         [xy_new, std_new, combined] = ...
            SRc.clusterSR(xy_old, std_old, Sigma_Reg);
         new(i + 1, j + 1) = size(xy_new, 1);

         if Plotting
            plot(xy_old(:, 1), xy_old(:, 2), 'c.', 'MarkerSize', 10);
            plot(xy_new(:, 1), xy_new(:, 2), 'b.', 'MarkerSize', 10);
            plot(center(:, 1), center(:, 2), 'ro', 'MarkerSize', 10);
         end

         total_clustered = 0;
         for l = 1 : length(combined)
            c = combined{l};
            n_clustered = length(c);
            total_clustered = total_clustered + n_clustered;
            if n_clustered == 1
               % This never really happens because of the nature of combined.
               %count(i + 1, j + 1, 1) = count(i + 1, j + 1, 1) + 1;
               fprintf('n_clustered == 1!\n');
            elseif n_clustered == 2
               count(i + 1, j + 1, 2) = count(i + 1, j + 1, 2) + 1;
               if Plotting
                  plot(xy_old(c, 1), xy_old(c, 2), 'k-');
               end
            else
               count(i + 1, j + 1, 3) = count(i + 1, j + 1, 3) + 1;
               if Plotting
                  %k = convhull(xy_old(c, 1), xy_old(c, 2));
                  k = boundary(xy_old(c, 1), xy_old(c, 2), ShrinkFactor);
                  plot(xy_old(c(k), 1), xy_old(c(k), 2), 'k-');
               end
            end
         end
         count(i + 1, j + 1, 1) = n_pts - total_clustered;
      end
   end
   if Plotting
      axis equal
      %xlabel('dy = 10, 20, ... (nm)');
      ylabel('independent trials');
      hold off
      if Saving
         saveas(h, sprintf('../studies/study_%d', m), 'fig');
      end
   end

   fprintf( ...
'N_points = %d, N_trials = %d, N_observ = %d, Sigma_Reg = [%d, %d]\n', ...
           n_points, n_trials, SD.N_observ, Sigma_Reg);
   if Saving
   fprintf(fid, ...
'N_points = %d, N_trials = %d, N_observ = %d, Sigma_Reg = [%d, %d]\n', ...
           n_points, n_trials, SD.N_observ, Sigma_Reg);
   end

   mean_new(m, :) = mean(new, 2)';
   stdv_new(m, :) = std(new');
   for i = 1 : 3
      mean_count(m, :, i) = mean(count(:, :, i), 2);
   end

   dy = 0;
   for i = 1 : n_studies
      dy = dy + 10;
      fprintf('dy = %3d: %.3f +- %.3f => (1) %.3f, (2) %.3f, (>2) %.3f\n', ...
              dy, mean_new(m, i), stdv_new(m, i), mean_count(m, i, :));
      if Saving
      fprintf(fid,                                                         ...
              'dy = %3d: %.3f +- %.3f => (1) %.3f, (2) %.3f, (>2) %.3f\n', ...
              dy, mean_new(m, i), stdv_new(m, i), mean_count(m, i, :));
      end
   end

end

if Saving
   fclose(fid);
end

color = ['k', 'g', 'b', 'm', 'c'];
studies = 10 : 10 : 10*n_studies;
leg = cell(M, 1);

h = figure();
hold on
for m = 1 : M
   switch study
   case 1
      leg{m} = sprintf('Sigma\\_Reg = %d', 10*(m - 1));
   case 2
      leg{m} = sprintf('mean(N\\_obs) = %d', 10*m);
   case 3
      leg{m} = sprintf('sigma\\_loc = %d', 20*m);
   case 4
      leg{m} = sprintf('multiplier = %d', m);
   end
   clr = color(SD.nMODm(m, 5));
   %plot(studies, mean_new(m, :), [clr, '-'], 'LineWidth', 3);
   %plot(studies, mean_new(m, :) - stdv_new(m, :), 'r--');
   %plot(studies, mean_new(m, :) + stdv_new(m, :), 'r--');
   plot(studies, mean_new(m, :), [clr, 'o'], 'MarkerSize', 10);
   if study ~= 4
   plot([studies, studies], ...
        [mean_new(m, :), mean_new(m, :) - stdv_new(m, :)], [clr, '--']);
   plot([studies, studies], ...
        [mean_new(m, :), mean_new(m, :) + stdv_new(m, :)], [clr, '--']);
   end
end
xlabel('dy (nm)');
ylabel('# of collapsed points');
legend(leg, 'Location', 'NorthEast');
hold off
if Saving
   saveas(h, '../studies/summary', 'fig');
end
