classdef ROITools < handle

% ROITools class written by Michael Wester (10/17/2017) <wester@math.unm.edu>
% The New Mexico Center for the Spatiotemporal Modeling of Cell Signaling
% University of New Mexico Health Sciences Center
% Albuquerque, New Mexico, USA   87131
% Copyright (c) 2017-2018 by Michael J. Wester and Keith A. Lidke

% =============================================================================
properties
% =============================================================================

   pixel2nm = 16000/150;   % Conversion factor from pixels to nm.
   % Box diameters used when clicking the left mouse button (nm).
   ROI_sizes = [3000, 3000];   % nm
   color = ['g', 'm', 'k', 'c', 'r', 'y'];   % Label colors.
   order = 1 : 6;                            % Plotting order.
   msize = 7;                                % Marker size.
   XYvsYX = true;                            % Coordinate order.

% =============================================================================
end % properties

methods
% =============================================================================

function [n_ROIs, RoI, Sigma_Reg] = getROI(obj, src, txt)
% Let the user select ROIs from
%    an SMA_SR or SR_demo data structure,
%    point (x, y) coordinates [provided as N x 2 matrices, or if x_STD and
%       y_STD are appended, then N x 4 matrices],
%    super-resolution image files containing an SMA_SR or SR_demo data
%       structure,
%    or directly from a data structure with fields X and Y.
%
% Inputs:
%    src              cell array of (x, y) coordinate sources, although it is
%                     permissible to omit the cell array for a single source
%    txt              [OPTIONAL] text to label the ROI figure
% Outputs (units in nm):
%    n_ROIs           number of ROIs created
%    RoI              cell array for each ROI of
%                        ROI            [xmin, xmax, ymin, ymax] of ROI
%                        X, Y           (x, y) coordinates of points inside
%                        X_STD, Y_STD   (x, y) localization errors for the
%                                       above
%    Sigma_Reg        cell array (per label) of image registration error

   x_size = obj.ROI_sizes(1);
   y_size = obj.ROI_sizes(2);

   if ~exist('txt', 'var')
      txt = '';
   end

   n_ROIs = 0;
   RoI = [];

   if ~iscell(src)
      [XY{1}, XY_STD{1}, Sigma_Reg{1}] = ROITools.import_XY(src, obj.pixel2nm);
   else
      n_labels = numel(src);
      for i = 1 : n_labels
         [XY{i}, XY_STD{i}, Sigma_Reg{i}] = ...
            ROITools.import_XY(src{i}, obj.pixel2nm);
      end
   end
   [n_ROIs, RoI] = obj.getROI_XY(XY, XY_STD, x_size, y_size, txt);

end

% -----------------------------------------------------------------------------

function [n_ROIs, RoI] = getROI_XY(obj, XY, XY_STD, x_size, y_size, txt)
% Let the user select ROIs from (x, y)-coordinates.
%
% Inputs:
%    XY               cell array of (x, y) coordinates of points labeled in
%                     image (nm); XY{i} is the coordinates for label i
%    XY_STD           cell array of (x, y) standard deviations of the
%                     coordinate values (nm)
%    x_size, y_size   box diameters used when clicking the left mouse button
%                     (nm)
%    txt              text to label the ROI figure
% Outputs:
%    n_ROIs           number of ROIs created
%    RoI              cell array for each ROI of
%                        ROI            [xmin, xmax, ymin, ymax] of ROI
%                        X, Y           (x, y) coordinates of points inside
%                        X_STD, Y_STD   (x, y) localization errors for the
%                                       above
%                     NOTE: X{i}{j} is the x-coordinates of label j in ROI i,
%                     etc.

   if obj.XYvsYX
      ix = 1;
      iy = 2;
   else
      ix = 2;
      iy = 1;
   end

   n_labels = numel(XY);

   for j = 1 : n_labels
      n = size(XY{j}, 1);
      X_STD{j} = NaN(n, 1);
      Y_STD{j} = NaN(n, 1);

      X{j} = XY{j}(:, ix);
      Y{j} = XY{j}(:, iy);
      if ~isempty(XY_STD{j})
         X_STD{j} = XY_STD{j}(:, ix);
         Y_STD{j} = XY_STD{j}(:, iy);
      end
   end

   [n_ROIs, ROI, index_ROI] = obj.get_ROI(X, Y, x_size, y_size, txt);
   RoI = cell(1, n_ROIs);
   for i = 1 : n_ROIs
      fprintf('ROI %d = %.3f %.3f %.3f %.3f\n', i, ROI{i});
      RoI{i}.ROI   = ROI{i};
      for j = 1 : n_labels
         k = index_ROI{i}{j};
         RoI{i}.X{j}     = X{j}(k);
         RoI{i}.Y{j}     = Y{j}(k);
         RoI{i}.X_STD{j} = X_STD{j}(k);
         RoI{i}.Y_STD{j} = Y_STD{j}(k);
      end
   end

end

% -----------------------------------------------------------------------------

function [n_ROIs, ROI, index_ROI] = get_ROI(obj, X, Y, x_size, y_size, txt)
% Use the mouse to select ROIs (regions of interest):
%    left click  chooses the center of a fixed size (x_size x y_size) region
%    right click chooses an adjustable rectangular size region
%    key press:
%       backspace or delete   deletes the previous region
%       anything else         terminates selection
%
% Inputs:
%    X                cell array of the x-coordinates of the labeled points in
%                     the entire image, where X{i} corresponds to the label i
%                     x-coordinates
%    Y                cell array of the y-coordinates of the labeled points in
%                     the entire image, where Y{i} corresponds to the label i
%                     y-coordinates
%    x_size, y_size   box diameters used when clicking the left mouse button
%    txt              text to label the ROI figure
% Outputs:
%    n_ROIs           number of ROIs created
%    ROI              cell array of [xmin, xmax, ymin, ymax] for each ROI
%    index_ROI        cell array of labeled point indices in each ROI, where
%                     index_ROI{i}{j} corresponds to label j in ROI i

   n_labels = numel(X);

   n_ROIs = 0;   ROI = [];   index_ROI = [];

   BS = char(8);   DEL = char(127);

   % selected = 0   terminate selection
   %            1   valid button press
   %            2   ignored or region delete
   h = figure();
   % An idea by Samantha Schwartz
   %set(h, 'Position', [200, 70, 900, 900*(x_size/y_size)]);
   hold on
   for i = 1 : n_labels
      j = obj.order(i);
      plot(X{j}, Y{j}, [obj.color(j), '.'], 'MarkerSize', obj.msize);
   end
   xlabel('x (nm)');
   ylabel('y (nm)');
   title(txt);
   done = false;
   while ~done
      clickval = waitforbuttonpress;
      if clickval == 0   % if a mouse button was pressed ...
         clickType = get(gcf, 'SelectionType');
         fprintf('%s: ', clickType);
         switch clickType
         case 'normal'   % left button press: draw a fixed rectangle
            selected = 1;
            curpt = get(gca, 'CurrentPoint');
            xmin = curpt(1, 1) - x_size/2;
            xmax = curpt(1, 1) + x_size/2;
            ymin = curpt(1, 2) - y_size/2;
            ymax = curpt(1, 2) + y_size/2;
         case 'alt'      % right button press: draw an adjustable rectangle
            selected = 1;
            rect = getrect;
            xmin = rect(1);
            xmax = rect(1) + rect(3);
            ymin = rect(2);
            ymax = rect(2) + rect(4);
         otherwise       % middle button press (extend) and double click (open)
            selected = 2;
         end
      else   % key was pressed
         charChoice = get(gcf, 'CurrentCharacter');
         % If a backspace or a delete, cancel the previous ROI
         if charChoice == BS | charChoice == DEL
            selected = 2;
            if n_ROIs > 0
               delete(r(n_ROIs));
               %delete(p(n_ROIs));
               delete(t(n_ROIs));

               ROI{n_ROIs} = [];
               index_ROI{n_ROIs} = [];

               fprintf('delete ROI %d\n', n_ROIs);
               n_ROIs = n_ROIs - 1;
            end
         else
            selected = 0;
         end
      end
      if selected == 1
         n_ROIs = n_ROIs + 1;
         fprintf('add ROI %d\n', n_ROIs);

         r(n_ROIs) = plot([xmin, xmin, xmax, xmax, xmin], ...
                          [ymin, ymax, ymax, ymin, ymin], ...
                          'r-', 'LineWidth', 3);

         ROI{n_ROIs} = [xmin, xmax, ymin, ymax];
         for i = 1 : n_labels
            index_ROI{n_ROIs}{i} = ...
               xmin <= X{i} & X{i} <= xmax & ymin <= Y{i} & Y{i} <= ymax;
         end

         not_empty = true;
         for i = 1 : n_labels
            not_empty = not_empty && any(index_ROI{n_ROIs}{i} == true);
         end
         if not_empty
            XX = [];   YY = [];
            for i = 1 : n_labels
               XX = [ XX; X{i}(index_ROI{n_ROIs}{i}) ];
               YY = [ YY; Y{i}(index_ROI{n_ROIs}{i}) ];
            end
            %p(n_ROIs) = plot(XX, YY, 'r.');   % color the selected points red
            t(n_ROIs) = text((xmin + xmax)/2, (ymin + ymax)/2, ...
                             int2str(n_ROIs));
            t(n_ROIs).Color = 'black';
            t(n_ROIs).FontWeight = 'bold';
         else
            fprintf('One label has no points in this ROI!  Deleting ...\n');
            delete(r(n_ROIs));

            ROI{n_ROIs} = [];
            index_ROI{n_ROIs} = [];

            fprintf('delete ROI %d\n', n_ROIs);
            n_ROIs = n_ROIs - 1;
         end
      elseif selected == 0
         done = true;
      end
   end
   hold off

   if length(ROI) > n_ROIs
      ROI = ROI(~cellfun('isempty', ROI));
      index_ROI = index_ROI(~cellfun('isempty', index_ROI));
   end

end

% =============================================================================
end % properties

methods(Static)
% =============================================================================

function [XY, XY_STD, Sigma_Reg] = import_XY(src, pixel2nm)
% Import N x 2 (x, y) coordinates and standard deviations from
%    an SMA_SR or SR_demo data structure,
%    point (x, y) coordinates [provided as N x 2 matrices, or if x_STD and
%       y_STD are appended, then N x 4 matrices],
%    super-resolution image files containing an SMA_SR or SR_demo data
%       structure,
%    or directly from a data structure with fields X and Y.
%
% Inputs:
%    src              (x, y) coordinate source
%    pixel2nm         conversion factor from pixels to nm
% Outputs (units in nm):
%    XY               array of (x, y) point coordinates [N x 2]
%    XY_STD           array of (x, y) standard deviations of the coordinate
%                     values [N x 2]
%    Sigma_Reg        image registration error in x and y [1 x 2]

   XY_STD = [];
   Sigma_Reg = [10, 10];   % default registration error

   if strcmp(class(src), 'SMA_SR')
      % src is an SR_demo data structure
      XY = [ double(src.SMR.X) .* pixel2nm, ...
             double(src.SMR.Y) .* pixel2nm ];
      XY_STD = [ double(src.SMR.X_SE) .* pixel2nm, ...
                 double(src.SMR.Y_SE) .* pixel2nm ];
      if isprop(src, 'DriftCorrect_XYShift') && ...
         ~isempty(src.DriftCorrect_XYShift)
         Sigma_Reg = std(src.DriftCorrect_XYShift) .* pixel2nm;
      end
   elseif strcmp(class(src), 'SR_demo')
      % src is an SR_demo data structure
      XY = [ double(src.Results.X) .* pixel2nm, ...
             double(src.Results.Y) .* pixel2nm ];
      XY_STD = [ double(src.Results.X_STD) .* pixel2nm, ...
                 double(src.Results.Y_STD) .* pixel2nm ];
      if isprop(src, 'DriftCorrect_XYShift') && ...
         ~isempty(src.DriftCorrect_XYShift)
         Sigma_Reg = std(src.DriftCorrect_XYShift) .* pixel2nm;
      end

   elseif ismatrix(src) && ~ischar(src) && ~isstruct(src)
      % src is a matrix
      n_cols = size(src, 2);
      if n_cols == 2 || n_cols == 4
         XY = src(:, 1:2) .* pixel2nm;
      else
         error('src matrix has %d rather than 2 or 4 columns!', n_cols);
      end
      if n_cols == 4
         XY_STD = src(:, 3:4) .* pixel2nm;
      end

   elseif ischar(src)
      % src is a filename
      load(src);
      if exist('SMASR', 'var')
         [XY, XY_STD, Sigma_Reg] = ROITools.import_XY(SMASR, pixel2nm);
      elseif exist('SMR', 'var')
         [XY, XY_STD, Sigma_Reg] = ROITools.import_XY(SMR, pixel2nm);
      elseif exist('SR', 'var')
         [XY, XY_STD, Sigma_Reg] = ROITools.import_XY(SR, pixel2nm);
      elseif exist('SRD', 'var')
         [XY, XY_STD, Sigma_Reg] = ROITools.import_XY(SRD, pixel2nm);
      elseif exist('SRtest', 'var')
         [XY, XY_STD, Sigma_Reg] = ROITools.import_XY(SRtest, pixel2nm);
      else
         error('No SMASR, SMR, SR, SRD or SRtest object found in %s!', src);
      end

   elseif isstruct(src)
      % src should be a data structure with fields X, Y and optionally
      % X_SE, Y_SE or X_STD, Y_STD
      if isfield(src, 'X') & isfield(src, 'Y')
         XY = [ double(src.X) .* pixel2nm, double(src.Y) .* pixel2nm ];
         if isfield(src, 'X_SE') & isfield(src, 'Y_SE')
            XY_STD = [ double(src.X_SE) .* pixel2nm, ...
                       double(src.Y_SE) .* pixel2nm ];
         elseif isfield(src, 'X_STD') & isfield(src, 'Y_STD')
            XY_STD = [ double(src.X_STD) .* pixel2nm, ...
                       double(src.Y_STD) .* pixel2nm ];
         end
      else
         error('Fields X, Y not found in src!');
      end

   else
      error('Incorrect argument type for src!');
   end

end

% =============================================================================
end % methods(Static)
% =============================================================================
end % classdef
