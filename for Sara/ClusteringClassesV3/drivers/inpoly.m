clear all
close all

x = rand(10, 1);
y = rand(10, 1);

[k, A] = convhull(x, y);
x_hull = x(k);
y_hull = y(k);

x_r = rand(100, 1);
y_r = rand(100, 1);

in = inpolygon(x_r, y_r, x_hull, y_hull);

figure();
hold on
plot(x, y, 'k.', 'MarkerSize', 20); % original points
plot(x_hull, y_hull, 'm-');         % convex hull
plot(x_r(in),  y_r(in),  'b.');     % inside points
plot(x_r(~in), y_r(~in), 'r.');     % outside points
title(sprintf('convhull area = %f', A));
hold off

shp = alphaShape(x, y, 'HoleThreshold', 1);
a1 = criticalAlpha(shp, 'one-region');
shp.Alpha = a1;

[k, v] = boundaryFacets(shp);
k = [k(:, 1); k(1, 1)];
x_hull = v(k, 1);
y_hull = v(k, 2);
in = inShape(shp, x_r, y_r);
A = area(shp);

figure();
hold on
plot(shp);
plot(x_hull, y_hull, 'm-', 'LineWidth', 2);
plot(x, y, 'k.', 'MarkerSize', 20); % original points
plot(x_r(in),  y_r(in),  'b.');     % inside points
plot(x_r(~in), y_r(~in), 'r.');     % outside points
title(sprintf('alphaShape area = %f', A));
hold off

% s = 0   convex hull
% s = 1   alpha shape
s = 0.5;
k = boundary(x, y, s);
x_hull = x(k);
y_hull = y(k);
in = inpolygon(x_r, y_r, x_hull, y_hull);
A = polyarea(x_hull, y_hull);

figure();
hold on
plot(x_hull, y_hull, 'm-');
plot(x, y, 'k.', 'MarkerSize', 20); % original points
plot(x_r(in),  y_r(in),  'b.');     % inside points
plot(x_r(~in), y_r(~in), 'r.');     % outside points
title(sprintf('boundary area = %f', A));
hold off
