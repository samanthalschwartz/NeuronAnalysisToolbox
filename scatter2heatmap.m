function [ga,N] = scatter2heatmap(x,y,bins)
if nargin<3
    bins = 100;
end
% Bin the data:
if ~isrow(x)
    x = x';
end
if ~isrow(y)
    y = y';
end
pts = linspace(min([x,y]), max([x,y]), bins);
N = histcounts2(y(:), x(:), pts, pts);

% Plot heatmap:
figure
imagesc(pts, pts, N);
axis equal;
ga = gca;
set(ga, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
colormap jet
end