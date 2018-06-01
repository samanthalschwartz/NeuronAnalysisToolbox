function J = green(m)
%Red    Variant of HSV
%   JET(M), a variant of HSV(M), is an M-by-3 matrix containing
%   the default colormap used by CONTOUR, SURF and PCOLOR.
%   The colors begin with dark blue, range through shades of
%   blue, cyan, green, yellow and red, and end with dark red.
%   JET, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   See also HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.

%   Copyright 1984-2004 The MathWorks, Inc.

if nargin < 1
   m = size(get(gcf,'colormap'),1);
end
minval = 1/(m-1);
R = zeros(m,1);
G = [0:minval:1]';
B = zeros(m,1);

J = [R,G,B];

end