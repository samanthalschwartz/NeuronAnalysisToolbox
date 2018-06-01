function perim = maskperim(mask,dist)
if nargin<2
    dist=1;
end
    maskdt = dt(mask);
    perim = maskdt==dist;
end
    