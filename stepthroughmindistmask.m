function track = stepthroughmindistmask(mindistmask,sinkdistancetransform,startpixel,mindistance)
% startpixel = [2,5] matlab matrix indexing
maxtracklength = floor(mindistance) + 1;
track = zeros(maxtracklength,2);
tracklength = 1;
track(1,:) = startpixel;
dist2sink = mindistance;
neighborstencil = [-1,-1;-1,0;-1,1;0,-1;0,1;1,-1;1,0;1,1];
curloc = startpixel;
while dist2sink>0
    neighbors = neighborstencil + repmat(curloc,size(neighborstencil,1),1);
    neighbordistance = sinkdistancetransform(neighbors);
    idxs = find(neighbordistance == min(neighbordistance));
    if numel(idxs)>1 %if more than 1 possible minimum distance neighbor, chose a neighbor randomly
        idx = idxs(randi(numel(idxs),1));
    else
        idx = idxs;
    end
    curloc = neighbors(idx);
    tracklength = tracklength + 1;
    track(tracklength,:) = curloc;
    dist2sink = neighbordistance(idx);
end

end