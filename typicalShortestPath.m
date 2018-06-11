function track = typicalShortestPath(sinkdistancetransform,startpixel,mindistance)
% if seed input is a mask: measure and take center of mass position


% startpixel = [2,5] matlab matrix indexing
maxtracklength = floor(mindistance) + 1;
track = nan(maxtracklength,2);
tracklength = 1;
track(1,:) = startpixel;
dist2sink = mindistance;
neighborstencil = [-1,-1;-1,0;-1,1;0,-1;0,1;1,-1;1,0;1,1];
coststencil = [sqrt(2),0,sqrt(2),0,0,sqrt(2),0,sqrt(2)];
curloc = startpixel;
dip_sinkdistancetransform = dip_image(sinkdistancetransform);
% make 1 pixel of nan around the edge

while dist2sink>0
    neighbors = neighborstencil + repmat(curloc,size(neighborstencil,1),1);
    neighbors(:,:) = max(neighbors(:,:),0);
    neighbors(:,:) = min(neighbors(:,:),size(dip_sinkdistancetransform,1)-1);
    neighbordistance = [];
    for ii = 1:size(neighbors,1)
    neighbordistance(ii) = dip_sinkdistancetransform(neighbors(ii,1),neighbors(ii,2));
    end
    neighbordistance = neighbordistance + coststencil;
    idxs = find(neighbordistance == min(neighbordistance));
    if numel(idxs)>1 %if more than 1 possible minimum distance neighbor, chose a neighbor randomly
        idx = idxs(randi(numel(idxs),1));
    else
        idx = idxs;
    end
    curloc = neighbors(idx,:);
    tracklength = tracklength + 1;
    track(tracklength,:) = curloc;
    dist2sink = neighbordistance(idx);
end
track(isnan(track(:,1)),:) = [];
% figure; hold on;
% for ii = 1:size(track,1)
%    scatter(track(ii,:),'*');   
% end
end