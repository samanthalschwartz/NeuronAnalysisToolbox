
% @matrix test
% @matrix MAX_test1
% @OUTPUT net.imagej.Dataset sinkDist
geoframe = bdilation(logical(MAX_test1));
sinkframe = bdilation(logical(test));
sinkDist = bwdistgeodesic(logical(geoframe),logical(sinkframe),'quasi-euclidean');
