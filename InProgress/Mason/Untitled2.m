% testing local thresholding

T = adaptthresh(image,'ForegroundPolarity','bright');
% dipshow(T);
BW = imbinarize(image,T);
dipshow(BW)

test = single(image)-single(T);
T2 = adaptthresh(test,'ForegroundPolarity','bright');
BW2 = imbinarize(test,T2);
dipshow(BW2)

%%

out = bclosing(currim,2,-3,0);
test = bwmorph3(logical(currim),'clean');
test = GeneralAnalysis.bwmorph_timeseries(currim,'bridge',50);
test2 = bclosing(logical(test),3,-3,0);