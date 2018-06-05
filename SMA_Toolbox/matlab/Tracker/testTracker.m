
N = 500;
nFrames = 100;
frameT = 0.05;
originX = 0;
originY = 1;
sizeX = 1000;
sizeY = 1100;

params = struct();
params.D = 1e3;
params.kon = 0.1;
params.koff = 0.11;
params.maxSpeed = -1;
params.maxGapCloseFrames = 10;
params.minGapCloseTrackLength = 1;
params.minFinalTrackLength = 1;
params.maxPositionDisplacementSigma = 5;
params.maxFeatureDisplacementSigma = 5;
params.featureVar=[1, 10];

trackMethod = 'LAPTrack';
frameIdx=1000+sort(unique(randi(3*N,1,N)));
if numel(frameIdx)<N
    frameIdx = [frameIdx, frameIdx(end)+1:frameIdx(end)+N-numel(frameIdx)];
end
assert(numel(frameIdx)==N && numel(unique(frameIdx)) && issorted(frameIdx));
xs = rand(N,1)*sizeX+originX;
ys = rand(N,1)*sizeY+originY;
pos = [xs, ys];
posSE = 0.1*ones(N,2);

params.rho = numel(frameIdx)/(sizeX*sizeY*nFrames);
tk=Tracker(trackMethod, params);
tk.initializeTracks(frameIdx,pos,posSE);
stats = tk.getStats();

tracks = tk.getTracks();
tk.linkF2F();
[curIdx, nextIdx, costMat, connections, conn_costs] = tk.debugF2F(frameIdx(1));

tracks = tk.getTracks();
[costMat,connections, conn_costs]=tk.debugCloseGaps();

fM=full(costMat);
fM(fM==0)=Inf;
[rowsol, cost, v, u, rMat] = lapjv(fM);
% tracks = tk.getTracks();

% [curIdx, nextIdx, costMat, connections, conn_costs]=tk.debugF2F(0);
