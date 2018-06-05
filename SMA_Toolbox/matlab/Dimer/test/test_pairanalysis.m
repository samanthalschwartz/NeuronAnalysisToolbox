exposureT = 0.002;
D_A = 0.100;
D_B = D_A;
D_AB = 0.125;
meanV = 0.020^2;
rho = 0.020;
N=25;
p=zeros(N,18);
p(:,1) = (0:N-1)'*exposureT;
[p(:,2),~,p(:,10)] = DEstimator.simulate1D(D_A, N, meanV, exposureT, exposureT, 'gamma');
[p(:,3),~,p(:,11)] = DEstimator.simulate1D(D_A, N, meanV, exposureT, exposureT, 'gamma');
[p(:,6),~,p(:,14)] = DEstimator.simulate1D(D_B, N, meanV, exposureT, exposureT, 'gamma');
[p(:,7),~,p(:,15)] = DEstimator.simulate1D(D_B, N, meanV, exposureT, exposureT, 'gamma');
% p(:,[6,7]) = p(:,[2,3]);
p(:,[10,11,14,15]) = sqrt(p(:,[10,11,14,15]));
dit=InteractionChangePoint(D_A, D_B, D_AB, exposureT, rho);
tic;
dit.initializePair(p);
fprintf('Comp time [N:%i] [Entries:%i] [Size:%iMb] - %gs\n',N, N^2-N, round(N^2*8/2^20), toc); 
dit.plotPairMatrixDistances();
dit.plotLengthNSegmentComparison(1);


hd=SimulateHomoDimer2D(25,5000);
% hd.koff = 0;
hd.simulateTracks();
hd.plotParticles();
iit=InteractionChangePoint(D_A, D_B, D_AB, exposureT, rho);
iit.initializePair(hd.observePair)

% dit.initializePair(p);
% dit.plotPairMatrixDistances()
% dit.plotSingleJumpComparison()



% freeLLH = it.computeFreePairLikelihood(pair_mat);
% freeLLH_dst = dit.computeFreePairLikelihoodDEstimator(p);
% % freeLLH_delta = freeLLH - freeLLH_dst;
% boundLLH = it.computeBoundPairLikelihood(pair_mat);
% 
% f=figure();
% surface(1:size(pair_mat,1),1:size(pair_mat,1),boundLLH-freeLLH)
% pbaspect([1,1,1])
% 

