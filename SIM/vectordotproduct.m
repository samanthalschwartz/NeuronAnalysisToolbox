function [Zdist, rdist] = vectordotproduct(L1,L2,p1)
% L1 pre-synaptic marker
% L2 post-synaptic marker
% p1 - abeta point
lvec = L2-L1;
pvec = p1-L1;
Zdist = dot(lvec,pvec)/sqrt(sum(lvec.^2));
rdist = sqrt(abs(sum(pvec.^2) - Zdist^2));
end