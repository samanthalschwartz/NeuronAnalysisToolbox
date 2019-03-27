function [Zdist, rdist] = vectordotproduct(lvec,pvec)
% L1 pre-synaptic marker
% L2 post-synaptic marker
% p1 - abeta point
Zdist = dot(lvec,pvec)/sqrt(sum(lvec.^2));
rdist = sqrt(abs(sum(pvec.^2) - Zdist^2));
end