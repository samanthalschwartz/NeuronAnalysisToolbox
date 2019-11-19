clear all
close all

XY = load('../data/threeD.csv');

c = Clustering();
c.Results = 'results';

P{1} = XY(1:2:end, :);
P{2} = XY(2:2:end, :);
H_nm = 10000;
V_nm = 10000;
D_nm = 10000;
base_name = 'threeD';
particle_types = {'A', 'B'};
c.cluster_stats('PairwiseDist',       P, base_name, particle_types, ...
                                      H_nm, V_nm, D_nm);
c.cluster_stats('PairwiseMutualDist', P, base_name, particle_types, ...
                                      H_nm, V_nm, D_nm);
c.cluster_stats('Hopkins',            P, base_name, particle_types, ...
                                      H_nm, V_nm, D_nm);
c.cluster_stats('Ripley',             P, base_name, particle_types, ...
                                      H_nm, V_nm, D_nm);
c.cluster_stats('BivariateRipley',    P, base_name, particle_types, ...
                                      H_nm, V_nm, D_nm);
c.cluster_stats('Dendrogram',         P, base_name, particle_types, ...
                                      H_nm, V_nm, D_nm);
