
points=rand(2,100)*10;
tree=KDTree(points);
q=tree.query([2 3], [7 8]);
min(q,[],2)
max(q,[],2)
