function coords = getCOMfromVerts(vertices,mask)
out = poly2mask(vertices(:,1)',vertices(:,2)',size(mask,1),size(mask,2));
lbch1 = label(mask,1);
if size(vertices,1)<3 
    disp('Not enough vertices to make a region');
    return;
end
if size(vertices,2) == 2
    labeledmask = lbch1.*out;
else
    labeledmask = lbch1(:,:,vertices(1,3)).*out;
end 
test = labeledmask;
if sum(test(:)) == 0
    coords=[];
    return;
end
test(test==0) = NaN;
id = mode(single(test(:)));
imtest = label(lbch1==id,1);
msr = measure(imtest,single(imtest),({'Gravity'}));
coords = msr.Gravity;
end
