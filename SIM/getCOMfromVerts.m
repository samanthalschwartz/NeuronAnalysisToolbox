function coords = getCOMfromVerts(vertices,mask)
out = poly2mask(vertices(:,1)',vertices(:,2)',size(mask,1),size(mask,2));
lbch1 = label(mask,1);
labeledmask = lbch1(:,:,vertices(1,3)).*out;
test = labeledmask;
test(test==0) = NaN;
id = mode(single(test(:)));
imtest = label(lbch1==id,1);
msr = measure(imtest,single(imtest),({'Gravity'}));
coords = msr.Gravity;
end
