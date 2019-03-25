image = logical(0*obj.ch1.mask);
for aa = 1:numel(obj.results.selectedROIs)
   vertices = obj.results.selectedROIs{aa};
   out = poly2mask(vertices(:,1)',vertices(:,2)',size(obj.ch1.mask,1),size(obj.ch1.mask,2));
   image = image | dip_image(out); 
end
%%
aa = 5;
vertices = obj.results.selectedROIs{aa};
out = poly2mask(vertices(:,1)',vertices(:,2)',size(obj.ch1.mask,1),size(obj.ch1.mask,2));
lbch1 = label(obj.ch1.mask,1);
labeledmask = lbch1(:,:,vertices(1,3)).*out;
test = labeledmask;
test(test==0) = NaN;
id = mode(single(test(:)));
imtest = label(lbch1==id,1);
msr = measure(imtest,single(imtest),({'Gravity'}));
%%
testim1= 0.*obj.ch1.mask;
testim1(round(msr.Gravity(1)),round(msr.Gravity(2)),round(msr.Gravity(3))) = 1;
joinchannels('rgb',imtest,testim1)