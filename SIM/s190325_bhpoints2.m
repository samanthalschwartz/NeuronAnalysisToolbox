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
%% Zdistance is incorrect
totab = sum(obj.abeta.COM_image);
plotvals = zeros(2,totab);

wb = waitbar(0);
cnt = 0;
for aa = 1:numel(obj.results.selectedROIs)
    vertices = obj.results.selectedROIs{aa};
    coords_pre = getCOMfromVerts(vertices,obj.ch1.mask);
    coords_post = getCOMfromVerts(vertices,obj.ch2.mask);
    
    lowbound = max(min(coords_pre,coords_post)-[20;20;3],[0; 0; 0]);
    maxbound = min(max(coords_pre,coords_post)+[20;20;3],(size(obj.abeta.image)-1)');
    
    abeta_im = obj.abeta.COM_image;
    labab = label(logical(abeta_im),1);
    cropped_labab = labab(floor(lowbound(1)):ceil(maxbound(1)),floor(lowbound(2)):ceil(maxbound(2)),floor(lowbound(3)):ceil(maxbound(3)));
    ids = unique(single(cropped_labab(:))); ids = ids(ids~=0);
    for ii =1 :numel(ids)
        cnt= cnt+1;
    c = findcoord(labab==ids(ii));
    p1 = c';
    
%     pre_img = 0*obj.ch1.mask;
%     pre_img( round(coords_pre(1)),round(coords_pre(2)),round(coords_pre(3)) ) = 1;
%     pre_imgsub = pre_img(floor(lowbound(1)):ceil(maxbound(1)),floor(lowbound(2)):ceil(maxbound(2)),floor(lowbound(3)):ceil(maxbound(3)))
% 
%     post_img = 0*obj.ch1.mask;
%     post_img( round(coords_post(1)),round(coords_post(2)),round(coords_post(3)) ) = 1;
%     post_imgsub = post_img(floor(lowbound(1)):ceil(maxbound(1)),floor(lowbound(2)):ceil(maxbound(2)),floor(lowbound(3)):ceil(maxbound(3)))

   
    pre = coords_pre;
    post = coords_post;
    [Zdist, rdist] = vectordotproduct(pre,post,p1);
    plotvals(:,cnt) = [rdist;Zdist];
    end
%     joinchannels('rgb',abetaCOM,pre_imgsub,post_imgsub)
    waitbar(aa/numel(obj.results.selectedROIs),wb);
end
close(wb);
%%
%- do this with 0 = center of pre/post alignment
% - plot histgram of pre/post distance
