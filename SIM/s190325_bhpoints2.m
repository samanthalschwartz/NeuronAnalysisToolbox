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
prepostdis_list = zeros(1,numel(obj.results.selectedROIs));
numabeta_aroundsynapse = zeros(1,numel(obj.results.selectedROIs));
wb = waitbar(0);
cnt = 0;
prepostdis_cnt = 0;
for aa = 1:numel(obj.results.selectedROIs)
    prepostdis_cnt = prepostdis_cnt+1;
    vertices = obj.results.selectedROIs{aa};
    zscale = obj.Zpxsize/obj.XYpxsize;
    coords_pre = getCOMfromVerts(vertices,obj.ch1.mask);
    coords_post = getCOMfromVerts(vertices,obj.ch2.mask);
    
    PrPo_vec = coords_post - coords_pre;
    PrPo_vec(3) = PrPo_vec(3)*zscale;
    Lmag = sqrt(sum(PrPo_vec).^2);
    prepostdis_list(prepostdis_cnt) = Lmag;
    
    lowbound = max(min(coords_pre,coords_post)-[20;20;3],[0; 0; 0]);
    maxbound = min(max(coords_pre,coords_post)+[20;20;3],(size(obj.abeta.image)-1)');
    
    abeta_im = obj.abeta.COM_image;
    labab = label(logical(abeta_im),1);
    cropped_labab = labab(floor(lowbound(1)):ceil(maxbound(1)),floor(lowbound(2)):ceil(maxbound(2)),floor(lowbound(3)):ceil(maxbound(3)));
    ids = unique(single(cropped_labab(:))); ids = ids(ids~=0);
    numabeta_aroundsynapse(prepostdis_cnt) = numel(ids);
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
        L1 = coords_pre;
        L2 = coords_post;
        c1 = (L1+L2)/2;
        L_vec = L2 - c1;
        L_vec(3) = L_vec(3)*zscale;
        pvec = p1 - c1;
        pvec(3) = pvec(3)*zscale;
        [Zdist, rdist] = vectordotproduct(Lvec,pvec);
        plotvals(:,cnt) = [rdist;Zdist];
    end
%     joinchannels('rgb',abetaCOM,pre_imgsub,post_imgsub)
    waitbar(aa/numel(obj.results.selectedROIs),wb);
end
close(wb);
%%
figure; histogram(prepostdis_list,10); title('Distance between Pre and Post');
figure; histogram(numabeta_aroundsynapse,10); title('Number of Abeta within Region');
figure; scatter(plotvals(1,:),plotvals(2,:),'*')

%% 
close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','G:\Hannah Dropbox SIM data\SIM_Files');
% wb = waitbar(0,'Looping Through Files');
prepostdis_list_all = cell(1,numel(filepath));
numabeta_aroundsynapse_all = cell(1,numel(filepath));
plotvals_all = cell(1,numel(filepath));
for ff= 1:numel(filepath)
    clear obj;
    load(filepath{ff});
    obj.abetaDensityAlongPrePost();
    obj.save(filepath{ff});

    prepostdis_list_all{ff} =  obj.results.prepostdis_list;
    numabeta_aroundsynapse_all{ff} = obj.results.numabeta_aroundsynapse;
    plotvals_all{ff} = obj.results.plotvals;
end
%- do this with 0 = center of pre/post alignment
% - plot histgram of pre/post distance
%%

close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','G:\Hannah Dropbox SIM data\SIM_Files');
prepostdis_list_all = cell(1,numel(filepath));
numabeta_aroundsynapse_all = cell(1,numel(filepath));
plotvals_all = cell(1,numel(filepath));
wb = waitbar(0,'Looping Through Run Files');
for ff= 1:numel(filepath)
    clear obj;
    load(filepath{ff});
    prepostdis_list_all{ff} =  obj.results.prepostdis_list;
    numabeta_aroundsynapse_all{ff} = obj.results.numabeta_aroundsynapse;
    plotvals_all{ff} = obj.results.plotvals;
    waitbar(ff/numel(filepath),wb);
end

numabeta = [];
all_dists = [];
plotsvals = [];
figure;hold on;
cols = lines(numel(plotvals_all));
for ii = 1:numel(plotvals_all)
   all_dists = [all_dists; prepostdis_list_all{ii}'];
   numabeta = [numabeta; numabeta_aroundsynapse_all{ii}'];
   plotsvals = cat(1,plotsvals,plotvals_all{ii}');
   scatter(plotvals_all{ii}(1,:),plotvals_all{ii}(2,:),'MarkerFaceColor',cols(ii,:),'MarkerEdgeColor',cols(ii,:))
end
figure; scatter(plotsvals(:,1),plotsvals(:,2),'*'); pbaspect([1 1 1])
ylabel('Distance from Center of Pre-Post')
xlabel('Distance Perpendicular to Synapse')
xlim([0 20]); ylim([-20 20])
hold on; plot([0:40],repmat(0,1,41),'w')

plotsvals(plotsvals(:,1)==0,:) = [];
[gca,N] = scatter2heatmap(plotsvals(:,1).*obj.XYpxsize,plotsvals(:,2).*obj.XYpxsize,30);

xlim([0 1]); ylim([-1 1]); pbaspect([1 1 1])
hold on; plot([0:40],repmat(0,1,41),'w','LineWidth',2)
ylabel('Distance from Center of Pre-Post')
xlabel('Distance Perpendicular to Synapse')

figure; histogram(all_dists.*obj.XYpxsize,20); title('Distance from Pre to Post');
xlabel('Distance in microns')
figure; histogram(numabeta,20); title('Number of Abeta COM within region Analyzed');
%%
goodids = plotsvals(:,1)<600;
alongZ = plotsvals(goodids,2).*obj.XYpxsize;
inc = 0.0667;
edges = -1:inc:1;
[N,edges] = histcounts(alongZ,edges);
figure; plot(N)
bins = -1+(inc/2):inc:1-(inc/2);
binsx = ones(1,numel(bins))
figure; plot(bins,N);
hold on; plot([0,0],[0,max(N)])
figure;imagesc(bins,binsx,N)


post = sum(plotsvals(:,2)>0);
pre = sum(plotsvals(:,2)<0);