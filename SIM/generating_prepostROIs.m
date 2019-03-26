close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','G:\Hannah Dropbox SIM data\SIM_Files');
for ff= 1:numel(filepath)
load(filepath{ff});
obj.selectPrePostROI;
obj.save;
end


%%
joinchannels('rgb',obj.abeta.mask,obj.ch1.mask,obj.ch2.mask)
for pp = 1:numel(obj.results.selectedROIs)
   vertices = obj.results.selectedROIs{pp};
   p = patch('Vertices',vertices,'EdgeColor',[1 0 0],'Faces',1:size(vertices,1),'FaceAlpha',0);
end
