topdir = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\SIM_files\072617';
ids = [3 1 2];
channelorderingstr = {'chABeta','PSD95i','Gephyrin'}; % channel abeta, channel 1, channel2
files = uipickfiles('FilterSpec',topdir,'Prompt',['Pick Files for ' channelorderingstr{1} ,...
    '-' channelorderingstr{2} '-' channelorderingstr{3}]);

allAB = [];
allch1 = [];
allch2 = [];
for ff = 1:numel(files)
   load(files{ff});
   obj.make_maskchAB;
   msrAB = measure(obj.abeta.labeled_mask,obj.abeta.image,{'size'});
   msrch1 = measure(obj.ch1.mask,obj.abeta.image,{'size'});
   msrch2 = measure(obj.ch2.mask,obj.abeta.image,{'size'});
   allAB = [allAB, msrAB.size.*obj.XYpxsize.*obj.XYpxsize.*obj.Zpxsize];
   allch1 = [allch1, msrch1.size.*obj.XYpxsize.*obj.XYpxsize.*obj.Zpxsize];
   allch2 = [allch2, msrch2.size.*obj.XYpxsize.*obj.XYpxsize.*obj.Zpxsize];
end

figure; histogram(allAB); title('Abeta Pre threshold')
allAB(allAB<(6.1860e-04)) = [];
figure; histogram(allAB,200); title('Abeta Post threshold'); 
figure; histogram(allch1,100); title('PSD95ib')
figure; histogram(allch2,100); title('Gephyrin')