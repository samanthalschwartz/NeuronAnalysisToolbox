testmask = zeros(size(aa.cleanedcargomask));
for ii = 1:size(aa.cleanedcargomask,3)
    testmask(:,:,ii) = aa.cleanedcargomask(:,:,ii).*ii;
end
testmasknan = testmask;
testmasknan(testmasknan==0) = 10^10;
tes2 = min(dip_image(testmasknan),[],3);
tes2(tes2==10^10)=0;

blackjet = flip(jet(255));
blackjet(1,:) = [0 0 0]; blackjet(end,:) = [1 1 1];


% get colorbar tick info
%            colorunit = size(labeledim,3)/255;
colorunit = (size(aa.cleanedcargomask,3)-6)/255;
numofcolbarval = 4;
colbarplace =[0:numofcolbarval]*255/4;
colbarval = [floor(colbarplace * colorunit)];

%-- now plot results
h = dipshow(tes2(:,:,5:end),blackjet);
dipmapping(h,[6 size(aa.cleanedcargomask,3)]);
diptruesize(h,100);


c = colorbar;
c.Location = 'WestOutside';
c.Ticks = colbarplace;
c.TickLabels = colbarval;
c.FontSize = 16;
c.Label.String = 'Time of First Appeareance After Release (min)';
c.Label.FontSize = 16;
h.OuterPosition = h.OuterPosition + [0 0 400 50];