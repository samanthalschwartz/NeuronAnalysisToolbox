filedir = 'C:\Users\schwsama\Documents\Data\SR\coverslip4-1a-647sr-rawdata';
filess = dir(fullfile(filedir,'*.tif'));

totframes = numel(filess);
breakdown = 5;
inc = ceil(totframes/breakdown);
increments = [1:inc:totframes,totframes];
movie = [];
frame1 = loadtiff(fullfile(filedir,filess(1).name));

for bb = 2:numel(increments)
    movsize = increments(bb)-increments(bb-1);
movie = zeros(size(frame1,1),size(frame1,2),movsize);
w = waitbar(0);
moviecnt = 1;
for ff = increments(bb-1):increments(bb)
    frame = loadtiff(fullfile(filedir,filess(ff).name));
%     frame = GeneralAnalysis.loadtiff_3ch(fullfile(filedir,filess(ff).name));
    movie(:,:,moviecnt) = squeeze(uint16(frame(:,:,1)));
    waitbar(moviecnt/movsize,w);
    moviecnt = moviecnt+1;
end
close(w)
sequence = movie;
savename = ['C:\Users\schwsama\Documents\Data\SR\coverslip4-1a-647srTry2_' num2str(increments(bb)) '.mat'];
save(savename,'sequence');
end