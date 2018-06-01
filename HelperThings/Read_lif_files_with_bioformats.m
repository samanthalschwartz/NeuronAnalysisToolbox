

%Load with bioformats
function [CO,A,B] = Read_lif_files_with_bioformats(filepath_name)
BFim=bfopen(filepath_name)
seriesCount = size(BFim, 1);
SZ=size(BFim{1,1}{1,1})
A=newim([SZ seriesCount]); %could make this
B=newim([SZ  seriesCount]);
for ii=0:seriesCount-1
    A(:,:,ii) = single(BFim{ii+1,1}{1,1});
    B(:,:,ii) = single(BFim{ii+1,1}{2,1});
end
CO=joinchannels('RGB',A,B);
end

