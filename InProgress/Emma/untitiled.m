files = uipickfiles('FilterSpec','C:\Users\sammy\Dropbox\Sam Kennedy Lab\Emma');
filecheck = ones(numel(files),1);
while sum(filecheck)~=0
    ff = 1;
    [FILEPATH,NAME,EXT] = fileparts(files{ff});
    strid = regexp(NAME,'\d_C\d');
    chnum = strid+3;
    basename = NAME(1:strid);
    ids = cellfun(@(x) strfind(x,basename);
    
  
end
    
    
    
    
for ff = 1:numel(files)
  [FILEPATH,NAME,EXT] = fileparts(files{ff});
  strid = regexp(NAME,'\d_C\d');
  chnum = strid+3;
  basename = NAME(1:strid);
  
    
    
    
    
end