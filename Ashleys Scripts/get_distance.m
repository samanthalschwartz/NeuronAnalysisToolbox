function distance = get_distance(city1, city2)
[num txt rawData] = xlsread(Distances.xls)
name1 = city1
name2 = city2
%look at all row elements in the first column
idx_1 = find(strcmp(rawData(:,1),name1))
idx_2 = find(strcmp(rawData(1,:),name2))
x = Aidx_1
y = idx_21
if isempty(idx_1)
    distance = -1
elseif isempty(idx_2)
    distance = -1
else
    distance = idx_2
end