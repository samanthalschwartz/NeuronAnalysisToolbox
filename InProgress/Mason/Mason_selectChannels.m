function [cellfillid, golgiid] =  Mason_selectChannels(image,cellfillid)
cellfillid.value = 1;
golgiid.value = 3;
numchannels = size(image,4);
x = 1:numchannels;
ch_opts = {string(x)};
fig = uifigure('CloseRequestFcn',@(fig, event)my_closereq(fig));
pannel_cellfill = uipanel('Parent',fig,'Title','Select Cell Fill Channel','Position', [ 20    20   260   221])
pannel_golgi = uipanel('Parent',fig,'Title','Select Golgi Channel','Position',[ 260    20   260   221])
ch_opts = cellmap(@int2str, num2cell(x));
dd = uidropdown(fig,...
    'Position',[50 150 100 22],...
    'Items',ch_opts,...
    'Value', ch_opts{1},...
    'ValueChangedFcn',@(dd,event) cellfillselection(dd,cellfillid));
cc = uidropdown(fig,...
    'Position',[280 150 100 22],...
    'Items',ch_opts,...
    'Value', ch_opts{3},...
    'ValueChangedFcn',@(cc,event) golgiselection(cc,golgiid));
end
selection(dd,p)
function cellfillid = cellfillselection(dd,cellfillid)
 val = dd.Value;
        p.Color = val;       
cellfillid.value = str2double(dd.Value);
end

function golgiid = golgiselection(cc,golgiid)
golgiid.value = str2double(cc.Value);
end

function my_closereq()
print('testing');

end