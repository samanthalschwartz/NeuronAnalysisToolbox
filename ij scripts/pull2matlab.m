% gui to calculate shift and then apply shift to any number of images
% select timeseries to drift correct shift
% 
function pull2matlab()
% Create figure and components:

fig = uifigure('CloseRequestFcn',@(fig, event)my_closereq(fig));

filenames = cell(MIJ.getListImages);
% Create drop-down component
dd = uidropdown(fig,...
    'Position',[430 210 100 22],...
    'Items',filenames,...
    'Value',filenames{1},...   
    'ValueChangedFcn',@(dd,event) selection(dd));
end
function my_closereq(fig)
delete(fig);
end
% Create ValueChangedFcn callback:
    function selection(dd)
       val = dd.Value;
       currim = MIJ.getImage(val);
       
       assignin('base','imIJ', MIJ.getImage(val));
    end