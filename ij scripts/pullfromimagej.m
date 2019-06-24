function pullfromimagej
fig = uifigure('CloseRequestFcn',@(fig,event) my_closereq(fig));
images = cell(MIJ.getListImages);
dd = uidropdown(fig, 'Items',images,...
    'Value', images{1},...
    'Position',[84 204 150 50],...
    'ValueChangedFcn',@(dd,event) selection(dd ));
end
% Create ValueChangedFcn callback:
function selection(dd)
        imgname =  dd.Value;
        output = MIJ.getImage(imgname);
        assignin('base','output', output);
end

function my_closereq(fig)
%         assignin('base','output', images{outputstr});
delete(fig); 
end

