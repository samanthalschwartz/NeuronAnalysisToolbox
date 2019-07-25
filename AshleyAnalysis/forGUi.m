% gui layout

%--- load files gui

%--- mask files layout
function sliderchanging
% Create figure window and components

fig = uifigure('Position',[100 100 350 275]);

cg = uigauge(fig,'Position',[100 100 120 120]);

sld = uislider(fig,...
               'Position',[100 75 120 3],...
               'ValueChangingFcn',@(sld,event) sliderMoving(event,cg));

end

% Create ValueChangingFcn callback
function sliderMoving(event,cg)
cg.Value = event.Value;
end