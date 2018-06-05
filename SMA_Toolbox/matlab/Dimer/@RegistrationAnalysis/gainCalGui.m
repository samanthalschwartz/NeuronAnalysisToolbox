function gainCalGui(obj, update_CB)
if nargin==1
    update_CB=[];
end
if ishandle(obj.gainCalGuiFig)
    figure(obj.gainCalGuiFig);
    return
end
gainCalFigs=[];

CCDBackgroundFrames=[]; %Loaded Background frames
CCDGainFrames=[]; %Loaded Gain Calibration Frames

font_size=10;
sp=3;
uH=25;
boarder=10;
fig_sz=[750 550];
but_sz=[140 uH];
label_sz=[250 uH*0.75];
edit_sz=[300 uH];
path_edit_sz=[fig_sz(1)-but_sz(1)-2*boarder-sp uH];
guiFig = figure('Units','pixels','Position',[300 300 fig_sz],'Resize','off',...
    'WindowStyle','normal','MenuBar','none','ToolBar','none','NumberTitle','off','Name','CCD Gain Calibration');
set(guiFig,'CloseRequestFcn',@close_CB);
obj.gainCalGuiFig=guiFig;

bot_row.pos=[boarder boarder fig_sz(1)-2*boarder uH]; 
bot_row.names={'Calibrate','ViewGainCal', 'ViewBackground','Clear'};
bot_row.CBs={@calibrate_CB,@viewGainCal_CB,@viewBackground_CB,@clear_CB};
GUIBuilder.buttonRow(guiFig, bot_row.pos, but_sz, bot_row.names, bot_row.CBs);
uicontrol('Parent',guiFig,'Style','pushbutton','String','Quit',...
         'Position',[fig_sz(1)-but_sz(1)-boarder boarder but_sz],'Callback',@close_CB);

handles.CCDBackgroundPath = uicontrol('Parent',guiFig,'Style','edit','String',obj.getFilePath('CCDBackground'),...
                            'Position',[boarder+but_sz(1)+sp boarder+uH+3*sp path_edit_sz],...
                            'Callback',@loadCCDBackground_CB);
handles.CCDBackgroundBut = uicontrol('Parent',guiFig,'Style','pushbutton','String','Load BG Image',...
                            'Position',[boarder boarder+uH+3*sp but_sz],...
                            'Callback',@loadCCDBackground_CB);
uicontrol('Parent',guiFig,'Style','text','String','CCD Background Calibration Sequence',...
                            'Position',[boarder boarder+2*(uH+sp)+2*sp fig_sz(1)-2*boarder 16],...
                            'HorizontalAlignment','left','FontSize',font_size);
handles.CCDGainCalPath = uicontrol('Parent',guiFig,'Style','edit','String',obj.getFilePath('CCDGainCal'),...
                            'Position',[boarder+but_sz(1)+sp boarder+3*(uH+sp)+2*sp path_edit_sz],...
                            'Callback',@loadCCDGainCal_CB);
handles.CCDGainCalBut = uicontrol('Parent',guiFig,'Style','pushbutton','String','Load Beads Image',...
                            'Position',[boarder boarder+3*(uH+sp)+2*sp but_sz],...
                            'Callback',@loadCCDGainCal_CB);
uicontrol('Parent',guiFig,'Style','text','String','CCD Out of focus Beads Calibration Sequence',...
                            'Position',[boarder boarder+4*(uH+sp)+2*sp fig_sz(1)-2*boarder 16],...
                            'HorizontalAlignment','left','FontSize',font_size);
annotation(guiFig, 'rectangle','Units','pixels','Position',[boarder-2*sp boarder+uH+2*sp fig_sz(1)-2*boarder+4*sp 4*(uH+sp)+sp])


input_loc=[boarder boarder+5*(uH+sp)+4*sp];
handles.CCDBackgroundEdit = labeled_editbox('CCD Background',obj.CCDBackground,input_loc);
handles.CCDBackgroundEdit.Callback = @setCCDBackground_CB;
handles.CCDGainEdit = labeled_editbox('CCD Gain [e-/ADU]',obj.CCDGain,input_loc+[0 (uH+sp)]);
handles.CCDGainEdit.Callback = @setCCDGain_CB;

axes_loc=[boarder+44,boarder+8*(uH+sp)+8*sp];
axs_sz=[250, 250];
handles.gainCalAxes= axes('Units', 'Pixels','Visible','off',...
                     'Position',[axes_loc axs_sz]);
updateAxes(handles.gainCalAxes, 'CCD Gain Cal Frames', CCDGainFrames);
handles.bgAxes= axes('Units', 'Pixels','Visible','off',...
                     'Position',[axes_loc+[360+2*sp 0] axs_sz]);
updateAxes(handles.bgAxes, 'Background Frames', CCDBackgroundFrames);

function editH=labeled_editbox(name,val,loc)
    %Make a labeled edit box
    uicontrol('Parent',guiFig,'Style','text','String',name,...
              'Position',[loc label_sz],...
              'HorizontalAlignment','left','FontSize',font_size);
    editH=uicontrol('Parent',guiFig,'Style','edit','String',num2str(val),...
                    'Position',[loc(1)+label_sz(1)+sp loc(2) edit_sz],...
                    'HorizontalAlignment','right','FontSize',font_size);
end

function frames=updateAxes(handle, plot_title, frames)
    axes(handle);
    imagesc(single(frames(:,:,1)));
    colormap(gray);
    title(plot_title);
    xlabel('X');
    ylabel('Y');
    set(handle,'Visible','on');
end

function calibrate_CB(~,~)
    if ~isempty(gainCalFigs)
        gainCalFigs=gainCalFigs(arrayfun(@ishandle,gainCalFigs));
        closeFigs=intersect(findobj('Type','figure'),gainCalFigs);
        if ~isempty(closeFigs)
            close(closeFigs);
        end
        gainCalFigs=[];
    end
    try
        gainPath = handles.CCDGainCalPath.String;
        bgPath = handles.CCDBackgroundPath.String;
        figHs = obj.calibrateCCD(gainPath, bgPath);
        gainCalFigs=[gainCalFigs figHs(:)'];
        obj.appendOpenFigs(figHs);
    catch err
        disp(getReport(err));
        errordlg(err.message,err.identifier);
        return
    end
    handles.CCDBackgroundEdit.String = num2str(obj.CCDBackground);
    handles.CCDGainEdit.String = num2str(obj.CCDGain);
    if ~isempty(update_CB)
        update_CB();
        figure(guiFig);
    end
end

function loadCCDBackground_CB(~,~)
    currPath = handles.CCDBackgroundPath.String;
    if isempty(currPath)
        currPath = handles.CCDGainCalPath.String;
        if isempty(currPath)
            currPath = obj.workingDir;
        else
            [currPath,~,~] = fileparts(currPath);
        end
    end
    [filename, pathname] = uigetfile(obj.LoadableRawDataFormats, 'Load CCD Background File', currPath);
    if filename
        newPath = fullfile(pathname,filename);
        S = load(newPath,'-mat');
        if isfield(S,'sequence')
            CCDBackgroundFrames = S.sequence;
        elseif isfield(S,'Data')
            CCDBackgroundFrames = S.Data;
        else
            error('gainCalGui:BadSequence','Unknown image format.');
        end
        handles.CCDBackgroundPath.String = newPath;
        updateAxes(handles.bgAxes, 'CCD Background', CCDBackgroundFrames);
    end
end



function loadCCDGainCal_CB(~,~)
    currPath=get(handles.CCDGainCalPath,'String');
    if isempty(currPath)
        currPath=get(handles.CCDBackgroundPath,'String');
        if isempty(currPath)
            currPath=obj.workingDir;
        else
            [currPath,~,~]=fileparts(currPath);
        end
    end
    [filename, pathname]=uigetfile(obj.LoadableRawDataFormats, 'Load CCD Gain File', currPath);
    if filename
        newPath=fullfile(pathname,filename);
        S = load(newPath,'-mat');
        if isfield(S,'sequence')
            CCDGainFrames = S.sequence;
        elseif isfield(S,'Data')
            CCDGainFrames = S.Data;
        else
            error('gainCalGui:BadSequence','Unknown image format.');
        end
        set(handles.CCDGainCalPath,'String',newPath);
        updateAxes(handles.gainCalAxes, 'CCD Gain Calibration (beads)', CCDGainFrames);
    end
end

function setCCDBackground_CB(H,~)
    val=str2double(H.String);
    if val>=0 && isfinite(val) && val~=obj.CCDBackground
        obj.CCDBackground=val;
        if obj.calibrated
            obj.calibrateCCD();
        end
        if ~isempty(update_CB)
            update_CB;
            figure(guiFig);
        end
    else
        H.String = num2str(obj.CCDBackground);
    end
end

function setCCDGain_CB(H,~)
    val = str2double(H.String);
    if val>=0 && isfinite(val) && val~=obj.CCDGain
        obj.CCDGain=val;
        if obj.calibrated
            obj.calibrateCCD();
        end
        if ~isempty(update_CB)
            update_CB;
            figure(guiFig);
        end
    else
        H.String = num2str(obj.CCDGain);
    end
end

function clear_CB(~,~)
    CCDBackgroundFrames = [];
    CCDGainFrames = [];
    set(handles.gainCalAxes,'Visible','off');
    cla(handles.gainCalAxes);
    set(handles.bgAxes,'Visible','off');
    cla(handles.bgAxes);
    set(handles.CCDGainCalPath,'String',[]);
    set(handles.CCDBackgroundPath,'String',[]);
end

function close_CB(~,~)
    closeFigs=obj.figHs( arrayfun(@ishandle,obj.figHs) );
    if ~isempty(closeFigs)
        close(closeFigs);
    end
    delete(guiFig);
    obj.gainCalGuiFig=[];
end

function viewBackground_CB(~,~)
    frames = CCDBackgroundFrames;
    if isempty(frames); return; end
    h=obj.viewMaximizedDipFig(frames);
    gainCalFigs(end+1)=h;
    obj.appendOpenFigs(h);
end

function viewGainCal_CB(~,~)
    frames = CCDGainFrames;
    if isempty(frames); return; end
    h=obj.viewMaximizedDipFig(frames);
    gainCalFigs(end+1)=h;
    obj.appendOpenFigs(h);
end

end %gainCalgui
