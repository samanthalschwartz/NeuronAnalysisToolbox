function varargout = FRAP(varargin)
% FRAP MATLAB code for FRAP.fig
%      FRAP, by itself, creates a new FRAP or raises the existing
%      singleton*.
%
%      H = FRAP returns the handle to a new FRAP or the handle to
%      the existing singleton*.
%
%      FRAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FRAP.M with the given input arguments.
%
%      FRAP('Property','Value',...) creates a new FRAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FRAP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FRAP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FRAP

% Last Modified by GUIDE v2.5 25-Jun-2014 12:52:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FRAP_OpeningFcn, ...
                   'gui_OutputFcn',  @FRAP_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FRAP is made visible.
function FRAP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FRAP (see VARARGIN)

% Choose default command line output for FRAP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FRAP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FRAP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in choosefile.
function choosefile_Callback(hObject, eventdata, handles)
% hObject    handle to choosefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = uigetfile('.tif','Choose .tif file to import');
imginfo = imfinfo([pathname filename]);
w = imginfo(1).Width; h = imginfo(1).Height;
frames = length(imginfo);
img = zeros(h,w,frames);
tifobj = Tiff([pathname filename],'r');
for a = 1:frames
    tifobj.setDirectory(a)
    img(:,:,a) = tifobj.read();
    warning('off','last');
end
tifobj.close();
setappdata(handles.axes1,'imgdata',img);
axes(handles.axes1);
set(gcf,'toolbar','figure');
imagesc(min(img,[],3)); colormap('grey');
axis('equal'); axis('off')
filenamedisplay = fullfile(pathname,filename);
set(handles.filedisplay,'string',filenamedisplay);

% --- Executes on button press in frapselect.
function frapselect_Callback(hObject, eventdata, handles)
% hObject    handle to frapselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tstep = str2num(get(handles.fps,'string')); %#ok<ST2NM>
rectfrap = round(getrect(handles.axes1));
axes(handles.axes1);
recthandle1 = findobj(handles.axes1,'Type','Rectangle','EdgeColor',[0 0 1]);
if ~isempty(recthandle1)
    delete(recthandle1)
end
rectangle('Position',rectfrap,'EdgeColor',[0 0 1]);
img = getappdata(handles.axes1,'imgdata');
frapcutout = img(rectfrap(2):(rectfrap(2)+rectfrap(4)),rectfrap(1):(rectfrap(1)+rectfrap(3)),:);
frapdata1 = squeeze(mean(mean(frapcutout)));
frapdata2 = frapdata1 - min(frapdata1);
frapdata3 = frapdata2/max(frapdata2);
zeropoint = find(frapdata3 == 0);
t = [(-1*(length(1:(zeropoint-1)))):-1 0:(length(frapdata3)-zeropoint)]*tstep;
axes(handles.axes2);
cla(handles.axes2);
plot(t,frapdata3,'bo'); hold on;
xlabel('Time (s)')
ylabel('Intensity (au)')
xlim([min(t) max(t)])
ylim([0 1.1]);
bgdata = ones(length(frapdata3),1);
frapdatacorr = frapdata3./bgdata;
setappdata(handles.frapselect,'frapdata',[t' frapdata3 frapdatacorr])

% --- Executes on button press in bgselect.
function bgselect_Callback(hObject, eventdata, handles)
% hObject    handle to bgselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rectbg = round(getrect(handles.axes1));
axes(handles.axes1);
recthandle2 = findobj(handles.axes1,'Type','Rectangle','EdgeColor',[0 1 0]);
if ~isempty(recthandle2)
    delete(recthandle2)
end
rectangle('Position',rectbg,'EdgeColor',[0 1 0]);
img = getappdata(handles.axes1,'imgdata');
bgcutout = img(rectbg(2):(rectbg(2)+rectbg(4)),rectbg(1):(rectbg(1)+rectbg(3)),:);
bgdata1 = squeeze(mean(mean(bgcutout)));
bgdata3 = bgdata1/max(bgdata1);
frapdata = getappdata(handles.frapselect,'frapdata');
frapdata(:,3) = frapdata(:,2)./bgdata3;
setappdata(handles.frapselect,'frapdata',frapdata);
axes(handles.axes2);
cla(handles.axes2);
t = frapdata(:,1);
plot(t,frapdata(:,2),'bo',t,frapdata(:,3),'ko');
xlim([min(t) max(t)])
ylim([0 1.1]);
legend('Original FRAP data','Corrected FRAP data');


% --- Executes on button press in multiexp.
function multiexp_Callback(hObject, eventdata, handles)
% hObject    handle to multiexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiexp


% --- Executes on button press in fitcurve.
function fitcurve_Callback(hObject, eventdata, handles)
% hObject    handle to fitcurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
multiexp = get(handles.multiexp,'value');
frapdata = getappdata(handles.frapselect,'frapdata');
t = frapdata(:,1);
zeropoint = find(t == 0);
t2 = t(zeropoint:end);
datafrap = frapdata(zeropoint:end,3);
if ~multiexp
    ft = fittype('A*(1-exp(-tau*t))','independent','t','coefficients',{'A','tau'});
    [fm,gof] = fit(t2,datafrap,ft,'Startpoint',[0.5,0.05],'Lower',[0 0],'Upper',[1 inf]);
    axes(handles.axes2);
    cla(handles.axes2);
    plot(t2,frapdata(zeropoint:end,2),'bo',t2,datafrap,'ko'); hold on;
    legend('Original Data','Corrected Data')
    plot(fm); xlabel('Time (s)'); ylabel('Intensity (au)');
    xlim([min(t) max(t)])
    ylim([0 1.1]);
    hold off;
    coeffs = coeffvalues(fm);
    confints = confint(fm);
    A1conf = diff(confints(:,1))/2;
    thalf1 = log(2)/coeffs(2);
    thalf1confl = log(2)/confints(2,2);
    thalf1confu = log(2)/confints(1,2);
    thalf1conf = abs(thalf1confu - thalf1confl)/2;
    set(handles.A1txt,'string',num2str(coeffs(1)));
    set(handles.A1conftxt,'string',num2str(A1conf));
    set(handles.t1txt,'string',num2str(thalf1));
    set(handles.t1conftxt,'string',num2str(thalf1conf));
    set(handles.rsqtxt,'string',num2str(gof.adjrsquare));
    set(handles.A2txt,'string',num2str(0));
    set(handles.A2conftxt,'string',num2str(0));
    set(handles.t2txt,'string',num2str(0));
    set(handles.t2conftxt,'string',num2str(0));
elseif multiexp
    ft = fittype('A1*(1-exp(-tau1*t)) + A2*(1-exp(-tau2*t))','independent','t','coefficients',{'A1','tau1','A2','tau2'});
    [fm,gof] = fit(t2,datafrap,ft,'Startpoint',[0.5 0.05 0.1 0.01],'Lower',[0 0 0 0],'Upper',[1 inf 1 inf]);
    axes(handles.axes2);
    cla(handles.axes2);
    plot(t2,frapdata(zeropoint:end,2),'bo',t2,datafrap,'ko'); hold on;
    legend('Original Data','Corrected Data')
    plot(fm); xlabel('Time (s)'); ylabel('Intensity (au)');
    xlim([min(t) max(t)])
    ylim([0 1.1]);
    hold off;
    coeffs = coeffvalues(fm);
    confints = confint(fm);
    A1conf = diff(confints(:,1))/2;
    thalf1 = log(2)/coeffs(2);
    thalf1confl = log(2)/confints(2,2);
    thalf1confu = log(2)/confints(1,2);
    thalf1conf = abs(thalf1confu - thalf1confl)/2;
    A2conf = diff(confints(:,3))/2;
    thalf2 = log(2)/coeffs(4);
    thalf2confl = log(2)/confints(2,4);
    thalf2confu = log(2)/confints(1,4);
    thalf2conf = abs(thalf2confu - thalf2confl)/2;
    set(handles.A1txt,'string',num2str(coeffs(1)));
    set(handles.A1conftxt,'string',num2str(A1conf));
    set(handles.t1txt,'string',num2str(thalf1));
    set(handles.t1conftxt,'string',num2str(thalf1conf));
    set(handles.A2txt,'string',num2str(coeffs(3)));
    set(handles.A2conftxt,'string',num2str(A2conf));
    set(handles.t2txt,'string',num2str(thalf2));
    set(handles.t2conftxt,'string',num2str(thalf2conf));
    set(handles.rsqtxt,'string',num2str(gof.adjrsquare));
end

% --- Executes on button press in clearbutton.
function clearbutton_Callback(hObject, eventdata, handles)
% hObject    handle to clearbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes2);
set(handles.A1txt,'string','0');
set(handles.A1conftxt,'string','0');
set(handles.t1txt,'string','0');
set(handles.t1conftxt,'string','0');
set(handles.A2txt,'string','0');
set(handles.A2conftxt,'string','0');
set(handles.t2txt,'string','0');
set(handles.t2conftxt,'string','0');
set(handles.rsqtxt,'string','0');
recthandle3 = findobj(handles.axes1,'Type','Rectangle');
if ~isempty(recthandle3)
    delete(recthandle3)
end
setappdata(handles.frapselect,'frapdata',0);


function fps_Callback(hObject, eventdata, handles)
% hObject    handle to fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fps as text
%        str2double(get(hObject,'String')) returns contents of fps as a double


% --- Executes during object creation, after setting all properties.
function fps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in exportdata.
function exportdata_Callback(hObject, eventdata, handles)
% hObject    handle to exportdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filesave,pathsave] = uiputfile('*.xls','Save Excel file as...');
xlssave = fullfile(pathsave,filesave);
[~,filename,ext] = fileparts(get(handles.filedisplay,'string'));
imagename = [filename ext];
datawrite1 = {'Image File:', imagename};
xlswrite(xlssave,datawrite1,1,'A1:B1');
frapdata = getappdata(handles.frapselect,'frapdata');
headers = [{'Time (s)'} {'FRAP (original)'} {'FRAP (corrected)'}];
xlswrite(xlssave,headers,1,'A3:C3');
xlswrite(xlssave,frapdata,1,'A4');
xlswrite(xlssave,[{'A1'} {'A1 (+/-)'}],1,'E3:F3');
xlswrite(xlssave,[{'thalf1'} {'thalf1 (+/-)'}],1,'E5:F5');
xlswrite(xlssave,[{'A2'} {'A2 (+/-)'}],1,'E7:F7');
xlswrite(xlssave,[{'thalf2'} {'thalf2 (+/-)'}],1,'E9:F9');
xlswrite(xlssave,{'R-squared'},1,'E11');
xlswrite(xlssave,str2double(get(handles.A1txt,'string')),1,'E4');
xlswrite(xlssave,str2double(get(handles.A1conftxt,'string')),1,'F4');
xlswrite(xlssave,str2double(get(handles.t1txt,'string')),1,'E6');
xlswrite(xlssave,str2double(get(handles.t1conftxt,'string')),1,'F6');
xlswrite(xlssave,str2double(get(handles.A2txt,'string')),1,'E8');
xlswrite(xlssave,str2double(get(handles.A2conftxt,'string')),1,'F8');
xlswrite(xlssave,str2double(get(handles.t2txt,'string')),1,'E10');
xlswrite(xlssave,str2double(get(handles.t2conftxt,'string')),1,'F10');
xlswrite(xlssave,str2double(get(handles.rsqtxt,'string')),1,'E12');












