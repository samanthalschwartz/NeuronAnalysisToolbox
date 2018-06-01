function varargout = PlotTrackLengthHists_analyzeplot(varargin)
% PLOTTRACKLENGTHHISTS_ANALYZEPLOT MATLAB code for PlotTrackLengthHists_analyzeplot.fig
%      PLOTTRACKLENGTHHISTS_ANALYZEPLOT, by itself, creates a new PLOTTRACKLENGTHHISTS_ANALYZEPLOT or raises the existing
%      singleton*.
%
%      H = PLOTTRACKLENGTHHISTS_ANALYZEPLOT returns the handle to a new PLOTTRACKLENGTHHISTS_ANALYZEPLOT or the handle to
%      the existing singleton*.
%
%      PLOTTRACKLENGTHHISTS_ANALYZEPLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTTRACKLENGTHHISTS_ANALYZEPLOT.M with the given input arguments.
%
%      PLOTTRACKLENGTHHISTS_ANALYZEPLOT('Property','Value',...) creates a new PLOTTRACKLENGTHHISTS_ANALYZEPLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PlotTrackLengthHists_analyzeplot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PlotTrackLengthHists_analyzeplot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PlotTrackLengthHists_analyzeplot

% Last Modified by GUIDE v2.5 24-Sep-2014 14:35:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PlotTrackLengthHists_analyzeplot_OpeningFcn, ...
                   'gui_OutputFcn',  @PlotTrackLengthHists_analyzeplot_OutputFcn, ...
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


% --- Executes just before PlotTrackLengthHists_analyzeplot is made visible.
function PlotTrackLengthHists_analyzeplot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PlotTrackLengthHists_analyzeplot (see VARARGIN)

% Choose default command line output for PlotTrackLengthHists_analyzeplot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PlotTrackLengthHists_analyzeplot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PlotTrackLengthHists_analyzeplot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Filelistbox1.
function Filelistbox1_Callback(hObject, eventdata, handles) %---Files
% hObject    handle to Filelistbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Filelistbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Filelistbox1


% --- Executes during object creation, after setting all properties.
function Filelistbox1_CreateFcn(hObject, eventdata, handles) %---Files
% hObject    handle to Filelistbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selectFilesbutton1.
function selectFilesbutton1_Callback(hObject, eventdata, handles) %---Select Files
% hObject    handle to selectFilesbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in FileGroupQueuelistbox2.
function FileGroupQueuelistbox2_Callback(hObject, eventdata, handles) %--- Groupings
% hObject    handle to FileGroupQueuelistbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FileGroupQueuelistbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileGroupQueuelistbox2


% --- Executes during object creation, after setting all properties.
function FileGroupQueuelistbox2_CreateFcn(hObject, eventdata, handles) %--- Groupings
% hObject    handle to FileGroupQueuelistbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function GroupNameedit1_Callback(hObject, eventdata, handles) %---Group Name
% hObject    handle to GroupNameedit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GroupNameedit1 as text
%        str2double(get(hObject,'String')) returns contents of GroupNameedit1 as a double


% --- Executes during object creation, after setting all properties.
function GroupNameedit1_CreateFcn(hObject, eventdata, handles) %---Group Name
% hObject    handle to GroupNameedit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GroupFilesbutton3.
function GroupFilesbutton3_Callback(hObject, eventdata, handles) %--- Group Files
% hObject    handle to GroupFilesbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Analyzebutton4.
function Analyzebutton4_Callback(hObject, eventdata, handles) %--- Analyze and Plot
% hObject    handle to Analyzebutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function PlotTitleedit2_Callback(hObject, eventdata, handles) %-----Plot Title and Save Name
% hObject    handle to PlotTitleedit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PlotTitleedit2 as text
%        str2double(get(hObject,'String')) returns contents of PlotTitleedit2 as a double


% --- Executes during object creation, after setting all properties.
function PlotTitleedit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotTitleedit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on button press in checkbox1. %------------ HistogramCheckbx
% function checkbox1_Callback(hObject, eventdata, handles)
% % hObject    handle to checkbox1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of checkbox1
% 
% 
% % --- Executes on button press in checkbox2. %---------% CumProbDistCheckbx
% function checkbox2_Callback(hObject, eventdata, handles)
% % hObject    handle to checkbox2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in AddFilesIndividuallybutton5.
function AddFilesIndividuallybutton5_Callback(hObject, eventdata, handles) %--- add files individually
% hObject    handle to AddFilesIndividuallybutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in clearFilesbutton6.
function clearFilesbutton6_Callback(hObject, eventdata, handles) %--- clear files from file list
% hObject    handle to clearFilesbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles) %--- clear groups/files from group listing
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ResetAllbutton8.
function ResetAllbutton8_Callback(hObject, eventdata, handles) %--- re-set All
% hObject    handle to ResetAllbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in GroupstoPlotlistbox3.
function GroupstoPlotlistbox3_Callback(hObject, eventdata, handles)
% hObject    handle to GroupstoPlotlistbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns GroupstoPlotlistbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GroupstoPlotlistbox3


% --- Executes during object creation, after setting all properties.
function GroupstoPlotlistbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GroupstoPlotlistbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadGroupsbutton9.
function LoadGroupsbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to LoadGroupsbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Plotbutton10.
function Plotbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to Plotbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in clearGroupsbutton11.
function clearGroupsbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to clearGroupsbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ShowFilesButton12.
function ShowFilesButton12_Callback(hObject, eventdata, handles)
% hObject    handle to ShowFilesButton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function saveGroupDiredit3_Callback(hObject, eventdata, handles)
% hObject    handle to saveGroupDiredit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of saveGroupDiredit3 as text
%        str2double(get(hObject,'String')) returns contents of saveGroupDiredit3 as a double


% --- Executes during object creation, after setting all properties.
function saveGroupDiredit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveGroupDiredit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selectSaveDirbutton13.
function selectSaveDirbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to selectSaveDirbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
