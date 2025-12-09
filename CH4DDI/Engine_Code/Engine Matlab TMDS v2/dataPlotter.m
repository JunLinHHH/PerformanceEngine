
function varargout = dataPlotter(varargin)
% DATAPLOTTER M-file for dataPlotter.fig
%      DATAPLOTTER, by itself, creates a new DATAPLOTTER or raises the existing
%      singleton*.
%
%      H = DATAPLOTTER returns the handle to a new DATAPLOTTER or the handle to
%      the existing singleton*.
%
%      DATAPLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATAPLOTTER.M with the given input arguments.
%
%      DATAPLOTTER('Property','Value',...) creates a new DATAPLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dataPlotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dataPlotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dataPlotter

% Last Modified by GUIDE v2.5 09-May-2020 18:07:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dataPlotter_OpeningFcn, ...
                   'gui_OutputFcn',  @dataPlotter_OutputFcn, ...
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


% --- Executes just before dataPlotter is made visible.
function dataPlotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dataPlotter (see VARARGIN)

set(0,'Units','pixels'); 
scnsize = get(0,'ScreenSize')
addpath(genpath('C:\Users\chris\OneDrive - UNSW\Desktop\Engine_Code\Engine Matlab TMDS v2'));
addpath(genpath('C:\Users\chris\OneDrive - UNSW\Desktop\Engine_Code\Engine Matlab TMDS v2\Required_toolbox'));
% Choose default command line output for dataPlotter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

function dataPlotter_OutputFcn(hObject, eventdata, handles, varargin)
% --- Executes during object creation, after setting all properties.
function ReadFilename1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ReadFilename1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in LoadDataButton.
function LoadDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles =guidata(hObject);
handles = SelectFileFunction(handles);
guidata(hObject,handles);

% --- Executes on button press in PlotDataButton.
function PlotDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles =guidata(hObject);
handles = Average_traces_TDMS(handles);
handles = Process_thermodynamics(handles);

% Save Data
DataSave(handles);
%Prepare plot screen
PreparePlot(handles);
uiwait;

guidata(hObject,handles);

function DataSave(handles)
filename = fullfile(handles.folder,handles.SetName,[handles.SetName,'_processed.mat']);
Traces = handles.traces;
Stats = handles.stats;
Thermo = handles.thermo;
TDCshift = handles.TDCshift;
save(filename,'Traces','Stats','Thermo','TDCshift');

% --- Executes during object creation, after setting all properties.
function input_Tw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_Tw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_Ti_Callback(hObject, eventdata, handles)
% hObject    handle to input_Ti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_Ti as text
%        str2double(get(hObject,'String')) returns contents of input_Ti as a double

global Ti
Ti = str2double(get(hObject,'String')); 


% --- Executes during object creation, after setting all properties.
function input_Ti_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_Ti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function input_RPM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_RPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function display_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SetNr_Callback(hObject, eventdata, handles)
% hObject    handle to SetNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SetNr as text
%        str2double(get(hObject,'String')) returns contents of SetNr as a double
handles = guidata(hObject);
if isfield(handles,'folder')==0
    handles.folder = '';
end
handles.LoadDataButton.Enable = 'on';
handles.PlotDataButton.Enable = 'off';
handles.SetName=['_Set',handles.SetNr.String,'_',handles.input_RPM.String,'_',handles.input_Ti.String,'_',handles.input_Tw.String];
handles.filelist = dir(fullfile(handles.folder,handles.SetName,[handles.SetName,'_Take*.tdms']));
if isempty(handles.filelist)==0
    handles = UpdateFilelist(handles);
    handles = LoadDataFunction(handles);
    handles.LoadDataButton.BackgroundColor = [0 0.8 0];
    handles.PlotDataButton.Enable = 'on';
else
    handles.LoadDataButton.BackgroundColor = [1 0 0];
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function SetNr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SetNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function InjT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InjT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function InjT_Callback(hObject, eventdata, handles)
% hObject    handle to InjT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InjT as text
%        str2double(get(hObject,'String')) returns contents of InjT as a double



function EditFolder_Callback(hObject, eventdata, handles)
% hObject    handle to EditFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditFolder as text
%        str2double(get(hObject,'String')) returns contents of EditFolder as a double


% --- Executes during object creation, after setting all properties.
function EditFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SelectFolder.
function SelectFolder_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.folder=handles.EditFolder.String;
handles.folder = uigetdir(handles.EditFolder.String);
handles.EditFolder.String = handles.folder;
if (max(handles.folder==0)) || exist(handles.folder,'dir')==0
    handles.LoadDataButton.Enable = 'off';
    handles.PlotDataButton.Enable = 'off';
    handles.LoadDataButton.BackgroundColor = [1 0 0];
else
    handles.LoadDataButton.Enable = 'on';
    handles.PlotDataButton.Enable = 'off';
    handles.SetName=['_Set',handles.SetNr.String,'_',handles.input_RPM.String,'_',handles.input_Ti.String,'_',handles.input_Tw.String];
    handles.filelist = dir(fullfile(handles.folder,handles.SetName,[handles.SetName,'_Take*.tdms']));
    if isempty(handles.filelist)==0
        handles = UpdateFilelist(handles);
        handles = LoadDataFunction(handles);
        handles.LoadDataButton.BackgroundColor = [0 0.8 0];
        handles.PlotDataButton.Enable = 'on';
    else
        handles.LoadDataButton.BackgroundColor = [1 0 0];
    end
    guidata(hObject,handles);
end


function handles = UpdateFilelist(handles)
    handles.ReadFilename1.String = handles.filelist(1).name;
    if length(handles.filelist)>1
        handles.ReadFilename2.String = handles.filelist(2).name;
    else
        handles.ReadFilename2.String = '';
    end
    if length(handles.filelist)>2
        handles.ReadFilename3.String = handles.filelist(3).name;
    else
        handles.ReadFilename3.String = '';
    end
    if length(handles.filelist)>3
        handles.ReadFilename4.String = handles.filelist(4).name;
    else
        handles.ReadFilename4.String = '';
    end
    if length(handles.filelist)>4
        handles.ReadFilename5.String = handles.filelist(5).name;
    else
        handles.ReadFilename5.String = '';
    end
    if length(handles.filelist)>5
        handles.ReadFilename6.String = handles.filelist(6).name;
    else
        handles.ReadFilename6.String = '';
    end
    if length(handles.filelist)>6
        handles.ReadFilename7.String = handles.filelist(7).name;
    else
        handles.ReadFilename7.String = '';
    end
    if length(handles.filelist)>7
        handles.ReadFilename8.String = handles.filelist(8).name;
    else
        handles.ReadFilename8.String = '';
    end
    if length(handles.filelist)>8
        handles.ReadFilename9.String = handles.filelist(9).name;
    else
        handles.ReadFilename9.String = '';
    end
    if length(handles.filelist)>9
        handles.ReadFilename10.String = handles.filelist(10).name;
    else
        handles.ReadFilename10.String = '';
    end

function handles = LoadDataFunction(handles)
    if isfield(handles,'raw_data')
        handles=rmfield(handles,'raw_data');
    end
    for i=1:length(handles.filelist)
        handles.raw_data(i)=TDMS_getStruct(fullfile(handles.folder,handles.SetName,handles.filelist(i).name));
    end

function handles = SelectFileFunction(handles)
   [handles.filelist,handles.folder] = uigetfile(fullfile(handles.folder,handles.SetName,'*.tdms'),'Select the files you wish to process','MultiSelect','on');
   handles.EditFolder.String = handles.folder;
   if iscell(handles.filelist)==0 && max(handles.filelist==0)
        handles=rmfield(handles,'filelist');
        handles.filelist = [];
    elseif iscell(handles.filelist)
        fl_copy = handles.filelist;
        handles=rmfield(handles,'filelist');
        for i=1:length(fl_copy)
            handles.filelist(i).name = fl_copy{i};
        end
    else
        filename =handles.filelist;
        handles=rmfield(handles,'filelist');
        handles.filelist(1).name = filename;
    end
    if isempty(handles.filelist)==0
        handles.AllegdedSetName = handles.SetName;
        handles.SetName = '';
        handles = UpdateFilelist(handles);
        handles = LoadDataFunction(handles);
        handles.LoadDataButton.BackgroundColor = [0 0.8 0];
        handles.PlotDataButton.Enable = 'on';
    else
        handles.LoadDataButton.BackgroundColor = [1 0 0];
        handles.PlotDataButton.Enable = 'off';
    end

function TDCshift_Callback(hObject, eventdata, handles)
% hObject    handle to TDCshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function TDCshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TDCshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Psens_Callback(hObject, eventdata, handles)
% hObject    handle to Psens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Psens as text
%        str2double(get(hObject,'String')) returns contents of Psens as a double


% --- Executes during object creation, after setting all properties.
function Psens_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Psens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
