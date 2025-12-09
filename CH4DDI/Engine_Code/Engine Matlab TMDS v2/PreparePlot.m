function varargout = PreparePlot(varargin)
% PREPAREPLOT MATLAB code for PreparePlot.fig
%      PREPAREPLOT, by itself, creates a new PREPAREPLOT or raises the existing
%      singleton*.
%
%      H = PREPAREPLOT returns the handle to a new PREPAREPLOT or the handle to
%      the existing singleton*.
%
%      PREPAREPLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPAREPLOT.M with the given input arguments.
%
%      PREPAREPLOT('Property','Value',...) creates a new PREPAREPLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PreparePlot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PreparePlot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PreparePlot

% Last Modified by GUIDE v2.5 09-May-2020 21:41:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PreparePlot_OpeningFcn, ...
                   'gui_OutputFcn',  @PreparePlot_OutputFcn, ...
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


% --- Executes just before PreparePlot is made visible.
function PreparePlot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PreparePlot (see VARARGIN)

% Choose default command line output for PreparePlot
handles.output = hObject;
Hin = varargin{1};

%Update indicators
handles.PeakP.String = num2str(Hin.stats.P_max,'%3.1f');
handles.CoV.String = num2str(Hin.stats.CoV,'%0.4f');
handles.CA10.String = num2str(Hin.thermo.CA10,'%2.1f');
handles.CA50.String = num2str(Hin.thermo.CA50,'%2.1f');
handles.CA90.String = num2str(Hin.thermo.CA90,'%2.1f');
handles.ID.String = num2str(Hin.thermo.IgnDelay,'%2.1f');
handles.BurnDur.String = num2str(Hin.thermo.BurnDur,'%2.1f');
handles.IMEP.String = num2str(Hin.stats.IMEP,'%2.3f');
handles.PeakPRR.String = num2str(Hin.stats.PRR_max,'%2.1f');
handles.TotalHRR.String = num2str(Hin.thermo.totHRR,'%4.0f');

%Plot PV diagram
axes(handles.axes3);
plot(Hin.traces.Vol*1000,Hin.traces.Pmean,'LineWidth',1);
hold on; grid on; box on;
title('PV diagram');
xlabel('Volume [dm3]'); ylabel('Pressure [bar]')

axes(handles.axes7);
plot(Hin.traces.Vol*1000,Hin.traces.Pmean,'LineWidth',1);
hold on;
xlim([0.025 0.05]);ylim([35,100]);

%Plot gas exchange
axes(handles.axes4);
plot(Hin.traces.Vol,Hin.traces.Pmean,'LineWidth',1);
hold on; grid on; box on;
title('Gas Exchange');ylim([0.6,1.8]);
xlabel('Volume [dm3]'); ylabel('Pressure [bar]')

%Plot pressure/injection
axes(handles.axes1);
plot(Hin.traces.CA,Hin.traces.Pmean,'LineWidth',1);
hold on; grid on; box on;
title('Pressure trace');xlim([-20,140]);
xlabel('CA [°]'); ylabel('Pressure [bar]')

axes(handles.axes2);
plot(Hin.traces.CA,Hin.traces.InjMean,'LineWidth',1,'COlor','r');
hold on;
xlim([-20,140]);
ylabel('Injection signal');
set(gca,'XTick',[],'Color','none','YAxisLocation','right');

axes(handles.axes8);
plot(Hin.traces.CA,Hin.traces.Pmean,'LineWidth',1);
hold on; grid on; box on;
xlim([-10,30]);


%Plot HRR/cumHRR
axes(handles.axes5);
plot(Hin.traces.CA(1:end-1),Hin.thermo.aHRR,'LineWidth',1);
hold on; grid on; box on;
title('HRR analysis');xlim([-20,140]);
xlabel('CA [°]'); ylabel('aHRR [J/°]');

axes(handles.axes6);
plot(Hin.traces.CA(1:end-1),Hin.thermo.cumHRR,'LineWidth',1,'COlor','r');
hold on; box on;
xlim([-20,140]);
ylabel('Cumulative HRR [J]');
set(gca,'XTick',[],'Color','none','YAxisLocation','right');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PreparePlot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PreparePlot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function PeakP_Callback(hObject, eventdata, handles)
% hObject    handle to PeakP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PeakP as text
%        str2double(get(hObject,'String')) returns contents of PeakP as a double


% --- Executes during object creation, after setting all properties.
function PeakP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PeakP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CoV_Callback(hObject, eventdata, handles)
% hObject    handle to CoV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CoV as text
%        str2double(get(hObject,'String')) returns contents of CoV as a double


% --- Executes during object creation, after setting all properties.
function CoV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CoV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CA10_Callback(hObject, eventdata, handles)
% hObject    handle to CA10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CA10 as text
%        str2double(get(hObject,'String')) returns contents of CA10 as a double


% --- Executes during object creation, after setting all properties.
function CA10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CA10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CA50_Callback(hObject, eventdata, handles)
% hObject    handle to CA50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CA50 as text
%        str2double(get(hObject,'String')) returns contents of CA50 as a double


% --- Executes during object creation, after setting all properties.
function CA50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CA50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CA90_Callback(hObject, eventdata, handles)
% hObject    handle to CA90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CA90 as text
%        str2double(get(hObject,'String')) returns contents of CA90 as a double


% --- Executes during object creation, after setting all properties.
function CA90_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CA90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ID_Callback(hObject, eventdata, handles)
% hObject    handle to ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ID as text
%        str2double(get(hObject,'String')) returns contents of ID as a double


% --- Executes during object creation, after setting all properties.
function ID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BurnDur_Callback(hObject, eventdata, handles)
% hObject    handle to BurnDur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BurnDur as text
%        str2double(get(hObject,'String')) returns contents of BurnDur as a double


% --- Executes during object creation, after setting all properties.
function BurnDur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BurnDur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IMEP_Callback(hObject, eventdata, handles)
% hObject    handle to IMEP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IMEP as text
%        str2double(get(hObject,'String')) returns contents of IMEP as a double


% --- Executes during object creation, after setting all properties.
function IMEP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IMEP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PeakPRR_Callback(hObject, eventdata, handles)
% hObject    handle to PeakPRR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PeakPRR as text
%        str2double(get(hObject,'String')) returns contents of PeakPRR as a double


% --- Executes during object creation, after setting all properties.
function PeakPRR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PeakPRR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TotalHRR_Callback(hObject, eventdata, handles)
% hObject    handle to TotalHRR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TotalHRR as text
%        str2double(get(hObject,'String')) returns contents of TotalHRR as a double


% --- Executes during object creation, after setting all properties.
function TotalHRR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TotalHRR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
