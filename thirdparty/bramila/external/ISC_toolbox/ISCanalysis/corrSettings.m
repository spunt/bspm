function varargout = corrSettings(varargin)
% CORRSETTINGS M-file for corrSettings.fig
%      CORRSETTINGS, by itself, creates a new CORRSETTINGS or raises the existing
%      singleton*.
%
%      H = CORRSETTINGS returns the handle to a new CORRSETTINGS or the handle to
%      the existing singleton*.
%
%      CORRSETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORRSETTINGS.M with the given input arguments.
%
%      CORRSETTINGS('Property','Value',...) creates a new CORRSETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before corrSettings_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to corrSettings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help corrSettings

% Last Modified by GUIDE v2.5 13-Jul-2010 20:46:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @corrSettings_OpeningFcn, ...
                   'gui_OutputFcn',  @corrSettings_OutputFcn, ...
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


% --- Executes just before corrSettings is made visible.
function corrSettings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to corrSettings (see VARARGIN)

handles.Pub = varargin{1};

set(handles.checkboxMean,'Value',handles.Pub.calcStandard)
set(handles.checkboxMedian,'Value',handles.Pub.calcStats)
set(handles.checkboxCorMat,'Value',handles.Pub.calcCorMatrices)
set(handles.editNrPerm,'String',num2str(handles.Pub.nrPermutations))
set(handles.editBatches,'String',num2str(handles.Pub.nrPermutationSets))
set(handles.editTotPerm,'String',...
    num2str(handles.Pub.nrPermutationSets*handles.Pub.nrPermutations),...
    'Enable','inactive')


% Choose default command line output for corrSettings
handles.output = handles.Pub;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes corrSettings wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = corrSettings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if(isfield(handles, 'output'))
    varargout{1} = handles.output;
    delete(handles.figure1)
else
    varargout{1} = [];
end

% --- Executes on button press in checkboxMean.
function checkboxMean_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.Pub.calcStandard = get(hObject,'Value');
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of checkboxMean
%if get(hObject,'Value')
%    set(handles.checkboxMedian,'Value',0)
%else
%    set(handles.checkboxMedian,'Value',1)
%end

% --- Executes on button press in checkboxMedian.
function checkboxMedian_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxMedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.Pub.calcStats = get(hObject,'Value');
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of checkboxMedian
%if get(hObject,'Value')
%    set(handles.checkboxMean,'Value',0)
%else
%    set(handles.checkboxMean,'Value',1)
%end

function editNrPerm_Callback(hObject, eventdata, handles)
% hObject    handle to editNrPerm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNrPerm as text
%        str2double(get(hObject,'String')) returns contents of editNrPerm as a double

val = str2double(get(hObject,'String'));
if isnan(val) || val <= 0
   set(hObject,'String',num2str(handles.Pub.nrPermutations))
else
    handles.Pub.nrPermutations = ceil(val);
    set(hObject,'String',num2str(val))
    set(handles.editTotPerm,'String',...
    num2str(handles.Pub.nrPermutations*handles.Pub.nrPermutationSets))
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editNrPerm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNrPerm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editBatches_Callback(hObject, eventdata, handles)
% hObject    handle to editBatches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBatches as text
%        str2double(get(hObject,'String')) returns contents of editBatches as a double
val = str2double(get(hObject,'String'));
if isnan(val) || val <= 0
   set(hObject,'String',num2str(handles.Pub.nrPermutationSets))
else
    handles.Pub.nrPermutationSets = ceil(val);
    set(hObject,'String',num2str(val))
    set(handles.editTotPerm,'String',...
    num2str(handles.Pub.nrPermutations*handles.Pub.nrPermutationSets))
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editBatches_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBatches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editTotPerm_Callback(hObject, eventdata, handles)
% hObject    handle to editTotPerm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTotPerm as text
%        str2double(get(hObject,'String')) returns contents of editTotPerm as a double


% --- Executes during object creation, after setting all properties.
function editTotPerm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTotPerm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonOK.
function pushbuttonOK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = handles.Pub;
guidata(hObject, handles);
uiresume(handles.figure1)

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1)

% --- Executes on button press in checkboxCorMat.
function checkboxCorMat_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxCorMat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxCorMat
handles.Pub.calcCorMatrices = get(hObject,'Value');
guidata(hObject, handles);

