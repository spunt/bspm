function varargout = WaveIDT1(varargin)
%WAVEIDT1 M-file for WaveIDT1.fig
%      WAVEIDT1, by itself, creates a new WAVEIDT1 or raises the existing
%      singleton*.
%
%      H = WAVEIDT1 returns the handle to a new WAVEIDT1 or the handle to
%      the existing singleton*.
%
%      WAVEIDT1('Property','Value',...) creates a new WAVEIDT1 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to WaveIDT1_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      WAVEIDT1('CALLBACK') and WAVEIDT1('CALLBACK',hObject,...) call the
%      local function named CALLBACK in WAVEIDT1.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WaveIDT1

% Last Modified by GUIDE v2.5 11-May-2011 10:40:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
% waveidtb_defaults;
global GUI_GLOBAL_DATA;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WaveIDT1_OpeningFcn, ...
                   'gui_OutputFcn',  @WaveIDT1_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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



% --- Executes just before WaveIDT1 is made visible.
function WaveIDT1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for WaveIDT1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes WaveIDT1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = WaveIDT1_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setupDenoising.
function setupDenoising_Callback(hObject, eventdata, handles)
% hObject    handle to setupDenoising (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% waveidtb_defaults;
hObject = waveidtb_getOutputDir(hObject);
global GUI_GLOBAL_DATA;
waveidtb_getData(hObject);
guidata( hObject , GUI_GLOBAL_DATA );

% --- Executes during object creation, after setting all properties.
function setupDenoising_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setupDenoising (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over setupDenoising.
function setupDenoising_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to setupDenoising (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in loadParameters.
function loadParameters_Callback(hObject, eventdata, handles)
% hObject    handle to loadParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% waveidtb_defaults;
hObject = waveidtb_loadParameters(hObject);
global GUI_GLOBAL_DATA;
guidata(hObject, GUI_GLOBAL_DATA);

% --- Executes during object creation, after setting all properties.
function loadParameters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over loadParameters.
function loadParameters_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to loadParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in runDenoising.
function runDenoising_Callback(hObject, eventdata, handles)
% hObject    handle to runDenoising (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global GUI_GLOBAL_DATA;
% guidata(hObject, GUI_GLOBAL_DATA);
waveidtb_startDenoising(hObject);
% guidata(hObject, GUI_GLOBAL_DATA);

% --- Executes during object creation, after setting all properties.
function runDenoising_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runDenoising (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over runDenoising.
function runDenoising_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to runDenoising (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function ritLogo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tigerLogo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate tigerLogo
axes(hObject);
imshow('RIT_Logo.jpg');

% --- Executes during object creation, after setting all properties.
function cisLogo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cisLogo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate cisLogo
axes(hObject);
imshow('cis_logonotag_bw.jpg');

% --- Executes during object creation, after setting all properties.
function mrnLogo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mrnLogo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate mrnLogo
axes(hObject);
imshow('header.jpg');

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function cisLogo_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to cisLogo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on mouse press over axes background.
function mrnLogo_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to mrnLogo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('               Thank you for using WaveIDT ver1a                 ');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('If you have questions, please email skhullar@mrn.org or vcalhoun@mrn.org.');clear ans;
% Hint: delete(hObject) closes the figure
delete(hObject);clear GUI_GLOBAL_DATA handles eventdata;


% --- Executes during object deletion, before destroying properties.
function ritLogo_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to ritLogo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function ritLogo_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ritLogo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function cisLogo_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to cisLogo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function mrnLogo_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to mrnLogo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function uipanel3_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function uipanel3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
