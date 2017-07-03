function varargout = ISCanalysis(varargin)
% ISCANALYSIS M-file for ISCanalysis.fig
%      ISCANALYSIS, by itself, creates a new ISCANALYSIS or raises the existing
%      singleton*.
%
%      H = ISCANALYSIS returns the handle to a new ISCANALYSIS or the handle to
%      the existing singleton*.
%
%      ISCANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ISCANALYSIS.M with the given input arguments.
%
%      ISCANALYSIS('Property','Value',...) creates a new ISCANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ISCanalysis_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ISCanalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ISCanalysis

% Last Modified by GUIDE v2.5 12-Nov-2013 13:09:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ISCanalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @ISCanalysis_OutputFcn, ...
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


% --- Executes just before ISCanalysis is made visible.
function ISCanalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ISCanalysis (see VARARGIN)

% Choose default command line output for ISCanalysis
handles.output = hObject;
set(gcf,'CloseRequestFcn','ISCclosereq');

if isempty(varargin)
    handles = initISCanalysis(handles);
elseif length(varargin) == 1
    handles.Pub = varargin{1}.PublicParams;
    handles.validFlag = false;
    handles.Priv = [];
    handles.ParamsValid = false;
else
   error('Only single parameter file must be given as input!') 
end

handles = setParamFields(handles,handles.Pub);
%disable gridcomputing checkbox if windows pc is detected or no grid 
%(SGE/slurm) is detected
if(ispc || isempty(testGrid))
    set(handles.checkboxdisableGrid, 'Enable', 'off')
else
    set(handles.checkboxdisableGrid, 'Enable', 'on')
end

set(handles.figure1,'Name','ISC toolbox 2.0: Start-up GUI')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ISCanalysis wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ISCanalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
if ~isempty(handles)
    varargout{1} = handles.output;
    delete(handles.figure1)
end

% --- Executes on button press in checkboxFreq.
function checkboxFreq_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxFreq

if get(hObject,'Value')
    set(handles.pushbuttonFreqSett,'Enable','on')
    set(handles.editBands,'Enable','on')
%    set(handles.editTR,'Enable','on')
    set(handles.editFreqPerm,'Enable','on')
    handles.Pub.nrFreqBands = 3;
    set(handles.removeFilterMaps, 'Enable','on')
else
    set(handles.pushbuttonFreqSett,'Enable','off')    
    set(handles.editBands,'Enable','off')
%    set(handles.editTR,'Enable','off')
    set(handles.editFreqPerm,'Enable','off')
    handles.Pub.nrFreqBands = 0;
    set(handles.removeFilterMaps, 'Enable','off')
end
set(handles.editBands,'String',handles.Pub.nrFreqBands)
guidata(hObject, handles);

% --- Executes on button press in checkboxTime.
function checkboxTime_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.editWinLen,'Enable','on')
    set(handles.editWinStep,'Enable','on')
else
    set(handles.editWinLen,'Enable','off')
    set(handles.editWinStep,'Enable','off')
end

handles.Pub.winOn = get(hObject,'Value');
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of checkboxTime

%set(hObject,'String', str2double(get(hObject,'String'))


function editWinLen_Callback(hObject, eventdata, handles)
% hObject    handle to editWinLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editWinLen as text
%        str2double(get(hObject,'String')) returns contents of editWinLen as a double
val = ceil(str2double(get(hObject,'String')));
if isnan(val) || val <= 0 || isinf(val)
   set(hObject,'String',num2str(handles.Pub.windowSize))
else
    handles.Pub.windowSize = val;
    set(hObject,'String',num2str(val))
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editWinLen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWinLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editWinStep_Callback(hObject, eventdata, handles)
% hObject    handle to editWinStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editWinStep as text
%        str2double(get(hObject,'String')) returns contents of editWinStep as a double
val = ceil(str2double(get(hObject,'String')));
if isnan(val) || val <= 0 || isinf(val)
   set(hObject,'String',num2str(handles.Pub.windowStep))
else
    handles.Pub.windowStep = val;
    set(hObject,'String',num2str(val))
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editWinStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWinStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxCor.
function checkboxCor_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxCor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxCor
if get(hObject,'Value')
     set(handles.pushbuttonCorSett,'Enable','on')
else
    set(handles.pushbuttonCorSett,'Enable','off')     
end
handles.Pub.corOn = get(hObject,'Value');
% handles.Pub.corOn = 1;
% set(hObject,'Value',1);
% set(handles.pushbuttonCorSett,'Enable','on')
guidata(hObject, handles);

% --- Executes on button press in checkboxKen.
function checkboxKen_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxKen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxKen
handles.Pub.kenOn = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in checkboxMI.
function checkboxMI_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Pub.nmiOn = get(hObject,'Value');
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of checkboxMI


% --- Executes on button press in checkboxSSI.
function checkboxSSI_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxSSI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxSSI
handles.Pub.ssiOn = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in pushbuttonCorSett.
function pushbuttonCorSett_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCorSett (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
outpt = corrSettings(handles.Pub);
if(~isempty(outpt))
    handles.Pub = outpt;
end
guidata(hObject, handles);


function editTR_Callback(hObject, eventdata, handles)
% hObject    handle to editTR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTR as text
%        str2double(get(hObject,'String')) returns contents of editTR as a double

val = str2double(get(hObject,'String'));
if isnan(val) || val <= 0
   set(hObject,'String',num2str(1/handles.Pub.samplingFrequency))
else
    handles.Pub.samplingFrequency = 1/val;
    set(hObject,'String',num2str(val))
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editTR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editBands_Callback(hObject, eventdata, handles)
% hObject    handle to editBands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBands as text
%        str2double(get(hObject,'String')) returns contents of editBands as a double
val = ceil(str2double(get(hObject,'String')));
if isnan(val) || val < 2 || val > 15
   set(hObject,'String',num2str(handles.Pub.nrFreqBands))
else
    handles.Pub.nrFreqBands = val;
    checkLen(handles.Pub.nrFreqBands)
    set(hObject,'String',num2str(val))
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editBands_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonFreqSett.
function pushbuttonFreqSett_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFreqSett (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Pub = freqSettings(handles.Pub);
guidata(hObject, handles);


function editSubj_Callback(hObject, eventdata, handles)
% hObject    handle to editSubj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSubj as text
%        str2double(get(hObject,'String')) returns contents of editSubj as a double

handles = setSubjectBox(handles);

handles = validateParams(handles,'subj');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editSubj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSubj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editDestin_Callback(hObject, eventdata, handles)
% hObject    handle to editDestin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDestin as text
%        str2double(get(hObject,'String')) returns contents of editDestin as a double

tPath = get(hObject,'String');
handles.Pub.dataDestination  = tPath;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editDestin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDestin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTemplates_Callback(hObject, eventdata, handles)
% % hObject    handle to editTemplates (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of editTemplates as text
% %        str2double(get(hObject,'String')) returns contents of editTemplates as a double
% 
% tPath = get(hObject,'String');
% handles.Pub.atlasPath  = tPath;
% handles.Pub.maskPath  = tPath;
% handles = validateParams(handles,'template');
% guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editTemplates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTemplates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editProject_Callback(hObject, eventdata, handles)
% hObject    handle to editProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tPath = get(hObject,'String');
handles.Pub.dataDescription  = tPath;
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of editProject as text
%        str2double(get(hObject,'String')) returns contents of editProject as a double


% --- Executes during object creation, after setting all properties.
function editProject_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuSession.
function popupmenuSession_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuSession contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuSession
val = get(hObject,'Value');
if val <= size(handles.Pub.subjectSource,1)
    for h = 1:size(handles.Pub.subjectSource,2)
        D{h} = handles.Pub.subjectSource{val,h};
    end
    set(handles.editSubj,'String',D);
else
    set(handles.editSubj,'String',{});
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenuSession_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNii_Callback(hObject, eventdata, handles)
% hObject    handle to editNii (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNii as text
%        str2double(get(hObject,'String')) returns contents of editNii as a double


% --- Executes during object creation, after setting all properties.
function editNii_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNii (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonRun.
function pushbuttonRun_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Params.PublicParams = handles.Pub;
handles.Params.PrivateParams = handles.Priv;
handles.validFlag = true;

%handles.ParamsValid
%handles.Params
%isequal(handles.ParamsValid,handles.Params)

if ~isequal(handles.ParamsValid,handles.Params)
    disp('You must succesfully validate parameters before running the analysis!')
    return
end

if handles.validFlag
    guidata(hObject, handles);
    runAnalysis(handles.Params);
end

% --- Executes on selection change in popupmenuFormat.
function popupmenuFormat_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuFormat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuFormat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuFormat
switch get(hObject,'Value')
    case 1
        handles.Pub.fileFormatSubj = 'nii';
    case 2     
        handles.Pub.fileFormatSubj = 'mat';
    otherwise
        error('Unknown file format!')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenuFormat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuFormat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editFreqPerm_Callback(hObject, eventdata, handles)
% hObject    handle to editFreqPerm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFreqPerm as text
%        str2double(get(hObject,'String')) returns contents of editFreqPerm as a double

val = ceil(str2double(get(hObject,'String')));
if isnan(val) || val <= 0 || isinf(val)
   set(hObject,'String',num2str(handles.Pub.permutFreqComp))
else
    handles.Pub.permutFreqComp = val;
    set(hObject,'String',num2str(val))
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function editFreqPerm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFreqPerm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function menuFile_Callback(hObject, eventdata, handles)
% hObject    handle to menuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuOpen_Callback(hObject, eventdata, handles)
% hObject    handle to menuOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.mat','Open Parameter File');
if FileName ~= 0
    try 
        load(fullfile(PathName,FileName)) %load Params
        Params = testParamsMat(Params);
        handles.Pub = Params.PublicParams;
        handles.Priv = Params.PrivateParams;
        handles = setParamFields(handles,handles.Pub);
    catch
        disp(lasterr)
       disp('Invalid parameter file.')
       return
    end
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function menuSave_Callback(hObject, eventdata, handles)
% hObject    handle to menuSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uiputfile('*.mat','Save Parameter File');
if file ~= 0
   Params.PublicParams = handles.Pub;
   Params.PrivateParams = handles.Priv;
   save(fullfile(path,file),'Params')
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function menuExit_Callback(hObject, eventdata, handles)
% hObject    handle to menuExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

user_response = confCloseModal('Title','Confirm Exit');
switch lower(user_response)
    case 'no'
        % take no action
    case 'yes'
        handles.output = handles;
        guidata(hObject, handles)
        uiresume(handles.figure1)
        %        delete(handles.figure1)
end


% --- Executes on button press in pushbuttonValid.
function pushbuttonValid_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonValid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Params.PublicParams = handles.Pub;
Params.PrivateParams = handles.Priv;

if isfield(handles,'ParamsValid')
    if isequal(handles.ParamsValid,Params)
        disp(' ')
        disp('Parameters already succesfully validated!')
        return
    end
end

handles = validateDataAndParams(handles);

if handles.validFlag
    guidata(hObject, handles);
end

% --- Executes on button press in checkboxTemplate.
function checkboxTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.radiobuttonTemplateMNI,'Value',1,'Enable','on')
    handles.Pub.useTemplate = 1;
    set(handles.editMask,'Enable','on','String',handles.Pub.atlasPath)
%    set(handles.editMask,'Enable','off','String','standard')
    set(handles.textMask,'Enable','on','String','Directory of standard templates','HorizontalAlignment','Left')
else
    set(handles.radiobuttonTemplateMNI,'Value',0,'Enable','off')    
    handles.Pub.useTemplate = 0;    
    set(handles.editMask,'Enable','on','String',handles.Pub.atlasPath,'HorizontalAlignment','Left')
    set(handles.textMask,'Enable','on','String','Binary mask file name (extension .nii or .mat)')
end
guidata(hObject, handles);

% --- Executes on button press in radiobuttonTemplateMNI.
function radiobuttonTemplateMNI_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonTemplateMNI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'Value',1)

% --- Executes on button press in checkboxPhase.
function checkboxPhase_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxPhase

handles.Pub.calcPhase = get(hObject,'Value');
guidata(hObject, handles);


function editMask_Callback(hObject, eventdata, handles)
% hObject    handle to editMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMask as text
%        str2double(get(hObject,'String')) returns contents of editMask as a double
handles.Pub.atlasPath = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editMask_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkboxdisableGrid.
function checkboxdisableGrid_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxdisableGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxdisableGrid

handles.Pub.disableGrid = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in checkboxRemoveMaps.
function checkboxRemoveMaps_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxRemoveMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxRemoveMaps
handles.Pub.removeMemmaps = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in removeFilterMaps.
function removeFilterMaps_Callback(hObject, eventdata, handles)
% hObject    handle to removeFilterMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of removeFilterMaps
handles.Pub.removeFiltermaps = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in checkboxFreqComp.
%function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxFreqComp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxFreqComp


% --- Executes on button press in checkboxSessionComp.
function checkboxSessionComp_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxSessionComp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxSessionComp
val = get(hObject,'Value');
handles.Pub.sessionCompOn = val;
if val
    set(handles.editSessionPerm,'Enable','on')
else
    set(handles.editSessionPerm,'Enable','off')    
end
guidata(hObject, handles);



function editZPFsession_Callback(hObject, eventdata, handles)
% hObject    handle to editSessionPerm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSessionPerm as text
%        str2double(get(hObject,'String')) returns contents of editSessionPerm as a double

val = ceil(str2double(get(hObject,'String')));
if isnan(val) || val <= 0 || isinf(val)
   set(hObject,'String',num2str(handles.Pub.permutFreqComp))
else
    handles.Pub.permutFreqComp = val;
    set(hObject,'String',num2str(val))
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editZPFsession_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSessionPerm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSessionPerm_Callback(hObject, eventdata, handles)
% hObject    handle to editSessionPerm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSessionPerm as text
%        str2double(get(hObject,'String')) returns contents of editSessionPerm as a double

val = ceil(str2double(get(hObject,'String')));
if isnan(val) || val <= 0 || isinf(val)
   set(hObject,'String',num2str(handles.Pub.permutSessionComp))
else
    handles.Pub.permutSessionComp = val;
    set(hObject,'String',num2str(val))
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editSessionPerm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSessionPerm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function checkboxSessionComp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkboxSessionComp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in checkboxFreqComp.
function checkboxFreqComp_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxFreqComp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxFreqComp
val = get(hObject,'Value');
handles.Pub.freqCompOn = val;
if val
    set(handles.editFreqPerm,'Enable','on')
else
    set(handles.editFreqPerm,'Enable','off')    
end
guidata(hObject, handles);



% --- Executes on button press in pushbuttonExport.
function pushbuttonExport_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt = {'Enter variable name:'};
title = 'Assign parameters to workspace';
lines = 1;
def = {'Params'};
answer = inputdlg(prompt,title,lines,def);

%T = createImageData(handles);
Params.PublicParams = handles.Pub;
Params.PrivateParams = handles.Priv;
%T.info = handles.info;
if ~isempty(answer)
    assignin('base',answer{1},Params);
end
%guidata(hObject, handles);




% --- Executes on button press in pushbuttonGUI.
function pushbuttonGUI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Params.PublicParams = handles.Pub;
Params.PrivateParams = handles.Priv;
ISCtool(Params);
