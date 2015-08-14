function varargout = LUI(varargin)
 % The Laterality User Interface (LUI) computes the lateral comparison
 % of nifti images.

try
  feature('JavaFigures', 0);
catch
end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LUI_OpeningFcn, ...
                   'gui_OutputFcn',  @LUI_OutputFcn, ...
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



% --- Executes just before LUI is made visible.
function LUI_OpeningFcn(hObject, eventdata, handles, varargin)

lui_spm_defaults;
%load to get left-right orientation

% Choose default command line pbPathOut for LUI
handles.pbPathOut = hObject;

% Update handles structure
guidata(hObject, handles);



%initialize input and pbPathOut tests
reset(hObject, eventdata, handles);

function reset(hObject, eventdata, handles)
%resets parameters for input
global isInputSet isOutputSet;
set(handles.edPathExport, 'String', 'Set Output Path');
set(handles.inputStatus, 'String', 'Status: None Selected');

isInputSet = false;
isOutputSet = false;

% --- Outputs from this function are returned to the command line.
function varargout = LUI_OutputFcn(hObject, eventdata, handles) 

% Get default command line pbPathOut from handles structure
varargout{1} = handles.pbPathOut;


% --- Executes on button press in input.
function input_Callback(hObject, eventdata, handles)

global isInputSet isOutputSet a;

a=lui_spm_select([1 inf],'image','Select input files...');
% a=lui_spm_select;

if a~=0
isInputSet = true;
set(handles.inputStatus, 'String', 'Status: Input Selected');
end
% --- Executes on button press in pbPathOut.
function pbPathOut_Callback(hObject, eventdata, handles)

global isInputSet isOutputSet;
% vRet = uigetdir(get(handles.edPathExport, 'String'), 'LUI choose directory');
vRet = lui_spm_select(1,'dir','Select Output Directory...');

if vRet ~= 0,
    set(handles.edPathExport, 'String', vRet);
    isOutputSet = true;
end

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)

global isInputSet isOutputSet a;

if isInputSet && isOutputSet

    %output path
    outPath = get(handles.edPathExport, 'String');
    % begin script    
    numfiles = size(a,1);

    %creating input folder file index list
        
    index = 0;
    prevPath = '';
    indexList = '';
    for j = 1:numfiles
        % if input file path is different than the previous, increment the index #
        [pathstr, name, ext, versn] = fileparts(a(j,:));
        if ~strcmp(pathstr,prevPath)
            
            index = index + 1;
            indexStr = num2str(index);
            
            while length(indexStr) < 3
                % prepending 0's for consistant string length
                indexStr = ['0',indexStr];
                
            end
            
            
            prevPath = pathstr;          
        end
        
        indexList = [indexList;indexStr];
    end
    
    for j = 1:numfiles,

      %GUI output status update
      set(handles.runStatus, 'String',['working on file #', num2str(j)]);
      drawnow;

       [V] = lui_spm_vol(a(j,:));
       data = lui_spm_read_vols(V);
       Lmsk = zeros(size(data));
       Rmsk = zeros(size(data));

       xdim = size(data,1);

       L=data(1:floor(xdim/2),:,:);
       R=data(end:-1:ceil(xdim/2)+1,:,:);

    %   data(1:floor(xdim/2),:,:) = (R-L)./(R+L+eps);
    %   data(1:floor(xdim/2),:,:) = atan(L./(R+eps))*180/pi;
       data(1:floor(xdim/2),:,:) = (L-R);

       data(end:-1:ceil(xdim/2)+1,:,:)=(R-L);

       if (xdim/2)~=round(xdim/2),
          data(floor(xdim/2)+1,:,:)=data(floor(xdim/2),:,:);
       end;

       [dr,name,ext,versn]=fileparts(a(j,:));
       
       % strip spm commas, etc... off of extensions
       ext=ext(1:4);
       
       fn = [name ext versn];

       V.fname = fullfile(outPath,['lat_' indexList(j,:) '_' fn]);
       lui_spm_write_vol(V,data);

       if (j==1),
          Lavg = L/numfiles;
          Ravg = R/numfiles;
       else,
          Lavg = Lavg+L/numfiles;
          Ravg = Ravg+R/numfiles;
       end;

       Lmsk(1:floor(xdim/2),:,:) = L;
       Rmsk(1:floor(xdim/2),:,:) = R;
       Lmsk(end:-1:ceil(xdim/2)+1,:,:)=L;
       Rmsk(end:-1:ceil(xdim/2)+1,:,:)=R;

       V.fname = fullfile(outPath,['L_' indexList(j,:) '_' fn]);
      %% V.dim(4)=16;
       lui_spm_write_vol(V,Lmsk);

       V.fname = fullfile(outPath,['R_' indexList(j,:) '_' fn]);
      %% V.dim(4)=16;
       lui_spm_write_vol(V,Rmsk);

    end;


    Lmsk(1:floor(xdim/2),:,:) = Lavg;
    Rmsk(1:floor(xdim/2),:,:) = Ravg;
    Lmsk(end:-1:ceil(xdim/2)+1,:,:)=Lavg;
    Rmsk(end:-1:ceil(xdim/2)+1,:,:)=Ravg;

    %masks
    V.fname = fullfile(outPath,['mask_Rp_Lp_' fn]);
    %V.dim(4)=16;
    lui_spm_write_vol(V,1*((Rmsk>0)&(Lmsk>0)));

    V.fname = fullfile(outPath,['mask_Rp_Ln_' fn]);
    %V.dim(4)=16;
    lui_spm_write_vol(V,1*((Rmsk>0)&(Lmsk<0)));

    V.fname = fullfile(outPath,['mask_Rn_Lp_' fn]);
    %V.dim(4)=16;
    lui_spm_write_vol(V,1*((Rmsk<0)&(Lmsk>0)));

    V.fname = fullfile(outPath,['mask_Rn_Ln_' fn]);
    %V.dim(4)=16;
    lui_spm_write_vol(V,1*((Rmsk<0)&(Lmsk<0)));

    V.fname = fullfile(outPath,['Lmean_' fn]);
    %V.dim(4)=16;
    lui_spm_write_vol(V,Lmsk);

    V.fname = fullfile(outPath,['Rmean_' fn]);
    %V.dim(4)=16;
    lui_spm_write_vol(V,Rmsk);

    %end script
    
    set(handles.runStatus, 'String','Computation Complete');
    %resetting parameters for next execution
    reset(hObject, eventdata, handles);
    
    drawnow;
   
else
    set(handles.runStatus,'String','Both input and output must be set');
end





function edPathExport_Callback(hObject, eventdata, handles)

%Output path Manually typed in
global isOutputSet;
if exist(get(handles.edPathExport, 'String'))
    isOutputSet = true;
else
    set(handles.edPathExport, 'String', 'Output Path not found');
    isOutputSet = false;
end


% --- Executes during object creation, after setting all properties.
function edPathExport_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


