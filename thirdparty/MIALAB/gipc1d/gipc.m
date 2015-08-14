function varargout = gipc(varargin)
% GIPC M-file for gipc.fig
%      GIPC, by itself, creates a new GIPC or raises the existing
%      singleton*.
%
%      H = GIPC returns the handle to a new GIPC or the handle to
%      the existing singleton*.
%
%      GIPC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GIPC.M with the given input arguments.
%
%      GIPC('Property','Value',...) creates a new GIPC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gipc_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gipc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gipc

% Last Modified by GUIDE v2.5 06-May-2008 14:18:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gipc_OpeningFcn, ...
                   'gui_OutputFcn',  @gipc_OutputFcn, ...
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

%% --- Executes just before gipc is made visible.
function gipc_OpeningFcn(hObject, eventdata, handles, varargin)
% 2 MATLAB AUTOMATED LINES
handles.output = hObject;
guidata(hObject, handles);

global gIpc
gipc_defaults;
ipctb_spm_defaults; %Load ipctb_spm_defaults for flip orientation of analyze images.

try %MRN MPAVO ce102107 fixes fonts, button color and other graphics
    %Unfortunatly it destroys UNIX ability to receive keyboard events.
    feature('JavaFigures', 0);
catch
end

% initiate form
set(handles.edQVal1Samp, 'String', num2str(gIpc.dFdrQ1Samp));
set(handles.edQVal2Samp, 'String', num2str(gIpc.dFdrQ2Samp));
set(handles.edPathExport, 'String', gIpc.sDirSave);

if gIpc.bFdr
    set(handles.la1Samp, 'String', 'FDR-q 1Samp:');
    set(handles.la2Samp, 'String', 'FDR-q 2Samp:')    
else
    set(handles.la1Samp, 'String', 'uncor-p 1Samp:');
    set(handles.la2Samp, 'String', 'uncor-p 2Samp:')        
end

if gIpc.bSkipStep1
    set(handles.pbPathGrp1, 'Enable', 'off');
    set(handles.pbPathGrp2, 'Enable', 'off');
    set(handles.stPathGrp1, 'String', '(browse off by gipc_defaults)');
    set(handles.stPathGrp2, 'String', '(browse off by gIpc.bSkipStep1)');    
else
    set(handles.pbPathGrp1, 'Enable', 'on');
    set(handles.pbPathGrp2, 'Enable', 'on');
    set(handles.stPathGrp1, 'String', 'None Selected');
    set(handles.stPathGrp2, 'String', 'None Selected');    
end

    
% % % % if ~bMinUi                                              %ce042408Takeaway
    % Render Logo Image
    h = handles.ipctb;
    image1 = imread('MindLogoIpc.png','png','BackgroundColor', [0 0 0]);
    handles.Image1 = image(image1,'Parent',handles.axLogo);
    set(handles.axLogo,'Visible', 'off');
% % % % else                                                    %ce042408Takeaway
% % % %     pbPathGrp1_Callback(hObject, eventdata, handles);   %ce042408Takeaway
% % % % %     load('/export/research/analysis/human/collaboration/olin/mialab/users/ceierud/gipc/gIpcUnix0417.mat');
% % % %     pbPathGrp2_Callback(hObject, eventdata, handles);   %ce042408Takeaway
% % % %     pbRun_Callback(hObject, eventdata, handles);        %ce042408Takeaway
% % % % end                                                     %ce042408Takeaway

% --- Outputs from this function are returned to the command line.
function varargout = gipc_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pbPathGrp1.
function pbPathGrp1_Callback(hObject, eventdata, handles)
    % Loads list of files manually for grp 1
    
    global gSave;
    [gSave.files1, designMatrix, gSave.numOfSub1, gSave.numOfSess1, dataSelMethod, gSave.diffTimePoints1, spmMatFlag] = ipctb_ica_dataSelection([], get(handles.stPathGrp1, 'String'), '', ['R_', 'I_'], 'real&imaginary');

    if length(deblank(gSave.files1(1,1).name(1,:))) <= 0
        % user never chose a file
        gSave.files1 = [];
        return
    end

    if gSave.numOfSub1 < 2
        msgbox(sprintf('At least 2 subjects are needed\n ---------------------------------------------------------------'));
        gSave.files1 = [];
        return
    end

    set(handles.stPathGrp1, 'String', [num2str(gSave.numOfSub1) ' subjects * ' num2str(gSave.numOfSess1) ' sessions selected']);
    
% --- Executes on button press in pbPathGrp2.
function pbPathGrp2_Callback(hObject, eventdata, handles)
    % Loads list of files manually for grp 2
    global gSave;
    [gSave.files2, designMatrix, gSave.numOfSub2, gSave.numOfSess2, dataSelMethod, gSave.diffTimePoints2, spmMatFlag] = ipctb_ica_dataSelection([], get(handles.stPathGrp1, 'String'), '', ['R_', 'I_'], 'real&imaginary');

    if length(deblank(gSave.files2(1,1).name(1,:))) <= 0
        % user never chose a file
        gSave.files2 = [];
        return
    end

    if gSave.numOfSub2 < 2
        msgbox(sprintf('At least 2 subjects are needed\n ---------------------------------------------------------------'));
        gSave.files2 = [];
        return
    end

    set(handles.stPathGrp2, 'String', [num2str(gSave.numOfSub2) ' subjects * ' num2str(gSave.numOfSess2) ' sessions selected']);    

%% --- Executes on button press in pbExit.
function pbExit_Callback(hObject, eventdata, handles)
% hObject    handle to pbExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    delete(get(0, 'children'));
    clear all;

%% --- Executes on button press in pbRun.
function pbRun_Callback(hObject, eventdata, handles)

    global gIpc gSave;
    
    % Create directory
    sDirOut = ipctb_backslash(get(handles.edPathExport, 'String'));

    if gIpc.bSkipStep1 % loads source file information
        % following load statment has to be at first (to be properly overwritten by more current values)
        load([sDirOut 'Settings.mat'], 'gSave');    
    end
    
    % Validation
    try, ntest = gSave.numOfSub1;
    catch, msgbox(sprintf('Have you imported Group 1?\n ---------------------------------------------------------------')); return;
    end
    try, ntest = gSave.numOfSub2;
    catch, msgbox(sprintf('Have you imported Group 2?\n ---------------------------------------------------------------'));return;
    end
    
    % Initialize variables based upon imported images
    gIpc.V = ipctb_spm_vol(gSave.files1(1,1).name(1, :));
    %These are the voxel dimensions of your fmri dataset.
    gIpc.xdim = gIpc.V.dim(1);
    gIpc.ydim = gIpc.V.dim(2);
    gIpc.zdim = gIpc.V.dim(3);
    gIpc.tdim = size(gSave.files1(1,1).name(:, 1),1);

    % ERROR VALIDATION STARTS HERE
    
    %Check onsets files
    for iSess = 1:2%gSave.numOfSess1 %loop through each session.      
        for iType = 1:size(gIpc.nSess(iSess).susOnset,1)
            if strcmp(gIpc.nSess(iSess).susOnset{iType,1}, gIpc.sOnsLockType)
                try, rodTrash = load([gIpc.sDirOnset gIpc.nSess(iSess).susOnset{iType,3}]);     
                catch, msgbox(sprintf('Are all onset files named correctly in gipc_defaults.m?\n ---------------------------------------------------------------'));return;
                end
            end
        end
    end         
    
    suTemp=ipctb_spm_vol(gIpc.sThreshMaskFile);
    if suTemp.dim~=gIpc.V.dim
        [sPath, sName, sExt, sVers] = fileparts(gIpc.sThreshMaskFile);
        msgbox(sprintf(['The mask file ' sName '.' sExt '\ndoes not match the voxels of your selected \nbrain images. Change mask or your input \ndata to match dimensions.\n ---------------------------------------------------------------']));
        return;
    end

    suTemp=ipctb_spm_vol(gIpc.sAnatTemp_dir);
    if suTemp.dim~=gIpc.V.dim
        [sPath, sName, sExt] = fileparts(gIpc.sAnatTemp_dir);
        msgbox(sprintf(['The anatomical file ' sName '.' sExt '\ndoes not match the voxels of your selected \nbrain images. Change mask or your input \ndata to match dimensions.\n ---------------------------------------------------------------']));
        return;
    end
    
    % error message in case of wrong nSplit value
    if mod(gIpc.xdim*gIpc.ydim*gIpc.zdim,gIpc.nSplit) ~= 0
        nNewSplit = 1;
        nNewSplit = nNewSplit + 1;
        while mod(gIpc.xdim*gIpc.ydim*gIpc.zdim,nNewSplit) ~= 0
            nNewSplit = nNewSplit + 1;
        end        

        sText = ['Brainvoxels are ' num2str(gIpc.xdim*gIpc.ydim*gIpc.zdim) ' and gIpc.nSplit = ' num2str(gIpc.nSplit) '.'];
        sText = [sText '\ngIpc.nSplit have to split the total number of brainvoxels evenly.'];
        sText = [sText '\nDo you want GIPC to automatically change gIpc.nSplit = ' num2str(nNewSplit) ' for this run or do you want to exit the program, change the gIpc.nSplit manually and restart GIPC from scratch?\n'];
        sButton = questdlg(sprintf(sText), ...
                       'GIPC', ...
                       'Automatically change gIpc.nSplit for this run','Exit and restart','Automatically change gIpc.nSplit for this run');        
        if strcmp(sButton,'Automatically change gIpc.nSplit for this run')
            gIpc.nSplit = nNewSplit;
            % just continue to the memory test
        else
            delete(get(0, 'children'));
            return
        end        
    end

    % Check memory and loop with new suggestions until user quits
    sRet = 'tryNewSplit';
    while strcmp(sRet,'tryNewSplit')
        sRet = fMemTest;
        if strcmp(sRet,'exit')
            delete(get(0, 'children'));
            return
        end
    end

    dRet = str2num(get(handles.edQVal1Samp, 'String'));
    if ~isempty(dRet) && isnumeric(dRet), gIpc.dFdrQ1Samp = dRet; else msgbox(sprintf('One sample FDR correction q-value/p-value has to be entered properly.\n ---------------------------------------------------------------')); return; end

    dRet = str2num(get(handles.edQVal2Samp, 'String'));
    if ~isempty(dRet) && isnumeric(dRet), gIpc.dFdrQ2Samp = dRet; else msgbox(sprintf('Two sample FDR correction q-value/p-value has to be entered properly.\n ---------------------------------------------------------------')); return; end

    if gIpc.nClusters>(size(gIpc.cmClust, 1)-1), msgbox(sprintf(['There are only colors for ' num2str(size(gIpc.cmClust, 1)-1) ' clusters.\n ---------------------------------------------------------------'])); return; end       
    
    % Create directories for files
    sPathIm = [sDirOut ipctb_backslash(gIpc.sImageDir)];
    bRet = ipctb_fCreateDir(sPathIm);
    if ~bRet % end if user wants to keep directory.
        return
    end        
    sPathGrp1 = [sDirOut 'Grp1Corr'];
    if ~gIpc.bSkipStep1 %creates grp directory
        bRet = ipctb_fCreateDir(sPathGrp1);
        if ~bRet % end if user wants to keep directory.
            return
        end    
    end
    sPathGrp2 = [sDirOut 'Grp2Corr'];   
    if ~gIpc.bSkipStep1 %creates grp directory
        bRet = ipctb_fCreateDir(sPathGrp2);
        if ~bRet % end if user wants to keep directory.
            return
        end  
    end
    % ERROR VALIDATION ENDS HERE

    if ~gIpc.bSkipStep1 %creates grp directory
        % create correlation map files
    % % % %     profile on;         %043008 take away
        bOk = ipctb_saveCorrMaps(sPathGrp1, gSave.files1, gSave.numOfSub1, gSave.numOfSess1);
    % % % %     p = profile('info');%043008 take away
    % % % %     save myprofiledata1 p;%043008 take away
    % % % %     clear p;            %043008 take away
    % % % %     profile off;        %043008 take away
    % % % %     profile on;         %043008 take away
        bOk = ipctb_saveCorrMaps(sPathGrp2, gSave.files2, gSave.numOfSub2, gSave.numOfSess2);
    % % % %     prof2 = profile('info');%043008 take away
    % % % %     save myprofiledata2 prof2;%043008 take away
    % % % %     clear prof2;        %043008 take away
    end
    disp(sprintf('\nAll correlation maps are written!'));
    
    [codMnG1 codSe2G1 dTMax1] = ipctb_ustatWrap(get(handles.edPathExport, 'String'), 'Grp1Corr', gIpc.sGrp1, gSave.numOfSub1);
    [codMnG2 codSe2G2 dTMax2] = ipctb_ustatWrap(get(handles.edPathExport, 'String'), 'Grp2Corr', gIpc.sGrp2, gSave.numOfSub2);
    disp(sprintf('\nAll U-statistics done!'));
    
    % Threshhold for 2 sample
    gIpc.dTMin2Sam = ipctb_P2T(gIpc.dFdrQ2Samp, gSave.numOfSub1 + gSave.numOfSub2 - 2);
    if gIpc.dTMin2Sam < 0, gIpc.dTMin2Sam = 0; end      
    
    % Two sided T-test
    if gIpc.b2SampT_Grp1Min2 %grp1-2
        codT2Samp = (codMnG1-codMnG2)./sqrt(  ((gSave.numOfSub1-1)*gSave.numOfSub1^2 + (gSave.numOfSub2-1)*gSave.numOfSub2^2) / (gSave.numOfSub1 + gSave.numOfSub2)  *  (1/gSave.numOfSub1 + 1/gSave.numOfSub2));   
    else%grp2-1
        codT2Samp = (codMnG2-codMnG1)./sqrt(  ((gSave.numOfSub1-1)*gSave.numOfSub1^2 + (gSave.numOfSub2-1)*gSave.numOfSub2^2) / (gSave.numOfSub1 + gSave.numOfSub2)  *  (1/gSave.numOfSub1 + 1/gSave.numOfSub2));   
    end
    codT2Raw = codT2Samp;
    % one dimensional t-values to 3D in accordance to brain
    codT2Raw=reshape(codT2Raw,gIpc.xdim, gIpc.ydim, gIpc.zdim);
    % save 2 sample raw T
    save([sDirOut gIpc.sFileT2Unthresh 'Extra.mat'], 'codT2Raw');
    % SAVE 2 sample raw T IN IMG-FORMAT (SPM AND NIFTI)
    gIpc.V.fname=[sDirOut gIpc.sFileT2Unthresh '.img'];
    sOldDir = pwd;
    cd(gIpc.sDirSpm);
    ipctb_spm_write_vol(gIpc.V, codT2Raw);
    cd(sOldDir);
    clear codT2Raw;
    
    %sort out neg T-values
    iPos = (codT2Samp>gIpc.dTMin2Sam);
    arTemp = zeros(gIpc.xdim*gIpc.ydim*gIpc.zdim,1);
    arTemp(iPos) = codT2Samp(iPos);
    codT2Samp = arTemp;
    clear arTemp;   

    %Mask for all brain data needed to discriminate nans
    arTemp=ipctb_spm_read_vols(ipctb_spm_vol(gIpc.sThreshMaskFile));
    arTemp(isfinite(arTemp) == 0) = 0;
    arTemp=reshape(arTemp,gIpc.xdim*gIpc.ydim*gIpc.zdim,1);
    iMask = (arTemp > mean(arTemp));
    
    ind = find(iPos & iMask);

    if gIpc.bFdr    
        % FDR correction
        [nSignif, iFdrMask] = ipctb_fdr(2*(1-ipctb_spm_Tcdf(codT2Samp(ind), gSave.numOfSub1+gSave.numOfSub2-2)), gIpc.dFdrQ2Samp);
        codFdrCorrected = zeros(gIpc.xdim*gIpc.ydim*gIpc.zdim, 1);
        if nSignif ~= 0
            codFdrCorrected(ind(iFdrMask), 1) = codT2Samp(ind(iFdrMask), 1);
        else
            disp('FDR sorted away all results')
        end
        codT2Samp = codFdrCorrected;
    end

    dTMax2Sam = max(codT2Samp);
    
    % one dimensional t-values to 3D in accordance to brain
    codT2Samp=reshape(codT2Samp,gIpc.xdim, gIpc.ydim, gIpc.zdim);
    % save 2 sample brain matrix
    save([sDirOut gIpc.sFileT2 'Extra.mat'], 'codT2Samp');
    % SAVE 2 sample brain matrix IN IMG-FORMAT (SPM AND NIFTI)
    gIpc.V.fname=[sDirOut gIpc.sFileT2 '.img'];
    sOldDir = pwd;
    cd(gIpc.sDirSpm);
    ipctb_spm_write_vol(gIpc.V, codT2Samp);
    cd(sOldDir);    
    disp(['Successfully saved two sample t-map ' gIpc.sFileT2]);

    if dTMax2Sam ~= 0
        % Run only if there is a T2 map
        % save the clusters for spm, slover or nifti
        bOk = ipctb_saveCluster(gIpc.sClust1, gSave.files1, sDirOut, gSave.numOfSub1, gSave.numOfSess1);
        bOk = ipctb_saveCluster(gIpc.sClust2, gSave.files2, sDirOut, gSave.numOfSub2, gSave.numOfSess2);
        bOk = ipctb_TimeLock([sDirOut gIpc.sClust1 '_tc.mat'], [sDirOut gIpc.sClust1 '_tcLock.mat']);
        bOk = ipctb_TimeLock([sDirOut gIpc.sClust2 '_tc.mat'], [sDirOut gIpc.sClust2 '_tcLock.mat']);    
        % put converse timecourses on the clustermaps as well
        bOk = ipctb_saveClusterTcConverse(gIpc.sClust1, gIpc.sClust2, sDirOut, gSave.numOfSess1);
        bOk = ipctb_saveClusterTcConverse(gIpc.sClust2, gIpc.sClust1, sDirOut, gSave.numOfSess2);   
        % Delete files only used for data processing
        delete([sDirOut gIpc.sClust1 'Copy.mat']);
        delete([sDirOut gIpc.sClust2 'Copy.mat']);
        bOk = ipctb_TimeLock([sDirOut gIpc.sClust1 '_tcConverse.mat'], [sDirOut gIpc.sClust1 '_tcLockConv.mat']);
        bOk = ipctb_TimeLock([sDirOut gIpc.sClust2 '_tcConverse.mat'], [sDirOut gIpc.sClust2 '_tcLockConv.mat']);    
    end
    disp(sprintf('\nClusters with timecourses done!'));

    dTMin1 = ipctb_P2T(gIpc.dFdrQ1Samp, gSave.numOfSub1-1);
    if dTMin1 < 0, dTMin1 = 0; end
    dTMin2 = ipctb_P2T(gIpc.dFdrQ1Samp, gSave.numOfSub2-1);
    if dTMin2 < 0, dTMin2 = 0; end

    % paint resulting image slices 
    sPathW = [sDirOut ipctb_backslash(gIpc.sImageDir)];
    hSlover = ipctb_spm_figure('Create','Graphics');%,'Graphics','on');%;%('Visible','on');
    
    % paint brains to *.png
    if dTMin1 > dTMax1, disp('Warning 1-sample have inverse thresholds'); end
    ipctb_sloverDisp('blob', [sDirOut gIpc.sGrp1 '.img'], [sPathW gIpc.sGrp1 '.png'], [dTMin1 dTMax1], hSlover);
    ipctb_sloverDisp('blob', [sDirOut gIpc.sGrp2 '.img'], [sPathW gIpc.sGrp2 '.png'], [dTMin2 dTMax2], hSlover);       
    if dTMax2Sam ~= 0
        % Run only if there is a T2 map
        if gIpc.dTMin2Sam > dTMax2Sam, disp('Warning 2-sample have inverse thresholds'); end
        ipctb_sloverDisp('blob', [sDirOut gIpc.sFileT2 '.img'], [sPathW gIpc.sFileT2 '.png'], [gIpc.dTMin2Sam dTMax2Sam], hSlover);
        ipctb_sloverDisp('clust', [sDirOut gIpc.sClust1 '.img'], [sPathW gIpc.sClust1 '.png'], [0 0], hSlover);
        ipctb_sloverDisp('clust', [sDirOut gIpc.sClust2 '.img'], [sPathW gIpc.sClust2 '.png'], [0 0], hSlover);
        % plot time courses to *.fig
        bOk = ipctb_TimeCourse([gIpc.sTitleGrp1 ', N = ' num2str(gSave.numOfSub1)], [sDirOut gIpc.sClust1 '_tcLock.mat'], [sPathW gIpc.sClust1 '_tcLock.fig']); 
        bOk = ipctb_TimeCourse([gIpc.sTitleGrp2 ' on clustermap of ' gIpc.sTitleGrp1 ', N = ' num2str(gSave.numOfSub2) ], [sDirOut gIpc.sClust1 '_tcLockConv.mat'], [sPathW gIpc.sClust1 '_tcLockConv.fig']);
        bOk = ipctb_TimeCourse([gIpc.sTitleGrp2 ', N = ' num2str(gSave.numOfSub2)], [sDirOut gIpc.sClust2 '_tcLock.mat'], [sPathW gIpc.sClust2 '_tcLock.fig']);
        bOk = ipctb_TimeCourse([gIpc.sTitleGrp1 ' on clustermap of ' gIpc.sTitleGrp2 ', N = ' num2str(gSave.numOfSub1)], [sDirOut gIpc.sClust2 '_tcLockConv.mat'], [sPathW gIpc.sClust2 '_tcLockConv.fig']);    
        disp(sprintf('\npng and fig files completed!'));
    else
        disp(sprintf('\nWarning no significant voxels were found for T2 map. Therefore no cluster files nor timecourse files were created.'));
    end
    close(hSlover);    
    
    save([sDirOut 'Settings.mat'], 'gIpc', 'gSave');
    msgbox('GIPC Completed!', 'GIPC');
    
%% Lets users browse for a directory
function pbPathOut_Callback(hObject, eventdata, handles)
try %Needed for save the dialog box to accept mouse clicks
    feature('JavaFigures', 1);
catch
end
vRet = uigetdir(get(handles.edPathExport, 'String'), 'GIPC choose directory');
if vRet ~= 0,
    set(handles.edPathExport, 'String', vRet);
end
try %MRN MPAVO ce102107 fixes fonts, button color and other graphics
    %Unfortunatly it destroys UNIX ability to receive keyboard events.
    feature('JavaFigures', 0);
catch
end

%% --- Executes on button press in pbAbout.
function pbAbout_Callback(hObject, eventdata, handles)
% hObject    handle to pbAbout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox(sprintf(['Interparticipant Correlation of fMRI Toolbox\n\n'...
'GIPC v1.0d\n\n'...
'Release Date: 05/12/2008\n\n'...
'Authors: Vince Calhoun, Dae Il Kim & Cyrus Eierud\n\n'...
'Organization: The Mind Research Network (www.mrn.org)\n\n'...
'Note: This IPC Toolbox uses code from SPM (Statistical Parametric Mapping)\n'...
'from (Wellcome Department of Imaging Neuroscience),\n'...
'GIFT from the Mind Research Network and\n'...
'and Piotr Dollars toolbox (http://vision.ucsd.edu/~pdollar).']));

% 'Website: http://gipc.sourceforge.net\n\n'...

%% memorytest
function sRet = fMemTest()
    global gIpc gSave;
    sRet = 'OkContinue';
    nGrpMax = max([gSave.numOfSub1 gSave.numOfSub2]);
    try
        adMemTest = zeros(prod([gIpc.xdim gIpc.ydim gIpc.zdim])/gIpc.nSplit, (nGrpMax^2-nGrpMax)/2   );
    catch
        clear adMemTest;        
        nNewSplit = gIpc.nSplit+1;
        while mod(gIpc.xdim*gIpc.ydim*gIpc.zdim,nNewSplit) ~= 0
            nNewSplit = nNewSplit + 1;
        end

        sText = ['Computer has not enough memory to run all chosen subjects with current gIpc.nSplit = ' num2str(gIpc.nSplit) '.'];
        sText = [sText '\nDo you want GIPC to automatically change gIpc.nSplit = ' num2str(nNewSplit) ' for this run or do you want to exit the program, change the gIpc.nSplit manually and restart GIPC from scratch?\n'];
        sButton = questdlg(sprintf(sText), ...
                       'GIPC', ...
                       'Automatically change gIpc.nSplit for this run','Exit and restart','Automatically change gIpc.nSplit for this run');        
        if strcmp(sButton,'Automatically change gIpc.nSplit for this run')
            gIpc.nSplit = nNewSplit;
            sRet = 'tryNewSplit';
            return;
        else
            sRet = 'exit';
            return;
        end
    end    
    clear adMemTest;