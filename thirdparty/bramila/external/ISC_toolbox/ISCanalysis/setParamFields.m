function handles = setParamFields(handles,PublicParams)

if PublicParams.nrFreqBands > 0;
    set(handles.checkboxFreq,'Value',1)
else
    set(handles.checkboxFreq,'Value',0)
end

set(handles.editTR,'Enable','on')
    
% set frequency band setting:
if get(handles.checkboxFreq,'Value') == 1
    set(handles.pushbuttonFreqSett,'Enable','on')
    set(handles.editBands,'Enable','on')
    if get(handles.checkboxFreqComp,'Value')
        set(handles.editFreqPerm,'Enable','on')
    else
        set(handles.editFreqPerm,'Enable','off')        
    end
%    PublicParams.nrFreqBands = 3;
else
    set(handles.pushbuttonFreqSett,'Enable','off')
    set(handles.editFreqPerm,'Enable','off')
    set(handles.editTR,'Enable','off')
    set(handles.editBands,'Enable','off')
%    PublicParams.nrFreqBands = 0;
end

set(handles.editBands,'String',num2str(PublicParams.nrFreqBands));
checkLen(PublicParams.nrFreqBands);

set(handles.editFreqPerm,'String',num2str(PublicParams.permutFreqComp));
set(handles.checkboxFreqComp,'Value',PublicParams.freqCompOn);

set(handles.editSessionPerm,'String',num2str(PublicParams.permutSessionComp));
set(handles.checkboxSessionComp,'Value',PublicParams.sessionCompOn);

if get(handles.checkboxSessionComp,'Value')
    set(handles.editSessionPerm,'Enable','on')
else
    set(handles.editSessionPerm,'Enable','off')
end

if get(handles.checkboxFreqComp,'Value')
    set(handles.editFreqPerm,'Enable','on')
else
    set(handles.editFreqPerm,'Enable','off')
end


set(handles.editTR,'String',num2str(1/PublicParams.samplingFrequency));

% set time window settings:

set(handles.checkboxTime,'Value',PublicParams.winOn)
set(handles.editWinLen,'String',num2str(PublicParams.windowSize));
set(handles.editWinStep,'String',num2str(PublicParams.windowStep));
if get(handles.checkboxTime,'Value') == 1
    set(handles.editWinLen,'Enable','on')
    set(handles.editWinStep,'Enable','on')
else
    set(handles.editWinLen,'Enable','off')
    set(handles.editWinStep,'Enable','off')
end

% set template settings:
set(handles.radiobuttonTemplateMNI,'Value',1)
if handles.Pub.useTemplate == 1
    set(handles.checkboxTemplate,'Value',1)
    set(handles.radiobuttonTemplateMNI,'Enable','on')
else
    set(handles.checkboxTemplate,'Value',0)
    set(handles.radiobuttonTemplateMNI,'Enable','off')    
end


% set similarity metric settings:

set(handles.checkboxCor,'Value',PublicParams.corOn)
set(handles.checkboxKen,'Value',PublicParams.kenOn,'Visible','off')
set(handles.checkboxMI,'Value',PublicParams.nmiOn,'Visible','off')
set(handles.checkboxSSI,'Value',PublicParams.ssiOn,'Visible','off')

if get(handles.checkboxCor,'Value') == 1
    set(handles.pushbuttonCorSett,'Enable','on')
else
    set(handles.pushbuttonCorSett,'Enable','off')    
end

set(handles.checkboxPhase,'Value',PublicParams.calcPhase)

% Set grid computation parameters:
set(handles.checkboxdisableGrid,'Value',PublicParams.disableGrid)

% set subject source files:

set(handles.popupmenuSession,'String','Session 1','Value',1)
for k1 = 1:size(PublicParams.subjectSource,2)
    D{k1} = PublicParams.subjectSource{1,k1};
end
set(handles.editSubj,'String',D)
handles = setSubjectBox(handles);
handles = validateParams(handles,'subj');
% set paths:

set(handles.editDestin,'String',PublicParams.dataDestination);
set(handles.editTemplates,'String',PublicParams.atlasPath);
set(handles.editProject,'String',PublicParams.dataDescription);
set(handles.editTemplates,'Visible','off');
set(handles.textTemplates,'Visible','off');

%handles = validateParams(handles,'template');

if PublicParams.useTemplate
    set(handles.checkboxTemplate,'Value',1)
    set(handles.radiobuttonTemplateMNI,'Value',1,'Enable','on')
    set(handles.editMask,'Enable','on','String',PublicParams.atlasPath)
    set(handles.textMask,'Enable','on','String','Directory of standard templates')
%    set(handles.editTemplates,'Visible','on');
%    set(handles.textTemplates,'Visible','on');
else
    set(handles.checkboxTemplate,'Value',0)
    set(handles.radiobuttonTemplateMNI,'Value',0,'Enable','off')
    set(handles.editMask,'Enable','on','String',PublicParams.atlasPath)
    set(handles.textMask,'Enable','on','String','Binary mask file name (extension .nii or .mat)')
%    set(handles.editTemplates,'Visible','off');
%    set(handles.textTemplates,'Visible','off');
end

% set source file type:

if strcmp(PublicParams.fileFormatSubj,'nii')
    set(handles.popupmenuFormat,'Value',1);
elseif strcmp(PublicParams.fileFormatSubj,'mat')
    set(handles.popupmenuFormat,'Value',2);
else
    error('Unknown file format!!')
end

