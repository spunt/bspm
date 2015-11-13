function xjview(varargin)
% xjview, version 8.13.1
%
% usage 1: xjview (no argument)
%           for displaying a result img file, or multiple image files,
%           (which will be loaded later) and changing p-value or t/f-value
% usage 2: xjview(imagefilename)
%           for displaying the result img file and changing p-value or
%           t-value
%           Example: xjView spmT_0002.img
%                    xjView('spmT_0002.img')
%                    xjView mymask.img
% usage 3: xjview(imagefilename1, imagefilename2, ...)
%           for displaying the result img files and changing p-value or
%           t/f-value
%           Example: xjView spmT_0002.img spmT_0005.img spmT_0007.img
%                    xjView('spmT_0002.img', 'spmT_0003.img', 'spmT_0006.img')
%                    xjView myMask1.img myMask2.img myMask3.img
% usage 4: xjview(mnicoord, intensity)
%           for displaying where are the mni coordinates
%           mnicoord: Nx3 matrix of mni coordinates
%           intensity: (optional) Nx1 matrix, usually t values of the
%           corresponding voxels
%           Example: xjView([20 10 1; -10 2 5],[1;2])
%                    xjView([20 10 1; -10 2 5])
%
% http://www.alivelearn.net/xjview
%
% by Xu Cui, Jian Li and Xiaowei Song
% http://www.alivelearn.net/xjview
% Xu Cui's personal blog: http://www.alivelearn.net
%
% last modified: 11/12/2015 at 16:30PM (fix stat bug)
% last modified: 11/09/2015 at 11:09AM (allow to change the minimum of color bar range, allow to resize window)
% last modified: 03/31/2014 (fix setting F-test df bug)
% last modified: 10/23/2012 (fix colorbar max bug)
% last modified: 7/4/2012 (fix overylay bug)
% last modified: 4/18/2012 (choice of color, set p-value individually)
% last modified: 10/21/2011 (compatible with MatLab 2011)
% last modified: 2/28/2011 (add small volume correction)
% last modified: 2/25/2011 (fix volume bug under SPM 8 r4010)
% last modified: 12/14/2009 (FDR and Feed)
% last modified: 06/24/2009 (fix bug in activation report)
% last modified: 05/24/2009 (slice view)
% last modified: 04/17/2009 (SPM8 compatible)
% last modified: 11/20/2007 (new database, all database from wfu_pickatlas, including aal
% last modified: 06/01/2007 (correct FDR corrected p-value list, change intensity to handles.intensity{1})
% last modified: 02/18/2007 (add colorbar max control)
% last modified: 11/16/2006 (keyboard shortcut for open image and open roi file)
% last modified: 06/16/2006 (spm5 compatible)
% last modified: 05/30/2006 (left/right flip, path of mask.img and templateFile.img)
% last modified: 05/08/2006 (debug CallBack_volumePush function, change handles.intensity{1} to intensity)
% last modified: 04/03/2006 (modify tr)
% last modified: 12/28/2005 (modify SPM process)
%
% Thank Sergey Pakhomov for sharing his database (MNI Space Utility).
% Thank Joseph Maldjian for WFU_PickAtlas
% Thank Yuval Cohen for the maximize figure function (maximize.m)
%

warnstate = warning;
warning off;

% pre-set values
% important! you need compare the display of xjview and spm. If you find
% xjview flipped the left/right, you need to set leftrightflip = 1;
% otherwise leave it to 0.
leftrightflip = 0; 

% You only need to change M and DIM when you want to use xjview under
% 'usage 4'.
M = ...
    [-4     0     0    84;...
     0     4     0  -116;...
     0     0     4   -56;...
     0     0     0     1];
DIM = [41 48 35]';
 M = [
     -2     0     0    92
      0     2     0  -128
      0     0     2   -74
      0     0     0     1];
 DIM = [91   109    91]';
M = ...
    [-4     0     0    82;...
     0     4     0  -116;...
     0     0     4   -54;...
     0     0     0     1];
DIM = [40 48 34]';

TR = 2; % you don't need to set TR is you only use the viewing part of xjview
XJVIEWURL = 'http://www.alivelearn.net/xjview';

% system settings
try
    spmdir = spm('dir');
    %spm('defaults', 'fmri');
catch
    disp('Please add spm path.');
    warning(warnstate(1).state);
    return
end

try
    spm('defaults', 'fmri');
catch
	[];
end


if ispc
    os = 'windows';
elseif isunix
    os = 'linux';
else
    warndlg('I don''t know what kind of computer you are using. I assumed it is unix.', 'What computer are you using?');
    os = 'linux';
end
screenResolution = get(0,'ScreenSize');

xjviewpath = fileparts(which('xjview'));

% pre-set values
pValue = 0.001;
intensityThreshold = 0;
clusterSizeThreshold = 5;


% Appearance Settings
figurePosition =                    [0.100,   0.050,    0.550,    0.880];
sectionViewPosition =               [0.5,0.61,0.45,0.45];
glassViewAxesPosition =             [0.000,   0.600,    0.464,    0.400];
if screenResolution(3) <= 1024
    figurePosition =                    [0.100,   0.050,    0.700,    0.900];
    sectionViewPosition =               [0.5,       0.61,   0.45,   0.46];
    glassViewAxesPosition =             [0.000,   0.600,    0.464,    0.400];    
end

left =                              0.01;
editBoxHeight =                     0.05;
editBoxWidth =                      0.200;
editBoxLeft =                       0.100;

controlPanelPosition =              [left,     0.080,  0.500,      0.500];
stretchMatrix =                     diag([controlPanelPosition(3),controlPanelPosition(4),controlPanelPosition(3),controlPanelPosition(4)]);
controlPanelOffset =                controlPanelPosition' .* [1,1,0,0]';
heightUnit  =                       0.055;

sliderPosition =                    stretchMatrix*[0.000,   0*heightUnit-0.01,    1.000,           editBoxHeight]' + controlPanelOffset;
pValueTextPosition =                stretchMatrix*[0.000,   0.8*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
pValueEditPosition =                stretchMatrix*[0.100,   1*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
fdrTextPosition =                   stretchMatrix*[0.350,   0.8*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
fdrEditPosition =                   stretchMatrix*[0.450,   1*heightUnit,    editBoxWidth/2,    editBoxHeight]' + controlPanelOffset;
fdrPushPosition     =               stretchMatrix*[0.300,   1*heightUnit,    editBoxWidth/2,    editBoxHeight]' + controlPanelOffset;
intensityThresholdTextPosition =    stretchMatrix*[0.55,   0.8*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
intensityThresholdEditPosition =    stretchMatrix*[0.65,   1*heightUnit,    editBoxWidth*3/4,    editBoxHeight]' + controlPanelOffset;
dfTextPosition =                    stretchMatrix*[0.840,   0.8*heightUnit,    0.8-0.74,    editBoxHeight]' + controlPanelOffset;
dfEditPosition =                    stretchMatrix*[0.900,   1*heightUnit,    editBoxWidth/2,    editBoxHeight]' + controlPanelOffset;
clusterSizeThresholdTextPosition =  stretchMatrix*[0.000,   1.8*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
clusterSizeThresholdEditPosition =  stretchMatrix*[0.150,   2*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
pickThisClusterPushPosition =       stretchMatrix*[0.400,   2*heightUnit,    editBoxWidth*1,    editBoxHeight]' + controlPanelOffset;
selectThisClusterPushPosition =     stretchMatrix*[0.600,   2*heightUnit,    editBoxWidth*1,    editBoxHeight]' + controlPanelOffset;
clearSelectedClusterPushPosition =   stretchMatrix*[0.800,   2*heightUnit,    editBoxWidth*1,    editBoxHeight]' + controlPanelOffset;
thisClusterSizeTextPosition =       stretchMatrix*[0.600,   2*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
thisClusterSizeEditPosition =       stretchMatrix*[0.800,   2*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
loadImagePushPosition =             stretchMatrix*[0.000,   3*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
imageFileEditPosition =             stretchMatrix*[0.200,   3*heightUnit,    1-editBoxWidth,    editBoxHeight]' + controlPanelOffset;
saveImagePushPosition =             stretchMatrix*[0.000,   4*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
saveImageFileEditPosition =         stretchMatrix*[0.200,   4*heightUnit,    1-editBoxWidth,    editBoxHeight]' + controlPanelOffset;
saveResultPSPushPosition =          stretchMatrix*[0.000,   5*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
resultPSFileEditPosition =          stretchMatrix*[0.200,   5*heightUnit,    1-editBoxWidth,    editBoxHeight]' + controlPanelOffset;
displayIntensityTextPosition =      stretchMatrix*[0.000,   3*heightUnit,    editBoxWidth+0.1,    editBoxHeight]' + controlPanelOffset;
allIntensityRadioPosition =         stretchMatrix*[0.250,   3*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
positiveIntensityRadioPosition =    stretchMatrix*[0.400,   3*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
negativeIntensityRadioPosition =    stretchMatrix*[0.550,   3*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
renderViewCheckPosition =           stretchMatrix*[0.730,   3*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
renderStylePopPosition =           stretchMatrix*[0.90,   3*heightUnit,    editBoxWidth/2,    editBoxHeight]' + controlPanelOffset;
hideControlPushPosition =           stretchMatrix*[0.200,   4*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
reportPushPosition =                stretchMatrix*[0.000,   4*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
volumePushPosition =                stretchMatrix*[0.200,   4*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
smallVolumePosition =          stretchMatrix*[0.400,   4*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
commonRegionPushPosition =          stretchMatrix*[0.600,   4*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;

reslicePushPosition =               stretchMatrix*[0.600,   4*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
sliceViewPushPosition =             stretchMatrix*[0.780,   4*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
displayPushPosition =               stretchMatrix*[0.200,   4*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
allinonePushPosition =              stretchMatrix*[0.400,   4*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
searchPushPosition =                stretchMatrix*[0.000,   6*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
searchContentEditPosition =         stretchMatrix*[0.200,   6*heightUnit,    editBoxWidth*2,    editBoxHeight]' + controlPanelOffset;
searchTextPosition =                stretchMatrix*[0.600,   6*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
searchEnginePopPosition =           stretchMatrix*[0.600,   6*heightUnit,    1-0.6,    editBoxHeight]' + controlPanelOffset;
overlayPushPosition =               stretchMatrix*[0.000,   5*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
overlayEditPosition =               stretchMatrix*[0.200,   5*heightUnit,    0.6-editBoxWidth,    editBoxHeight]' + controlPanelOffset;
overlayPopPosition =                stretchMatrix*[0.600,   5*heightUnit,    1-0.6,    editBoxHeight]' + controlPanelOffset;
helpPosition =                      stretchMatrix*[0.800,   5*heightUnit,    editBoxWidth,    editBoxHeight]' + controlPanelOffset;
infoTextBoxPosition =               stretchMatrix*[0.000,   8*heightUnit,    1,           editBoxHeight*9]' + controlPanelOffset;
%xjViewPosition =                    stretchMatrix*[0.400,   13*heightUnit,    editBoxWidth*2.5,    editBoxHeight*3]' + controlPanelOffset;

sectionViewListboxPosition =        [sectionViewPosition(1)+0.4, sectionViewPosition(2)+0.02, 0.1, 0.14];
sectionViewMoreTargetPushPosition = [sectionViewListboxPosition(1),sectionViewListboxPosition(2)-0.02,0.10,0.02];
xHairCheckPosition =                [sectionViewListboxPosition(1),sectionViewListboxPosition(2)+0.14,0.15,0.02];
setTRangeEditPosition =             [sectionViewListboxPosition(1),sectionViewListboxPosition(2)-0.06,0.05,0.02];
setTRangeEditPosition2 =             [sectionViewListboxPosition(1)+0.05,sectionViewListboxPosition(2)-0.06,0.05,0.02];
setTRangeTextPosition =             [sectionViewListboxPosition(1),sectionViewListboxPosition(2)-0.04,0.10,0.02];

getStructurePushPosition =          [glassViewAxesPosition(1), glassViewAxesPosition(2)-0.06, editBoxWidth/2, editBoxHeight/2];
structureEditPosition =             [getStructurePushPosition(1), getStructurePushPosition(2), 1, getStructurePushPosition(4)];
framePosition =                     (controlPanelPosition - controlPanelOffset')*1.05 + 0.95*controlPanelOffset';

% draw figure and control
figureBKcolor=[176/255 252/255 188/255];
figureBKcolor=get(0,'Defaultuicontrolbackgroundcolor');
f = figure('unit','normalized','position',figurePosition,'Color',figureBKcolor,'defaultuicontrolBackgroundColor', figureBKcolor,...
        'Name','xjView', 'NumberTitle','off','resize','on','CloseRequestFcn', {@CallBack_quit, warnstate(1).state}, 'visible','off', 'DoubleBuffer','on');
handles = guihandles(f);
 
% databases
try
    X = load('TDdatabase');
    handles.DB = X.DB;
    handles.wholeMaskMNIAll = X.wholeMaskMNIAll;
catch
    errordlg('I can''t find TDdatabase.mat','TDdatabase not found');
end

handles.figure = f;
handles.instructionText = uicontrol(handles.figure,'style','text',...
        'unit','normalized','position',[0.02 0.7 0.4 0.2],...
        'string','Click menu File | Open Images','horizontal','left', 'fontSize',40);     
handles.frame = uicontrol(handles.figure,'style','frame',... 
        'unit','normalized',...
        'position',framePosition,...
        'Visible','off');  
handles.slider = uicontrol(handles.figure,'style','slider',... 
        'unit','normalized',...
        'position',sliderPosition,...
        'max',1,'min',0,...
        'sliderstep',[0.01,0.10],...
        'callback',@CallBack_slider,...
        'value',0,'Visible','on');
handles.pValueTextPosition = uicontrol(handles.figure,'style','text',...
        'unit','normalized','position',pValueTextPosition,...
        'string','pValue=','horizontal','left');    
handles.pValueEdit = uicontrol(handles.figure,'style','edit',...
        'unit','normalized','position',pValueEditPosition,...
        'horizontal','left',...
        'String', num2str(pValue),...
        'BackgroundColor', 'w',...
        'callback',@CallBack_pValueEdit);
handles.fdrTextPosition = uicontrol(handles.figure,'style','text',...
        'unit','normalized','position',fdrTextPosition,...
        'string','FDR p=','horizontal','left');    
handles.fdrEdit = uicontrol(handles.figure,'style','edit',...
        'unit','normalized','position',fdrEditPosition,...
        'horizontal','left',...
        'String', '',...
        'BackgroundColor', 'w',...
        'callback',@CallBack_fdrEdit);    
handles.intensityThresholdText = uicontrol(handles.figure,'style','text',...
        'unit','normalized','position',intensityThresholdTextPosition,...
        'string',' intensity=','horizontal','left');        
handles.intensityThresholdEdit = uicontrol(handles.figure,'style','edit',...
        'unit','normalized','position',intensityThresholdEditPosition,...
        'horizontal','left',...
        'BackgroundColor', 'w',...
        'String', num2str(intensityThreshold),...
        'callback',@CallBack_intensityThresholdEdit);
handles.dfText = uicontrol(handles.figure,'style','text',...
        'unit','normalized','position',dfTextPosition,...
        'string','df= ','horizontal','right');        
handles.dfEdit = uicontrol(handles.figure,'style','edit',...
        'unit','normalized','position',dfEditPosition,...
        'horizontal','left',...
        'BackgroundColor', 'w',...
        'String', '',...
        'callback',@CallBack_dfEdit);    
handles.clusterSizeThresholdText = uicontrol(handles.figure,'style','text',...
        'unit','normalized','position',clusterSizeThresholdTextPosition,...
        'string','cluster size >=','horizontal','left');
handles.clusterSizeThresholdEdit = uicontrol(handles.figure,'style','edit',...
        'unit','normalized','position',clusterSizeThresholdEditPosition,...
        'horizontal','left',...
        'BackgroundColor', 'w',...
        'String', num2str(clusterSizeThreshold),...
        'callback',@CallBack_clusterSizeThresholdEdit);
handles.thisClusterSizeText = uicontrol(handles.figure,'style','text',...
        'unit','normalized','position',thisClusterSizeTextPosition,...
        'string','size= ','horizontal','right', 'visible', 'off');
handles.thisClusterSizeEdit = uicontrol(handles.figure,'style','edit',...
        'unit','normalized','position',thisClusterSizeEditPosition,...
        'horizontal','left',...
        'Enable', 'inactive',...        
        'String', '','visible','off');   

handles.imageFileEdit = uicontrol(handles.figure,'style','edit',...
        'unit','normalized','position',imageFileEditPosition,...
        'horizontal','left',...
        'String', '',...
        'BackgroundColor', 'w',...
        'callback',@CallBack_imageFileEdit,...
        'visible','off');
handles.saveImageFileEdit = uicontrol(handles.figure,'style','edit',...
        'unit','normalized','position',saveImageFileEditPosition,...
        'horizontal','left',...
        'BackgroundColor', 'w',...
        'String', 'myMask.img',...
        'callback',@CallBack_saveImageFileEdit,...
        'visible','off');    
handles.saveResultPSEdit = uicontrol(handles.figure,'style','edit',...
        'unit','normalized','position',resultPSFileEditPosition,...
        'horizontal','left',...
        'BackgroundColor', 'w',...
        'String', 'myResult.ps',...
        'callback',@CallBack_saveResultPSEdit,...
        'visible','off');    
handles.loadImagePush = uicontrol(handles.figure,'style','push',...
        'unit','normalized','position',loadImagePushPosition,...
        'string','Load Image','callback',@CallBack_loadImagePush,...
        'visible','off');    
handles.saveImagePush = uicontrol(handles.figure,'style','push',...
        'unit','normalized','position',saveImagePushPosition,...
        'string','Save Image','callback',@CallBack_saveImagePush,...
        'visible','off');  
handles.saveResultPSPush = uicontrol(handles.figure,'style','push',...
        'unit','normalized','position',saveResultPSPushPosition,...
        'string','Save Result','callback',@CallBack_saveResultPSPush,...
        'visible','off');      
handles.getStructurePush = uicontrol(handles.figure,'style','push',...
        'unit','normalized','position',getStructurePushPosition,...
        'string','Get Structure','callback',@CallBack_getStructurePush,'visible','off');      
handles.structureEdit = uicontrol(handles.figure,'style','edit',...
        'unit','normalized','position',structureEditPosition,...
        'horizontal','center',...
        'enable', 'on',...
        'UserData',struct(...
		'hReg',	[],...
		'M',	M,...
		'D',	DIM,...
		'xyz',	[0 0 0]	));
handles.pickThisClusterPush = uicontrol(handles.figure,'style','push',...
        'unit','normalized','position',pickThisClusterPushPosition,...
        'string','Pick Cluster/Info','callback',@CallBack_pickThisClusterPush);          
handles.selectThisClusterPush = uicontrol(handles.figure,'style','push',...
        'unit','normalized','position',selectThisClusterPushPosition,...
        'string','Select Cluster','callback',@CallBack_selectThisClusterPush);          
handles.clearSelectedClusterPush = uicontrol(handles.figure,'style','push',...
        'unit','normalized','position',clearSelectedClusterPushPosition,...
        'string','Clear Selection','callback',@CallBack_clearSelectedClusterPush);          
    
handles.displayIntensityText = uicontrol(handles.figure,'style','text',...
        'unit','normalized','position',displayIntensityTextPosition,...
        'string','display intensity','horizontal','left'); 
handles.allIntensityRadio = uicontrol(handles.figure,'style','radio',...
        'unit','normalized',...
        'string','All',...
        'position',allIntensityRadioPosition,...
        'value', 1,...
        'callback',@CallBack_allIntensityRadio);
handles.positiveIntensityRadio = uicontrol(handles.figure,'style','radio',...
        'unit','normalized',...
        'string','Only +',...
        'position',positiveIntensityRadioPosition,...
        'callback',{@CallBack_allIntensityRadio,'+'});
handles.negativeIntensityRadio = uicontrol(handles.figure,'style','radio',...
        'unit','normalized',...
        'string','Only -',...
        'position',negativeIntensityRadioPosition,...
        'callback',{@CallBack_allIntensityRadio,'-'});
handles.renderViewCheck = uicontrol(handles.figure,'style','checkbox',... 
        'unit','normalized',...
        'string','Render View' ,...
        'horizontal', 'right',...
        'position',renderViewCheckPosition,...
        'callback', @CallBack_renderViewCheck);   
handles.renderStylePop = uicontrol(...
		'Units','normalized', ...
		'ListboxTop',0, ...
		'Position',renderStylePopPosition, ...
		'String',{'new';'old'}, ...
		'Style','popupmenu', ...
		'value',1,...
        'callback', @CallBack_renderStylePop);     
handles.sectionViewListbox = uicontrol(handles.figure,'style','listbox',...
        'unit','normalized',...
        'String', {'single T1','avg152PD','avg152T1','avg152T2','avg305T1','ch2','ch2bet','aal','brodmann'}, ...
        'value',3,...
        'position',sectionViewListboxPosition,...
        'callback',@CallBack_sectionViewListbox);
    
handles.xHairCheck = uicontrol(handles.figure,'style','checkbox',... 
        'unit','normalized',...
        'string','XHairs Off' ,...
        'horizontal','left',...
        'position',xHairCheckPosition,...
        'callback',@CallBack_xHairCheck);  
handles.sectionViewMoreTargetPush = uicontrol(handles.figure,'style','push',...
        'unit','normalized','position',sectionViewMoreTargetPushPosition,...
        'string','other ...','callback',@CallBack_sectionViewMoreTargetPush);   
handles.setTRangeEdit = uicontrol(handles.figure,'style','edit',...
        'unit','normalized','position',setTRangeEditPosition,'BackgroundColor', 'w',...
        'string','auto','callback',@CallBack_setTRangeEdit);
handles.setTRangeEdit = uicontrol(handles.figure,'style','edit',...
        'unit','normalized','position',setTRangeEditPosition2,'BackgroundColor', 'w',...
        'string','auto','callback',@CallBack_setTRangeEdit2);    
handles.setTRangeText = uicontrol(handles.figure,'style','text',...
        'unit','normalized','position',setTRangeTextPosition,...
        'string','colorbar max&min');    
handles.hideControlPush = uicontrol(handles.figure, 'style', 'push',...
        'unit','normalized',...
        'String', '<', 'position', hideControlPushPosition,...
        'visible','off');
handles.reportPush = uicontrol(handles.figure, 'style', 'push',...
        'unit','normalized',...
        'String', 'report', ...
        'position', reportPushPosition,...
        'callback', @CallBack_reportPush);    
handles.volumePush = uicontrol(handles.figure, 'style', 'push',...
        'unit','normalized',...
        'String', 'volume', ...
        'position', volumePushPosition,...
        'callback', @CallBack_volumePush);    
handles.commonRegionPush = uicontrol(handles.figure, 'style', 'push',...
        'unit','normalized',...
        'String', 'common region', ...
        'position', commonRegionPushPosition,...
        'callback', @CallBack_commonRegionPush);       
% handles.fdrPush = uicontrol(handles.figure, 'style', 'push',...
%         'unit','normalized',...
%         'String', 'FDR', ...
%         'position', fdrPushPosition,...
%         'callback', @CallBack_fdrPush);      
handles.reslicePush = uicontrol(handles.figure, 'style', 'push',...
        'unit','normalized',...
        'String', 'reslice', ...
        'visible', 'off', ...
        'position', reslicePushPosition,...
        'callback', @CallBack_reslicePush);           

handles.sliceViewPush = uicontrol(handles.figure, 'style', 'push',...
        'unit','normalized',...
        'String', 'slice view', ...
        'position', sliceViewPushPosition,...
        'callback', @CallBack_sliceViewPush);       
    
handles.displayPush = uicontrol(handles.figure, 'style', 'push',...
        'unit','normalized',...
        'String', 'display', ...
        'position', displayPushPosition,...
        'callback', @CallBack_displayPush,...
        'visible','off');    
handles.allinonePush = uicontrol(handles.figure, 'style', 'push',...
        'unit','normalized',...
        'String', 'all in one', ...
        'position', allinonePushPosition,...
        'callback', @CallBack_allinonePush,...
        'visible','off');        
handles.searchPush = uicontrol(handles.figure, 'style','push',...
        'unit','normalized','position',searchPushPosition,...
        'String', 'search','callback',@CallBack_searchPush, 'ForeGroundColor',[0 0 1]);
handles.searchContentEdit = uicontrol(handles.figure, 'style','edit',...
        'unit','normalized','position',searchContentEditPosition,...
        'ForeGroundColor',[0 0 1],...
        'BackgroundColor', 'w',...
        'horizontal','left');
handles.searchText = uicontrol(handles.figure, 'style','text',...
        'unit','normalized','position',searchTextPosition,...
        'string', '  in',...
        'horizontal','left',...
        'visible','off');
handles.searchEnginePop = uicontrol(...
		'Units','normalized', ...
		'ListboxTop',0, ...
		'Position',searchEnginePopPosition, ...
		'String',{'in xBrain';'in google scholar';'in pubmed';'in wiki'}, ...
		'Style','popupmenu', ...
		'value',1);    
handles.overlayPush = uicontrol(handles.figure, 'style','push',...
        'unit','normalized','position',overlayPushPosition,...
        'String', 'overlay','callback',@CallBack_overlayPush);
handles.overlayEdit = uicontrol(handles.figure, 'style','edit',...
        'unit','normalized','position',overlayEditPosition,...
        'BackgroundColor', 'w',...
        'horizontal','left',...
        'callback', @CallBack_overlayEdit);
handles.overlayPop = uicontrol(handles.figure, 'style','popupmenu',...
        'unit','normalized','position',overlayPopPosition,...
        'string', sort(fieldnames(handles.wholeMaskMNIAll)),...
        'horizontal','left',...
        'callback', @CallBack_overlayPop);
    
handles.helpPush = uicontrol(handles.figure, 'style','push',...
        'unit','normalized','position',helpPosition,...
        'String', 'help','callback',['web ' XJVIEWURL],'ForeGroundColor',[0 0 1],...
        'horizontal','left', ...
        'visible','off');    
handles.infoTextBox = uicontrol(handles.figure, 'style','edit',...
        'unit','normalized','position',infoTextBoxPosition,...
        'String', 'Welcome to xjView 8','ForeGroundColor','k', 'BackgroundColor', 'w',...
        'horizontal','left', ...'fontname','times',...
        'max',2, 'min',0); 
handles.getCurrentPosition = uicontrol(handles.figure, 'style','push',...
        'unit','normalized','position',[0.6 0.4 0.2 0.05],...
        'String', 'get xyz','callback',@CallBack_getCurrentPosition,'ForeGroundColor',[0 0 1],...
        'horizontal','left', ...
        'visible','off');      
handles.getCurrentPosition = uicontrol(handles.figure, 'style','push',...
        'unit','normalized','position',smallVolumePosition,...
        'String', 'Small volume','callback',@CallBack_smallVolume,...
        'horizontal','left');     
    

    
try
    urlread([XJVIEWURL '/stat.php']);    
	s = urlread([XJVIEWURL '/toUser.txt']);
    set(handles.infoTextBox, 'String', s);
end

% feed
try
    jObject = com.mathworks.mlwidgets.html.HTMLBrowserPanel;
    [browser,container] = javacomponent(jObject, [], f);
    set(container, 'Units','norm', 'Pos',[0.55,0.05,0.44,0.45]);
    browser.setCurrentLocation('http://www.alivelearn.net/xjview/ad.php');
    handles.feed = container;
%     s = urlread(['http://www.alivelearn.net/xjview8/getFeed.php']);
%     [feeds, links]=strread(s,'%s%s', 'delimiter','\t');
%     
%     handles.feed = uicontrol(handles.figure,'style','listbox',...
%         'unit','normalized',...
%         'String', feeds, ...
%         'userdata',links, ...
%         'value',1,...
%         'fontsize',12,...
%         'position',[0.55 0.05 0.44 0.25]);
%     handles.refreshclick  = uicontrol(handles.figure,'style','push',...
%         'unit','normalized',...
%         'String', 'Refresh', ...
%         'position',[0.55 0.02 0.1 0.03],...
%         'fontsize',12, ...
%         'callback',@refreshFeeds);
%     handles.feedclick  = uicontrol(handles.figure,'style','push',...
%         'unit','normalized',...
%         'String', 'Read', ...
%         'position',[0.88 0.02 0.1 0.03],...
%         'fontsize',12, ...
%         'callback',@openURL);
%     handles.postfeedclick  = uicontrol(handles.figure,'style','push',...
%         'unit','normalized',...
%         'String', 'Post', ...
%         'position',[0.65 0.02 0.05 0.03],...
%         'fontsize',12, ...
%         'callback','web(''http://www.alivelearn.net/xjview8/feed.php'')');    
end

handles.glassViewAxes = axes('unit','normalized','position',glassViewAxesPosition,'XTick',[],'YTick',[],'visible','off');


handles.testEdit = uicontrol(handles.figure,'style','edit',...
        'unit','normalized','position',[0.1 0.4 0.2 0.05],...
        'horizontal','left',...
        'callback',@test,...
        'visible','off');

% menu
cSHH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')
hMenuFile = findobj(get(handles.figure,'Children'),'flat','Label','&File');
if ~isempty(hMenuFile)
	hMenuFileOpen = findobj(get(handles.figure,'Children'),'Label','&Open...');
	set(hMenuFileOpen, 'label', 'Open Figure...');
	hMenuFileSave = findobj(get(handles.figure,'Children'),'Label','&Save');
	set(hMenuFileSave, 'label', 'Save Figure ...');
	hMenuFileSaveAs = findobj(get(handles.figure,'Children'),'Label','Save &As...');
	set(hMenuFileSaveAs, 'label', 'Save Figure As ...');
else
    hMenuFile = uimenu(handles.figure, 'label', '&File');
end

set(hMenuFile,'ForegroundColor',[0 0 1]);
set(findobj(hMenuFile,'Position',1),'Separator','on');
uimenu('Parent',hMenuFile,'Position',1,'ForegroundColor',[0 0 1],...
	'Label','Open Images (*.img) ...',...
	'CallBack',@CallBack_loadImagePush, 'Accelerator', 'o');
uimenu('Parent',hMenuFile,'Position',2,'ForegroundColor',[0 0 1],...
	'Label','Save Current Image (*.img) ...',...
	'CallBack',{@CallBack_saveImagePush, '', 0});
uimenu('Parent',hMenuFile,'Position',3,'ForegroundColor',[0 0 1],...
	'Label','Save Current Image as Mask (*.img) ...',...
	'CallBack',{@CallBack_saveImagePush, '', 1});
% uimenu('Parent',hMenuFile,'Position',4,'ForegroundColor',[0 0 1],...
% 	'Label','Save Result (*.ps/pdf) ...',...
% 	'CallBack',@CallBack_saveResultPSPush);
    
hMenuHelp = findobj(get(handles.figure,'Children'),'flat','Label','&Help');
if isempty(hMenuHelp)
    hMenuHelp = uimenu(handles.figure, 'label', 'xjView &Help');
end
set(hMenuHelp,'ForegroundColor',[0 0 1]);
uimenu('Parent',hMenuHelp,'Position',1,'ForegroundColor',[0 0 1],...
	'Label','xBrain.org: brain mapping database',...
	'CallBack','web http://www.xbrain.org');
uimenu('Parent',hMenuHelp,'Position',2,...
	'Label','xjview help','ForegroundColor',[0 0 1],...
	'CallBack',['web ' XJVIEWURL]);
set(findobj(hMenuHelp,'Position',3),'Separator','on');
set(0,'ShowHiddenHandles',cSHH)

if exist('cuixuBOLDretrieve')
    hMenuAnalyze = uimenu('label','&Analyze','ForegroundColor',[0 0 1],'visible','on');
    hMenuPreprocess = uimenu(hMenuAnalyze,'label','Preprocess CIBSR','ForegroundColor',[0 0 1],'callback', @CallBack_preprocess_cibsr);    
    hMenuPreprocess = uimenu(hMenuAnalyze,'label','Preprocess','ForegroundColor',[0 0 1],'callback', @CallBack_preprocess);
    hMenuProcess = uimenu(hMenuAnalyze,'label','Process (GLM estimation)','ForegroundColor',[0 0 1],'callback',@CallBack_process);
    hMenuSPMProcess = uimenu(hMenuAnalyze,'label','SPMProcess (GLM using SPM)','ForegroundColor',[0 0 1],'callback',@CallBack_SPMProcess);
    hMenuGLMPeak = uimenu(hMenuAnalyze,'label','GLM on peak BOLD','ForegroundColor',[0 0 1],'callback',@CallBack_GLMPeak);
    hMenuContrast = uimenu(hMenuAnalyze,'label','Contrast','ForegroundColor',[0 0 1],'callback',@CallBack_contrast);
    hMenuFDR = uimenu(hMenuAnalyze,'label','FDR','ForegroundColor',[0 0 1],'callback',@CallBack_fdr);
    hMenuROI = uimenu(hMenuAnalyze,'label','ROI: retrieve signal','ForegroundColor',[0 0 1],'callback',@CallBack_timeSeries);
    hMenuROIPlot = uimenu(hMenuAnalyze,'label','ROI: plot','ForegroundColor',[0 0 1],'callback',@CallBack_plotROI, 'Accelerator', 'm');
    hMenuROIIndividualPlot = uimenu(hMenuAnalyze,'label','ROI: individual plot','ForegroundColor',[0 0 1],'callback',@CallBack_plotIndividualROI);   
    hMenuROIIndividualPlotWithBehavior = uimenu(hMenuAnalyze,'label','ROI: individual plot with behavior','ForegroundColor',[0 0 1],'callback',@CallBack_plotIndividualROIWithBehavior);
    hMenuROICorrelationPlot = uimenu(hMenuAnalyze,'label','ROI: correlation plot','ForegroundColor',[0 0 1],'callback',@CallBack_plotCorrelationROI);
    %hMenuWholeBrainCorrelation = uimenu(hMenuAnalyze,'label','Whole brain correlation','ForegroundColor',[0 0 1],'callback',@CallBack_wholeBrainCorrelation);
    hMenuBehaviorAnalysis = uimenu(hMenuAnalyze,'label','Behavior analsyis','ForegroundColor',[0 0 1],'callback',@CallBack_behaviorAnalysis);
	hMenuHeadMovementAnalysis = uimenu(hMenuAnalyze,'label','Other analysis (head motion, gsr, physio etc)','ForegroundColor',[0 0 1],'callback',@CallBack_headMovementAnalysis);
    hMenuModelComparison = uimenu(hMenuAnalyze,'label','Linear model comparison','ForegroundColor',[0 0 1],'callback',@CallBack_modelComparison);
    
    hMenuHNLOnly = uimenu('label','For H&NL Only','ForegroundColor',[0 0 1],'visible','on');
    hMenuNew2Old = uimenu(hMenuHNLOnly,'label','Format convert','ForegroundColor',[0 0 1],'callback',@CallBack_new2old);
    hMenuPreprocessCluster = uimenu(hMenuHNLOnly,'label','Preprocess (on cluster)','ForegroundColor',[0 0 1],'callback',{@CallBack_preprocess, 'cluster'});
    hMenuProcessCluster = uimenu(hMenuHNLOnly,'label','Process (GLM, on cluster)','ForegroundColor',[0 0 1],'callback',{@CallBack_process, 'cluster'});
    hMenuSPMProcessCluster = uimenu(hMenuHNLOnly,'label','SPM Process (GLM using SPM, on cluster)','ForegroundColor',[0 0 1],'callback',{@CallBack_SPMProcess, 'cluster'});
    hMenuGLMPeakCluster = uimenu(hMenuHNLOnly,'label','GLM on peak BOLD (on cluster)','ForegroundColor',[0 0 1],'callback',{@CallBack_GLMPeak, 'cluster'});
    
    hMenuNIRS = uimenu('label','&NIRS','ForegroundColor',[0 0 1],'visible','on');
    hMenuNIRS_SPM = uimenu(hMenuNIRS,'label','NIRS-SPM','ForegroundColor',[0 0 1],'callback', @CallBack_nirs_spm);        
    hMenuNIRS_TOPO = uimenu(hMenuNIRS,'label','xTOPO','ForegroundColor',[0 0 1],'callback', @CallBack_nirs_topo);        
    hMenuNIRS_nirs = uimenu(hMenuNIRS,'label','nirs','ForegroundColor',[0 0 1],'callback', @CallBack_nirs);        
    
    hMenuOther = uimenu('label','&Other','ForegroundColor',[0 0 1],'visible','on');
    hMenuOther_signal = uimenu(hMenuOther,'label','Signal processing toolbox','ForegroundColor',[0 0 1],'callback', @CallBack_signal_processing);        
    hMenuOther_ge2sn = uimenu(hMenuOther,'label','ge2sn','ForegroundColor',[0 0 1],'callback', @CallBack_ge2sn);            
end

figurecm = uicontextmenu;
uimenu(figurecm,'label','Red','callback','set(gcf,''color'',''r'')');
uimenu(figurecm,'label','White','callback','set(gcf,''color'',''w'')');
uimenu(figurecm,'label','Gray','callback','set(gcf,''color'',[0.925, 0.914,  0.847])');
%set(handles.figure,'uicontextmenu',figurecm);
set(handles.figure,'WindowButtonDownFcn',@figureMouseUpFcn);
set(handles.figure,'KeyPressFcn',@KeyPressFcn);

set(handles.figure,'visible','on');

% save pre-set values
handles.system = os;
handles.spmdir = spmdir;
handles.screenResolution = screenResolution;
handles.pValue = pValue;
handles.intensityThreshold = intensityThreshold;
handles.clusterSizeThreshold = clusterSizeThreshold;
handles.sectionViewPosition = sectionViewPosition;
handles.sectionViewTargetFile = getSectionViewTargetFile(spmdir, 'avg152T1');

guidata(f, handles);

% global variables for rotation matrix M and dimension
global M_;
global DIM_;
global TR_;
global LEFTRIGHTFLIP_;
global TMAX_; % colorbar max to display in section view
global TMIN_; % colorbar min to display in section view
global XJVIEWURL_;

M_ = M;
DIM_ = DIM;
TR_ = TR;
LEFTRIGHTFLIP_ = leftrightflip;
TMAX_ = 'auto';
TMIN_ = 'auto';
XJVIEWURL_ = XJVIEWURL;

% check input arguments
if length(varargin) == 0
    [];
elseif isstr(varargin{1})
    CallBack_loadImagePush(handles.loadImagePush, [], varargin);
else
    mniCoord = varargin{1};
    if length(varargin) < 2
        intensity = ones(size(mniCoord,1),1);
    else
        intensity = varargin{2};
    end
    thisStruct.mni = mniCoord;
    thisStruct.intensity = intensity;
    thisStruct.M = M;
    thisStruct.DIM = DIM';    
    CallBack_loadImagePush(handles.loadImagePush, [], thisStruct);
end



function test(hObject, eventdata)
handles = guidata(gcbo);
set(hObject, 'String', num2str(handles.pValue));
vars = evalin('base','who');
vars
x = evalin('base',vars{1});
x
handles.DIM

function CallBack_getCurrentPosition(hObject, eventdata)
handles = guidata(hObject);
disp(handles.currentxyz)


function KeyPressFcn(src,evt)
handles = guidata(src);
m = handles.M{1};
m = abs(diag(m));
if(strcmp(evt.Key, 'leftarrow'))
    xyz = spm_XYZreg('GetCoords',handles.hReg);
    xyz(1) = xyz(1)-m(1);
    spm_XYZreg('SetCoords', xyz, handles.hReg);
elseif(strcmp(evt.Key, 'rightarrow'))
    xyz = spm_XYZreg('GetCoords',handles.hReg);
    xyz(1) = xyz(1)+m(1);
    spm_XYZreg('SetCoords', xyz, handles.hReg);
elseif(strcmp(evt.Key, 'pageup'))
    xyz = spm_XYZreg('GetCoords',handles.hReg);
    xyz(2) = xyz(2)+m(2);
    spm_XYZreg('SetCoords', xyz, handles.hReg);    
elseif(strcmp(evt.Key, 'pagedown'))
    xyz = spm_XYZreg('GetCoords',handles.hReg);
    xyz(2) = xyz(2)-m(2);
    spm_XYZreg('SetCoords', xyz, handles.hReg);    
elseif(strcmp(evt.Key, 'uparrow'))
    xyz = spm_XYZreg('GetCoords',handles.hReg);
    xyz(3) = xyz(3)+m(3);
    spm_XYZreg('SetCoords', xyz, handles.hReg);        
elseif(strcmp(evt.Key, 'downarrow'))
    xyz = spm_XYZreg('GetCoords',handles.hReg);
    xyz(3) = xyz(3)-m(3);
    spm_XYZreg('SetCoords', xyz, handles.hReg);            
end
 

function refreshFeeds(hObject, eventdata)
handles = guidata(hObject);
s = urlread(['http://www.alivelearn.net/xjview8/getFeed.php']);
[feeds, links]=strread(s,'%s%s', 'delimiter','\t');
set(handles.feed, 'String', feeds);    
set(handles.feed, 'userdata', links);


function openURL(hObject, eventdata)
handles = guidata(hObject);
value = get(handles.feed,'value');
urls  = get(handles.feed,'userdata');
if(~isempty(urls{value}))
    web(urls{value});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NIRS-SPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_nirs_spm(hObject, eventdata)
try
    nirs_spm;
catch
    disp('Please add NIRS-SPM path.');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NIRS-TOPO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_nirs_topo(hObject, eventdata)
try
    warndlg('You need to run xTopo in command line. Type topo and press enter.');
catch
    disp('Please add NIRS-SPM path.');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NIRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_nirs(hObject, eventdata)
try
    nirs;
catch
    disp('Please add NIRS-SPM path.');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NIRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_signal_processing(hObject, eventdata)
try
    signalPreprocess;
catch
    disp('Please add NIRS-SPM path.');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NIRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_ge2sn(hObject, eventdata)
try
    xge2sn;
catch
    disp('Please add NIRS-SPM path.');
    return
end

  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FDR edit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_fdrEdit(hObject, eventdata)
handles = guidata(hObject);
set(handles.infoTextBox, 'string', {'Right now FDR only works for T-test image.'});
if length(handles.imageFileName) > 1
    set(handles.infoTextBox, 'string', {'FDR can only work for a single image file. You opened multiple files.'});
    return
end

q = get(handles.fdrEdit,'String');
q = str2num(q);
if get(handles.allIntensityRadio, 'Value')
    set(handles.infoTextBox, 'string', {'FDR only works on positive or negative direction, not both.'});
    return;
elseif get(handles.positiveIntensityRadio, 'Value')
    positive = 1;
elseif get(handles.negativeIntensityRadio, 'Value')
    positive = -1;
end

VspmSv = spm_vol(handles.imageFileName{1});
Ts = spm_read_vols(VspmSv);
Ts(Ts==0) = [];
Ts(isnan(Ts)) = [];
Ts = positive*Ts;
Ts = sort(Ts(:));
if handles.TF{1} ~= 'P', Ts = flipud(Ts); end

ps = t2p(Ts, handles.df{1}, handles.TF{1});

intensity  = spm_uc_FDR(q, [1 handles.df{1}],handles.TF{1},1,ps,0);

pvalue = t2p(intensity, handles.df{1}, handles.TF{1});
set(handles.pValueEdit,'string',pvalue);
CallBack_pValueEdit(handles.pValueEdit, eventdata);
%pause(1)
%CallBack_pValueEdit(handles.pValueEdit, eventdata);

% set(handles.intensityThresholdEdit,'string',intensity);
% pause(2)
% CallBack_intensityThresholdEdit(handles.intensityThresholdEdit,
% eventdata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FDR push button
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_fdrPush(hObject, eventdata)
handles = guidata(hObject);
set(handles.infoTextBox, 'string', {'Right now FDR only works for T-test image.'});
if length(handles.imageFileName) > 1
    set(handles.infoTextBox, 'string', {'I can only work for a single image file. You opened multiple files.'});
    return
end

q = get(handles.pValueEdit,'String');
q = str2num(q);
if get(handles.allIntensityRadio, 'Value')
    set(handles.infoTextBox, 'string', {'FDR only works on positive or negative direction, not both.'});
    return;
elseif get(handles.positiveIntensityRadio, 'Value')
    positive = 1;
elseif get(handles.negativeIntensityRadio, 'Value')
    positive = -1;
end

VspmSv = spm_vol(handles.imageFileName{1});
Ts = spm_read_vols(VspmSv);
Ts(Ts==0) = [];
Ts(isnan(Ts)) = [];
Ts = positive*Ts;
Ts = sort(Ts(:));
if handles.TF{1} ~= 'P', Ts = flipud(Ts); end

ps = t2p(Ts, handles.df{1}, handles.TF{1});

intensity  = spm_uc_FDR(q, [1 handles.df{1}],handles.TF{1},1,ps,0);

pvalue = t2p(intensity, handles.df{1}, handles.TF{1});
set(handles.pValueEdit,'string',pvalue);
%CallBack_pValueEdit(handles.pValueEdit, eventdata);
pause(1)
%CallBack_pValueEdit(handles.pValueEdit, eventdata);

% set(handles.intensityThresholdEdit,'string',intensity);
% pause(2)
% CallBack_intensityThresholdEdit(handles.intensityThresholdEdit, eventdata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FDR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_fdr(hObject, eventdata)
handles = guidata(hObject);
set(handles.infoTextBox, 'string', {'Right now FDR only works for T-test image.'});
if length(handles.imageFileName) > 1
    set(handles.infoTextBox, 'string', {'I can only work for a single image file. You opened multiple files.'});
    return
end

q = get(handles.pValueEdit,'String');
q = str2num(q);
if get(handles.allIntensityRadio, 'Value')
    positive = 1;
elseif get(handles.positiveIntensityRadio, 'Value')
    positive = 1;
elseif get(handles.negativeIntensityRadio, 'Value')
    positive = -1;
end
xjviewpath = fileparts(which('xjview'));
maskImageFile = fullfile(xjviewpath, 'mask.img');
[threshold, pvalue] = fdr(handles.imageFileName{1}, q, positive, maskImageFile);
set(handles.pValueEdit,'string',pvalue);
CallBack_pValueEdit(handles.pValueEdit, eventdata);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% model comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_modelComparison(hObject, eventdata)
answer = getVariable({'y','x (full model)', 'x (reduced model)'});

if isempty(answer)
    return;
end

if ~isempty(answer{1})
    y = evalin('base',answer{1});
else
    return;
end
if ~isempty(answer{2})
    xf = evalin('base',answer{2}); % full model
else
    return;
end
if ~isempty(answer{3})
    xr = evalin('base',answer{3}); % reduced model
else
    xr = [];
end

% make vector column vector
[r, c] = size(y);
if r==1; y = y'; end
[r, c] = size(xf);
if r<c; xf = xf'; end
[r, c] = size(xr);
if r<c; xr = xr'; end

format long;
disp('------------------------------------------')
TotalVariance = var(y)
n = size(y, 1);
dfTotal = n - 1

disp('------------------------------------------')
disp('Full Model:');
beta = linearregression(y,xf)
predictedy = [xf ones(size(xf, 1),1)] * beta;
residual = y - predictedy;
r2 = 1 - (var(residual) / var(y))
ResidualVarianceFullModel =  var(residual)
ModelVarianceFullModel =  var(predictedy);
dfResidualFull = n - size(xf, 2) - 1

if isempty(xr)
    return;
end

disp('------------------------------------------')
disp('Reduced Model:');
beta = linearregression(y,xr)
predictedy = [xr ones(size(xr, 1),1)] * beta;
residual = y - predictedy;
r2 = 1 - (var(residual) / var(y))
ResidualVarianceReducedModel =  var(residual)
ModelVarianceReducedModel =  var(predictedy);
dfResidualReduced = n - size(xr, 2) - 1

disp('------------------------------------------')
disp('Model Comparison (is the full model significantly better than the reduced model?):');
dfDenominator = dfResidualFull;
dfNumerator = dfResidualReduced - dfResidualFull;
F = (ResidualVarianceReducedModel - ResidualVarianceFullModel)/dfNumerator / (ResidualVarianceFullModel / dfDenominator);
disp(['f(' num2str(dfNumerator) ',' num2str(dfDenominator) ')=' num2str(F)])
pvalue = 1-spm_Fcdf(F, [dfNumerator dfDenominator])
% 
% F = ModelVarianceFullModel/size(xf,2) / (ModelVarianceReducedModel/size(xr,2));
% disp(['f(' num2str(size(xf,2)) ',' num2str(size(xr,2)) ')=' num2str(F)])
% pvalue = 1-spm_Fcdf(F, [size(xf,2) size(xr,2)])

format short;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% behavior analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_behaviorAnalysis(hObject, eventdata)
global TR_;
answer = getVariable({'event data file (optional)','correlator (optional)', 'head motion, gsr, physio, eye tracking (optional)'});

if isempty(answer)
    return;
end

if ~isempty(answer{1})
    Pe = evalin('base',answer{1}); % event files
    numsubj = size(Pe, 1);
else
    Pe = [];
end
if ~isempty(answer{2})
    Pc = evalin('base',answer{2}); % correlator files
    numsubj = size(Pc, 1);
else
    Pc = [];
end
if ~isempty(answer{3})
    P = evalin('base',answer{3}); % headmovment files
    numsubj = size(P, 1);
else
    P = [];
end

if ~isempty(P)
    % read what data in there (P) from the 1st dataset
    tmp = load(deblank(P(1,:)));
    tmp1 = fields(tmp);
    eval(['M = tmp.' tmp1{1} ';']);
    dataNames = fields(M);
    % prompt for TR (s)
    promt = dataNames;
    defaultanswer = repmat({num2str(TR_)}, 1, length(dataNames));
    name='Sampling interval (TR) in seconds';
    options.Resize = 'on';
    TRs=inputdlg(promt,name,1,defaultanswer, options);
    if isempty(TRs);
        return;
    end
    TRs = str2double(TRs);
end
    
colors=[1 0 0;
    0 0 1;
    0 0 0;
    1 0 1;
    0 1 1;
    0 1 0;
    1 1 0;
    1 0.5 0.5;
    .5 .5 1;
    .5 .5 .5;
    1 .5 0;
    1 0 .5;
    0 1 .5;
    .5 1 0;
    0 .5 1;
    .5 0 1];
colors = repmat(colors, 10, 1);

% allow user to select which subject to plot
prompt=['I find you have ' num2str(numsubj) ' subjects.  Please let me know which subject(s) you want me to use to generate the plots.'];
name='Which subject(s) to plot?';
numlines=1;
defaultanswer={['[1:'  num2str(numsubj) ']']};
 
whichtoplot=inputdlg(prompt,name,numlines,defaultanswer);
if isempty(whichtoplot)
    return;
end
whichtoplot = eval(whichtoplot{1});

if isempty(Pe) && isempty(Pc) && ~isempty(P)
    figure;
    Maximize(gcf);
    for ii=whichtoplot

        kk = mod(ii,8);
        if kk == 0
            kk = 8;
        end
        subplot(4,2,kk)

        tmp = load(deblank(P(ii,:)));
        tmp2 = fieldnames(tmp);
        hm = eval(['tmp.' tmp2{1}]);
        hmcell = struct2cell(hm);
        hmname = fieldnames(hm);
        for jj=1:length(hmname)
            plot(TRs(jj)*[0:length(hmcell{jj})-1], hmcell{jj}, 'color', colors(jj,:));
            hold on
        end
        legend(hmname)
        title(['subject ' num2str(ii)])
        xlabel('time (s)')

        if mod(ii,8) == 0
            figure;
            Maximize(gcf);
        end    

    end

    return
end

%% event
if ~isempty(Pe)
    figure;
    Maximize(gcf);
    for ii=whichtoplot
        kk = mod(ii,4);
        if kk == 0
            kk = 4;
        end
        subplot(4,2,2*kk-1)
        
        hmmin = 0;
        hmmax = 1;        
        if ~isempty(P) % if there is head movement data, gsr data, eye tracking data, physio data etc
            tmp = load(deblank(P(ii,:)));
            tmp2 = fieldnames(tmp);
            Ohm = eval(['tmp.' tmp2{1}]);
            Ohmcell = struct2cell(Ohm);
            Ohmname = fieldnames(Ohm);
            hmmin = inf;
            hmmax = -inf;
            for jj=1:length(Ohmname)
                mn = min(Ohmcell{jj});
                mx = max(Ohmcell{jj});
                if hmmin > mn
                    hmmin = mn;
                end
                if hmmax < mx
                    hmmax = mx;
                end                
            end        
        end

        tmp = load(deblank(Pe(ii,:)));
        tmp2 = fieldnames(tmp);
        hm = eval(['tmp.' tmp2{1}]);
        hmcell = struct2cell(hm);
        hmname = fieldnames(hm);
        for jj=1:length(hmname)
            toplot1 = [];
            toplot2 = [];
            for mm = 1:length(hmcell{jj})
                toplot1 = [toplot1 hmcell{jj}(mm) hmcell{jj}(mm) nan];
                toplot2 = [toplot2 hmmin hmmax nan];
            end
            plot(toplot1, toplot2, 'color', colors(jj,:));
            %ylim([0 4])
            hold on
        end
        legend(hmname)
        title(['subject ' num2str(ii)])
        xlabel('time in s');
        
        if ~isempty(P) % if there is head movement data, gsr data, eye tracking data, physio data etc
            for jj=1:length(Ohmname)
                plot(TRs(jj)*[0:length(Ohmcell{jj})-1],Ohmcell{jj},'color', colors(jj,:));
                hold on
            end        
            legend([hmname; Ohmname])
        end
        
        %% prepare interval data for later plot
        alltimes{ii} = [];
        for jj=1:length(hmname)
            alltimes{ii} = union(alltimes{ii}, hmcell{jj});
        end

        subplot(4,2,2*kk)
        for jj=1:length(hmname)
            roundhmcell{jj} = round(hmcell{jj}*1);
            signal = zeros(1, max(roundhmcell{jj})+1);
            signal(roundhmcell{jj}+1) = 1;
            signal = signal - mean(signal);
            y = fft(signal);
            plotlength = round(length(y))/1;
            plot([1:plotlength]/plotlength, abs(y(1:plotlength)).^2, 'color', colors(jj,:));
            hold on
        end
        legend(hmname)
        xlabel('frequency');
        ylabel('fft power');
        title(['subject ' num2str(ii)])

        if mod(ii,4) == 0 && ii~=size(Pe, 1)
            figure;
            Maximize(gcf);
        end     
    end

    %% now I plot correlogram of intervals
    figure;
    kk = 1;
    for ii=whichtoplot
        alltimes{ii} = sort(alltimes{ii});
        intervals = diff(alltimes{ii});
        corrlag(intervals,intervals,[max(-10, -length(intervals)+2) : min(10, length(intervals)-2)], colors(ii,:));
        hold on;
        legendname{kk} = ['sub ' num2str(ii)];
        kk = kk+1;
    end
    legend(legendname)
    xlabel('lag');
    ylabel('correlation coefficient');
    title('correlogram of intervals')
end

%% correlator

if ~isempty(Pc)
    
	tmp = load(deblank(Pc(1,:)));
	tmp2 = fieldnames(tmp);
	hm = eval(['tmp.' tmp2{1}]);

    C = hm;
    cenames = fields(C);
    for jj=1:length(cenames)
        eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
        for kk=1:length(correlatornames{jj})
            eval(['C.' cenames{jj} '.' correlatornames{jj}{kk} '= [];']);
        end
    end
    
    %list all plot names in command window for selection
    promt{1} = 'Please select which two correlators do you want to plot. The first one will be x and the second will be y.  The two have to be of the same length.';
    promt{2} = 'index   event    correlator';
    count = 1;
    for jj=1:length(cenames)
        eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
        for kk=1:length(correlatornames{jj})
            promt{count+2} = sprintf('%d         %s         %s', count,  cenames{jj}, correlatornames{jj}{kk});
            count = count + 1;
        end
    end
    promt{count+2} = '';

    if count == 1
        return
    end
    promt = char(promt);
    name='Please select which two correlator do you want to plot. The first one will be x and the second will be y.';
    numlines=1;
    if count == 1
        defaultanswer={['You don''t have more than one correlators. Click Cancel']};
    else
        defaultanswer={['[1 '  num2str(count-1) ']']};
    end

    whichplot=inputdlg(promt,name,numlines,defaultanswer);
    if isempty(whichplot)
        return;
    end    
    whichplot = eval(whichplot{1});
    
    figure;
    Maximize(gcf);
    for ii=whichtoplot
        tmp = load(deblank(Pc(ii,:)));
    	tmp2 = fieldnames(tmp);
        hm = eval(['tmp.' tmp2{1}]);
        count = 1;
        for jj=1:length(cenames)
            eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
            for kk=1:length(correlatornames{jj})
                tmp2 = eval(['hm.' cenames{jj} '.' correlatornames{jj}{kk} ';']);
                if size(tmp2,1) == 1
                    tmp2 = tmp2';
                end
                eval(['C.' cenames{jj} '.' correlatornames{jj}{kk} '= [C.' cenames{jj} '.' correlatornames{jj}{kk} '; tmp2];']);
                
                if whichplot(1) == count
                    x = tmp2;
                    xlabeltext = [correlatornames{jj}{kk} ' at ' cenames{jj}];
                end
                if whichplot(2) == count
                    y = tmp2;
                    ylabeltext = [correlatornames{jj}{kk} ' at ' cenames{jj}];
                end
                count = count + 1;
            end
        end
        
        kk = mod(ii,16);
        if kk == 0
            kk = 16;
        end
        subplot(4,4,kk)
        
        plot2(x,y, 'b', 1);
        xlabel(xlabeltext)
        ylabel(ylabeltext)
        title(['subject ' num2str(ii)])    
        try
            % convert to colum vector
            if size(x,1)==1
                x = x';
            end
            if size(y,1)==1
                y = y';
            end

            % remove NaN
            tmpposx = find(isnan(x) | isinf(x));
            tmpposy = find(isnan(y) | isinf(y));
            x([tmpposx; tmpposy]) = [];
            y([tmpposx; tmpposy]) = [];

            b = linearregression(y,x);

            totalvar = var(y);
            residvar = var(y - x*b(1) - b(2));
            modelvar = totalvar - residvar;
            dfm = 1;
            dft = length(x) - 1;
            dfr = dft - dfm;
            mmodelvar = modelvar/dfm;
            mresidvar = residvar/dfr;
            F = mmodelvar/mresidvar;
            pvalue = 1-spm_Fcdf(F, [dfm dfr]);

            title(sprintf('sub %d, pValue=%g', ii, pvalue))
            xx=[min(x):(max(x)-min(x))/10:max(x)];
            yy=b(1) * xx + b(2);
            hold on;
            plot(xx,yy,'g');
        end
        
        if mod(ii,16) == 0 && ii~=size(Pc, 1)
            figure;
            Maximize(gcf);
        end    
    end

    % plot the average
    count = 1;
    for jj=1:length(cenames)
        eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
        for kk=1:length(correlatornames{jj})
            if whichplot(1) == count
                eval(['x = C.' cenames{jj} '.' correlatornames{jj}{kk} ';']);
                xlabeltext = [correlatornames{jj}{kk} ' at ' cenames{jj}];
            end
            if whichplot(2) == count
                eval(['y = C.' cenames{jj} '.' correlatornames{jj}{kk} ';']);
                ylabeltext = [correlatornames{jj}{kk} ' at ' cenames{jj}];
            end
            count = count + 1;
        end
    end
    
    figure
    plot2(x,y, 'b', 1);
    xlabel(xlabeltext)
    ylabel(ylabeltext)
    try
        % convert to colum vector
        if size(x,1)==1
            x = x';
        end
        if size(y,1)==1
            y = y';
        end

        % remove NaN
        tmpposx = find(isnan(x) | isinf(x));
        tmpposy = find(isnan(y) | isinf(y));
        x([tmpposx; tmpposy]) = [];
        y([tmpposx; tmpposy]) = [];

        b = linearregression(y,x);

        totalvar = var(y);
        residvar = var(y - x*b(1) - b(2));
        modelvar = totalvar - residvar;
        dfm = 1;
        dft = length(x) - 1;
        dfr = dft - dfm;
        mmodelvar = modelvar/dfm;
        mresidvar = residvar/dfr;
        F = mmodelvar/mresidvar;
        pvalue = 1-spm_Fcdf(F, [dfm dfr]);

        title(sprintf('pValue=%g; slope=%g; constant=%g', pvalue, b(1), b(2)))
        xx=[min(x):(max(x)-min(x))/10:max(x)];
        yy=b(1) * xx + b(2);
        hold on;
        plot(xx,yy,'g');
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% head motion, physio, gsr, eyetracking etc analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_headMovementAnalysis(hObject, eventdata)

global TR_;
handles = guidata(hObject);

answer = getVariable({'head motion, gsr, physio, eye tracking data', 'Event time data', 'Correlator (optional)'});

if isempty(answer)
    return;
end
if isempty(answer{1})
    return;
end

[filename, pathname] = uiputfile('*.mat', 'Save result to', '');
if isequal(filename,0) | isequal(pathname,0)
   return
else
   thisfilename = fullfile(pathname, filename);
end


P = evalin('base',answer{1}); % data
try
    E = evalin('base',answer{2}); % event
catch
    E = [];
end
try
    C = evalin('base',answer{3}); % correlator
catch
    C = [];
end

% read what data in there from the 1st dataset
tmp = load(deblank(P(1,:)));
tmp1 = fields(tmp);
eval(['M = tmp.' tmp1{1} ';']);
dataNames = fields(M);

% prompt for TR (s)
promt = dataNames;
defaultanswer = repmat({num2str(TR_)}, 1, length(dataNames));
name='Sampling interval (TR) in seconds';
options.Resize = 'on';
TRs=inputdlg(promt,name,1,defaultanswer, options);
if isempty(TRs);
    return;
end
TRs = str2double(TRs);

% % prompt for Moving Average Detrend window size (s)
% promt = dataNames;
% defaultanswer = repmat({num2str(40)}, 1, length(dataNames));
% name='Moving average detrending window size in seconds';
% options.Resize = 'on';
% MWs=inputdlg(promt,name,1,defaultanswer, options);
% if isempty(MWs);
%     return;
% end
% MWs = str2double(MWs);
% MWs = MWs./TRs;

h = waitbar(0,'Please wait...');
    
for nn=1:length(dataNames)

    for ii=1:size(P,1)
        tmp = load(deblank(P(ii,:)));
        tmp1 = fields(tmp);
        eval(['Ms = tmp.' tmp1{1} ';']);
        eval(['M = Ms.' dataNames{nn} ';']);        
        if isnumeric(M)
            M = double(M);
        end
        if ~isempty(E)
            tmp = load(deblank(E(ii,:)));
            tmp1 = fields(tmp);
            eval(['event{ii} = tmp.' tmp1{1} ';']);        
            if ~isstruct(event{ii})
                error('I need you specify the event in a structure, not a cell array!');
                return
            end
        else
            event = {};
        end
     
        [eventResponse{ii}, time{ii}, wholeOriginal{ii}] = cuixuSignalRetrieve(M, event{ii}, 50, -20, 0, TRs(nn));

        if ~isempty(C)
            tmp = load(deblank(C(ii,:)));
            tmp1 = fields(tmp);
            eval(['correlator{ii} = tmp.' tmp1{1} ';']);        
            if ~isstruct(correlator{ii})
                error('I need you specify the correlator in a structure, not a cell array!');
                return
            end

            names = fields(correlator{ii});
            for kk=1:length(names)
                eval(['names2 = fields(correlator{ii}.' names{kk} ');']);
                for mm=1:length(names2)
                    eval(['correlator{ii}.' names{kk} '.' names2{mm} ' = correlator{ii}.' names{kk} '.' names2{mm} '(remainedEvent{ii}.' names{kk} ');']);;
                end
            end
        end
        waitbar(ii/size(P,1),h,[dataNames{nn} ': finished ' num2str(ii) ' of ' num2str(size(P,1))]);
        disp(['finished ' num2str(ii) ' of ' num2str(size(P,1))]);
    end
        

    if ~exist('correlator')
        correlator = {};
    end

    TR = TRs(nn);
    [a,b,c] = fileparts(thisfilename);
    save([fullfile(a,b) '_' dataNames{nn}], 'eventResponse', 'time', 'wholeOriginal', 'correlator', 'event', 'P', 'TR');
end
close(h)

%CallBack_plotROI(hObject, eventdata, thisfilename);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% contrast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_contrast(hObject, eventdata)

answer = getVariable({'beta files', 'contrast vector'});

if isempty(answer)
    return;
end
if isempty(answer{1}) | isempty(answer{2}) 
    return;
end

[filename, pathname] = uiputfile('*.img', 'Save t-test image as', '');
if isequal(filename,0) | isequal(pathname,0)
   return
else
   thisfilename = fullfile(pathname, filename);
end


P1 = evalin('base',answer{1}); % beta image files
c = evalin('base',answer{2}); % contrast vector

if iscell(P1)
    [];
else
    disp('format wrong');
    return;
end


h = waitbar(0,'Please wait...');
V=spm_vol(deblank(P1{1}(1,:)));
M1={spm_read_vols(V)};
contr = zeros(size(M1{1},1), size(M1{1},2), size(M1{1},3), size(P1{1},1));

for ii=1:size(P1{1},1)
    for jj=1:length(P1)
        if c(jj)~=0
            V=spm_vol(deblank(P1{jj}(ii,:)));
            M1{jj}=spm_read_vols(V);
            contr(:,:,:,ii) = contr(:,:,:,ii)+c(jj)*double(M1{jj});
        end
    end

    waitbar(ii/size(P1{1},1),h,['finished ' num2str(ii) ' of ' num2str(size(P1{1},1))]);
end
close(h);


df = size(P1{1},1)-1;

T = squeeze(t(permute(contr,[4 1 2 3])));   % return a 41*48*35 matrix
xjviewpath = fileparts(which('xjview'));
V=spm_vol(fullfile(xjviewpath, 'mask.img'));
m=spm_read_vols(V);
T = T.*m;

targetfilename = thisfilename;
global M_;
global DIM_;
V.mat = M_;
V.dim = [DIM_(1) DIM_(2) DIM_(3) 16];
V.fname = targetfilename;
V.descrip = ['SPM{T_[' num2str(df) '.0]} ' num2str(c)];
V = spm_write_vol(V,T);

handles = guidata(hObject);
[tmp,tmp,ext] = fileparts(thisfilename);
if isempty(ext)
    ext = '.img';
elseif isequal(lower(ext), '.img')
    ext = '';
end
CallBack_loadImagePush(handles.loadImagePush, [], {[thisfilename ext]});

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPM process (use SPM processing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_SPMProcess(hObject, eventdata, singleOrCluster)

answer = getVariable({'Image files', 'Event time data', 'Modulator (optional)','Other (head motion, gsr, physio, eye track, etc, optional)'});

if isempty(answer)
    return;
end
if isempty(answer{1}) | isempty(answer{2})
    return;
end

P = evalin('base',answer{1}); % subject directories
E = evalin('base',answer{2}); % event time
try
    Modulator = evalin('base',answer{3}); % modulator
catch
    Modulator = [];
end
try
    headmovement = evalin('base',answer{4}); % other regressor such as headmovement
catch
    headmovement = [];
end

if nargin < 3
    singleOrCluster = 'single';
end


prompt='I will create a folder under each subject''s directory and put the SPM output files (including beta images, SPM.mat etc) for that subject into that folder. Please give a name to the folder (example: FixedEffect). No space or strange characters allowed.';
name='Give a name to SPM output folder';
numlines=1;
defaultanswer={'FixedEffect'};
 
wheretosave=inputdlg(prompt,name,numlines,defaultanswer);

if isempty(wheretosave)
    return
else
    wheretosave = wheretosave{1};
end

currentdir = pwd;

if strcmp(singleOrCluster, 'cluster')    
    handles = guidata(hObject);
    set(handles.infoTextBox, 'string', {'This will do SPM GLM estimation using cluster.', 'Check out xjviewtmp directory'});
    
    mkdir('xjviewtmp');
    cd('xjviewtmp');
    system('rm *');

    fid = fopen('processALL.pbs','wt');
    fprintf(fid,'#!/bin/bash\n\n');
    fprintf(fid,'for ((s=1;s<=%d;s++))\n', size(P,1));
    fprintf(fid,'do\n\tsleep 3.5\n\tqsub %s/process.pbs -v "s=$s" -N SPM_$s\ndone\n', pwd);
    fclose(fid);

    save P P;
    save E E;
    save Modulator Modulator;
    save headmovement headmovement;
    fid = fopen('process.pbs','wt');
    fprintf(fid,'#!/bin/bash\n\n');
    fprintf(fid,'#PBS -l nodes=1:ppn=1,walltime=04:00:00\n');
    fprintf(fid,'#PBS -N process\n');
    fprintf(fid,'#PBS -e %s\n', pwd);
    fprintf(fid,'#PBS -o %s\n\n', pwd);
    fprintf(fid,'sleep 3.5\n\n');
    fprintf(fid,'PATH=/usr/local/bin:$PATH\n\n');
    fprintf(fid,'matlab -nojvm -nodisplay \\\n');
    fprintf(fid,'-logfile "%s/process.$s.$PBS_JOBID.log" \\\n', pwd);


    fprintf(fid,['-r "addpath ~/gang/xjview;'...
        'wheretosave = ''%s'';'...
        'cd %s; '...
        'load P; '...
        'load E; '...
        'load Modulator;'... 
        'load headmovement;'...
        'tmp=load(deblank(P($s,:)));'...
        'fd=fields(tmp);'...
        'eval([''Pi=tmp.'' fd{1} '';'']);'...
        'N = size(Pi,1);'...
        'tmp=load(deblank(E($s,:))); '...
        'fd=fields(tmp); '...
        'eval([''EE=tmp.'' fd{1} '';'']); '...
        'thisModulator = [];'...
        'if ~isempty(Modulator);'...
            'tmp = load(deblank(Modulator($s,:)));  '...
            'fd=fields(tmp); '...
            'eval([''thisModulator=tmp.'' fd{1} '';'']);'...
        'end;'...
        '[otherregressor{1}.C, otherregressor{1}.name]=struct2matrix(N, EE, thisModulator);'...
        'if ~isempty(headmovement);'...
            'tmp = load(deblank(headmovement($s,:)));  '...
            'fd=fields(tmp); '...
            'eval([''thisheadmovement=tmp.'' fd{1} '';'']);'...
            'tmphd = struct2cell(thisheadmovement);'...
            'tmpname = fieldnames(thisheadmovement);'...
            'for kk=1:length(tmphd);'...
                '[r,c] = size(tmphd{kk});'...
                'if r==1; tmphd{kk} = tmphd{kk}''; end;'...
                'otherregressor{1}.C = [otherregressor{1}.C tmphd{kk}];'...
                'otherregressor{1}.name = [otherregressor{1}.name {tmpname{kk}}];'...
            'end;'...
        'end;'...
        'subDir = fileparts(deblank(Pi(1,:))); '...
        'for kk=1:length(subDir); '...
            'if subDir(end-kk+1)==filesep;'...
                'pos=kk;'...
                'break;'...
            'end;'...
        'end; '...
        'subDir(end-pos+1:end)=[]; '...
        'cd(subDir); '...
        'try;' ...
            'rmdir(wheretosave, ''s'');'...
        'end;'...
        'mkdir(wheretosave);'...
        'cd(wheretosave);'...
        'spm_defaults;'...
        'disp(subDir);'...
        'cuixuprocess(N, [0], {}, {}, {}, [], otherregressor, Pi, 16);'...
        'disp(''done!'');'...
        'exit;"'], wheretosave, pwd);
    
    fclose(fid);

    system(['ssh cluster.hnl.bcm.tmc.edu  PATH=$PATH:/usr/local/pbs/i686/bin bash ' pwd '/processALL.pbs']);
    cd ..
elseif  strcmp(singleOrCluster, 'single')
    handles = guidata(hObject);
    set(handles.infoTextBox, 'string', {'This will do SPM GLM estimation.'});
    for ii=1:size(P,1)
        tmp=load(deblank(P(ii,:)));
        fd=fields(tmp);
        eval(['Pi=tmp.' fd{1} ';']);
        N = size(Pi,1);
        tmp=load(deblank(E(ii,:))); 
        fd=fields(tmp); 
        eval(['EE=tmp.' fd{1} ';']); 
        thisModulator = [];
        if ~isempty(Modulator)
            tmp = load(deblank(Modulator(ii,:)));  
            fd=fields(tmp); 
            eval(['thisModulator=tmp.' fd{1} ';']);
        end
        [otherregressor{1}.C, otherregressor{1}.name]=struct2matrix(N, EE, thisModulator);
        if ~isempty(headmovement)
            tmp = load(deblank(headmovement(ii,:)));  
            fd=fields(tmp); 
            eval(['thisheadmovement=tmp.' fd{1} ';']);
            
            tmphd = struct2cell(thisheadmovement);
            tmpname = fieldnames(thisheadmovement);
            for kk=1:length(tmphd)
                [r,c] = size(tmphd{kk});
                if r==1; tmphd{kk} = tmphd{kk}'; end
                otherregressor{1}.C = [otherregressor{1}.C tmphd{kk}];
                otherregressor{1}.name = [otherregressor{1}.name {tmpname{kk}}];
            end
        end
        subDir = fileparts(deblank(Pi(1,:))); 
        for kk=1:length(subDir); 
            if subDir(end-kk+1)==filesep;
                pos=kk;
                break;
            end;
        end; 
        subDir(end-pos+1:end)=[]; 
        cd(subDir); 
        try
            rmdir(wheretosave, 's');
        end
        mkdir(wheretosave);
        cd(wheretosave);

        spm_defaults;
        disp(subDir);
        cuixuprocess(N, [0], {}, {}, {}, [], otherregressor, Pi, 16);
        
        tmp = get(handles.infoTextBox, 'string');
        set(handles.infoTextBox, 'string', [tmp; {['SPM Processing subject ' num2str(ii) ' / ' num2str(size(P,1)) ' ...']}]); 
        disp(['----finished ' num2str(ii) ' out of ' num2str(size(P,1)) ' subjects----']);
    end
	tmp = get(handles.infoTextBox, 'string');
	set(handles.infoTextBox, 'string', [tmp; {'Done!'}]); 
end

cd(currentdir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GLMPeak (only do regression on peak BOLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_GLMPeak(hObject, eventdata, singleOrCluster)

answer = getVariable({'Image data', 'Event time data', 'Modulator'});

if isempty(answer)
    return;
end
if isempty(answer{1}) | isempty(answer{2}) | isempty(answer{3})
    return;
end

directoryname = uigetdir(pwd, 'Pick a Directory to Save Your Beta Files...');
if isequal(directoryname, 0)
    return;
end

peakpoint = inputdlg({'What points in time do you considered as peak BOLD?'}, 'Peak point',1,{'[4 6]'});
if isempty(peakpoint); return; end
peakpointstr = peakpoint{1};
peakpoint = str2num(peakpoint{1});
    
    
P = evalin('base',answer{1}); % subject directories
E = evalin('base',answer{2}); % event time
Modulator = evalin('base',answer{3}); % modulator


if nargin < 3
    singleOrCluster = 'single';
end

if strcmp(singleOrCluster, 'cluster')
    mkdir('xjviewtmp');
    cd('xjviewtmp');
    system('rm *');

    fid = fopen('processALL.pbs','wt');
    fprintf(fid,'#!/bin/bash\n\n');
    fprintf(fid,'for ((s=1;s<=%d;s++))\n', size(P,1));
    fprintf(fid,'do\n\tsleep 3.5\n\tqsub %s/process.pbs -v "s=$s" -N pro_$s\ndone\n', pwd);
    fclose(fid);

    save P P;
    save E E;
    save Modulator Modulator;
    fid = fopen('process.pbs','wt');
    fprintf(fid,'#!/bin/bash\n\n');
    fprintf(fid,'#PBS -l nodes=1:ppn=1,walltime=08:00:00\n');
    fprintf(fid,'#PBS -N process\n');
    fprintf(fid,'#PBS -e %s\n', pwd);
    fprintf(fid,'#PBS -o %s\n\n', pwd);
    fprintf(fid,'sleep 3.5\n\n');
    fprintf(fid,'PATH=/usr/local/bin:$PATH\n\n');
    fprintf(fid,'matlab -nojvm -nodisplay \\\n');
    fprintf(fid,'-logfile "%s/process.$s.$PBS_JOBID.log" \\\n', pwd);

    fprintf(fid,['-r "addpath ~/gang/xjview; '...
            'cd %s; '...
            'load P; '...
            'load E; '...
            'load Modulator; '...
            'tmp=load(deblank(P($s,:))); '...
            'fd=fields(tmp); '...
            'eval([''SS=tmp.'' fd{1} '';'']); '...
            'tmp=load(deblank(E($s,:))); '...
            'fd=fields(tmp); '...
            'eval([''EE=tmp.'' fd{1} '';'']); '...
            'tmp = load(deblank(Modulator($s,:)));  fd=fields(tmp);  eval([''MM=tmp.'' fd{1} '';'']);'...
                    'names = fieldnames(EE);'...
            'for jj=1:length(names);'...
                'eval([''MM.'' names{jj} ''.time = EE.'' names{jj} '';'']);'...
            'end;'...
            'xjviewpath = fileparts(which(''xjview''));'...
            'cuixuGLMpeak(MM, SS, ''%s'', [''subject'' num2str($s)], fullfile(xjviewpath, ''mask.img''), %s, 100, 2, fullfile(xjviewpath, ''templateFile.img''));'...
            'exit;"'], pwd, directoryname, peakpointstr);

    fclose(fid);

    system(['ssh cluster.hnl.bcm.tmc.edu  PATH=$PATH:/usr/local/pbs/i686/bin bash ' pwd '/processALL.pbs']);
    cd ..
elseif strcmp(singleOrCluster, 'single')
    xjviewpath = fileparts(which('xjview'));
    handles = guidata(hObject);
    set(handles.infoTextBox, 'string', {'This will do simple GLM estimation on peak BOLD.'});
    for ii = 1:size(P,1)
        tmp = get(handles.infoTextBox, 'string');
        set(handles.infoTextBox, 'string', [tmp; {['Processing subject ' num2str(ii) ' / ' num2str(size(P,1)) ' ...']}]); 
        tmp=load(deblank(P(ii,:)));
        fd=fields(tmp);
        eval(['SS=tmp.' fd{1} ';']);
        tmp=load(deblank(E(ii,:)));
        fd=fields(tmp);
        eval(['EE=tmp.' fd{1} ';']);
        tmp=load(deblank(Modulator(ii,:)));
        fd=fields(tmp);
        eval(['MM=tmp.' fd{1} ';']);
        % convert Modulator and E to be recognized by cuixuGLMpeak
        names = fieldnames(EE);
        for jj=1:length(names)
            eval(['MM.' names{jj} '.time = EE.' names{jj} ';']);
        end
        cuixuGLMpeak(MM, SS, directoryname, ['subject' num2str(ii)], fullfile(xjviewpath, 'mask.img'), peakpoint, 100, 2, fullfile(xjviewpath, 'templateFile.img'));
        disp(['----finished ' num2str(ii) ' out of ' num2str(size(P,1)) ' subjects----']);
    end
    tmp = get(handles.infoTextBox, 'string');
    set(handles.infoTextBox, 'string', [tmp; {'Done!'}]); 
end
  
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% process (use my own process code, cuixuAGLM, estimation on detrended relative BOLD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_process(hObject, eventdata, singleOrCluster)

answer = getVariable({'Image data', 'Event time data', 'Modulator'});

if isempty(answer)
    return;
end
if isempty(answer{1}) | isempty(answer{2})
    return;
end

directoryname = uigetdir(pwd, 'Pick a Directory to Save Your Beta Files...');
if isequal(directoryname, 0)
    return;
end

P = evalin('base',answer{1}); % subject directories
E = evalin('base',answer{2}); % event time
try
    Modulator = evalin('base',answer{3}); % modulator
catch
    Modulator = [];
end


if nargin < 3
    singleOrCluster = 'single';
end

if strcmp(singleOrCluster, 'cluster')
    mkdir('xjviewtmp');
    cd('xjviewtmp');
    system('rm *');

    fid = fopen('processALL.pbs','wt');
    fprintf(fid,'#!/bin/bash\n\n');
    fprintf(fid,'for ((s=1;s<=%d;s++))\n', size(P,1));
    fprintf(fid,'do\n\tsleep 3.5\n\tqsub %s/process.pbs -v "s=$s" -N pro_$s\ndone\n', pwd);
    fclose(fid);

    save P P;
    save E E;
    save Modulator Modulator;
    fid = fopen('process.pbs','wt');
    fprintf(fid,'#!/bin/bash\n\n');
    fprintf(fid,'#PBS -l nodes=1:ppn=1,walltime=08:00:00\n');
    fprintf(fid,'#PBS -N process\n');
    fprintf(fid,'#PBS -e %s\n', pwd);
    fprintf(fid,'#PBS -o %s\n\n', pwd);
    fprintf(fid,'sleep 3.5\n\n');
    fprintf(fid,'PATH=/usr/local/bin:$PATH\n\n');
    fprintf(fid,'matlab -nojvm -nodisplay \\\n');
    fprintf(fid,'-logfile "%s/process.$s.$PBS_JOBID.log" \\\n', pwd);

    fprintf(fid,['-r "addpath ~/gang/xjview; '...
            'cd %s; '...
            'load P; '...
            'load E; '...
            'load Modulator; '...
            'tmp=load(deblank(P($s,:))); '...
            'fd=fields(tmp); '...
            'eval([''M=tmp.'' fd{1} '';'']); '...
            'if isnumeric(M); M=double(M); elseif isstr(M); V=spm_vol(M); M =spm_read_vols(V); else; error(''I don not understand the fileformat''); end;' ...
            'tmp=load(deblank(E($s,:))); '...
            'fd=fields(tmp); '...
            'eval([''EE=tmp.'' fd{1} '';'']); '...
            'if ~isempty(Modulator); tmp = load(deblank(Modulator($s,:)));  fd=fields(tmp);  eval([''thisModulator=tmp.'' fd{1} '';'']); else; thisModulator=[]; end;'...
            'M = permute(M,[4 1 2 3]); '...
            '[tmp,tmp,M]=mAveDetrend(M,100); '...
            'M = permute(M,[2 3 4 1]); '...
            '[a,b,c]=fileparts(deblank(P($s,:)));'...
            'cuixuAGLM(M,EE,{''%s'', b}, thisModulator); '...
            'exit;"'], pwd, directoryname);

    fclose(fid);

    system(['ssh cluster.hnl.bcm.tmc.edu  PATH=$PATH:/usr/local/pbs/i686/bin bash ' pwd '/processALL.pbs']);
    cd ..
elseif strcmp(singleOrCluster, 'single')
    handles = guidata(hObject);
    set(handles.infoTextBox, 'string', {'This will do simple GLM estimation.'});
    for ii=1:size(P,1)
        tmp = get(handles.infoTextBox, 'string');
        set(handles.infoTextBox, 'string', [tmp; {['Processing subject ' num2str(ii) ' / ' num2str(size(P,1)) ' ...']}]); 
        tmp=load(deblank(P(ii,:)));
        fd=fields(tmp);
        eval(['M=tmp.' fd{1} ';']);
        if isnumeric(M); M=double(M); elseif isstr(M); V=spm_vol(M); M =spm_read_vols(V); else; error('I don''t understand the fileformat'); end;
        tmp=load(deblank(E(ii,:)));
        fd=fields(tmp);
        eval(['EE=tmp.' fd{1} ';']);
        if ~isempty(Modulator); tmp = load(deblank(Modulator(ii,:)));  fd=fields(tmp);  eval(['thisModulator=tmp.' fd{1} ';']); else; thisModulator=[]; end;
        M = permute(M,[4 1 2 3]);
        [tmp,tmp,M]=mAveDetrend(M,100);
        M = permute(M,[2 3 4 1]);
        [a,b,c]=fileparts(deblank(P(ii,:)));
        cuixuAGLM(M,EE,{directoryname, b}, thisModulator);
        disp(['----finished ' num2str(ii) ' out of ' num2str(size(P,1)) ' subjects----']);
    end
	tmp = get(handles.infoTextBox, 'string');
	set(handles.infoTextBox, 'string', [tmp; {'Done!'}]); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% convert new data format to old
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_new2old(hObject, eventdata)
answer = getVariable({'subjects directories'});
if isempty(answer)
    return;
end
if isempty(answer{1})
    return;
end
P = evalin('base',answer{1}); % subject directories
new2old(P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pre-process_cibsr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_preprocess_cibsr(hObject, eventdata)
answer = getVariable({'subjects directories'});
if isempty(answer)
    return;
end
if isempty(answer{1})
    return;
end
P = evalin('base',answer{1}); % subject directories

warndlg('You need to addpath ml7spm2devel2')

cuixuPreprocessCIBSR(P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_preprocess(hObject, eventdata, singleOrCluster)
answer = getVariable({'subjects directories'});
if isempty(answer)
    return;
end
if isempty(answer{1})
    return;
end
P = evalin('base',answer{1}); % subject directories

if nargin < 3
    singleOrCluster = 'single';
end

if strcmp(singleOrCluster, 'cluster')
    mkdir('xjviewtmp');
    cd('xjviewtmp');
    system('rm *');

    fid = fopen('preprocessALL.pbs','wt');
    fprintf(fid,'#!/bin/bash\n\n');
    fprintf(fid,'for ((s=1;s<=%d;s++))\n', size(P,1));
    fprintf(fid,'do\n\tsleep 3.5\n\tqsub %s/preprocess.pbs -v "s=$s" -N pre_$s\ndone\n', pwd);
    fclose(fid);

    save P P;
    fid = fopen('preprocess.pbs','wt');
    fprintf(fid,'#!/bin/bash\n\n');
    fprintf(fid,'#PBS -l nodes=1:ppn=1,walltime=08:00:00\n');
    fprintf(fid,'#PBS -N preprocess\n');
    fprintf(fid,'#PBS -e %s\n', pwd);
    fprintf(fid,'#PBS -o %s\n\n', pwd);
    fprintf(fid,'sleep 3.5\n\n');
    fprintf(fid,'PATH=/usr/local/bin:$PATH\n\n');
    fprintf(fid,'matlab -nojvm -nodisplay \\\n');
    fprintf(fid,'-logfile "%s/preprocess.$s.$PBS_JOBID.log" \\\n', pwd);
    fprintf(fid,['-r "addpath ~/gang/xjview; cd ' pwd '; load P; cuixuSmartPreprocess(P, $s);exit;"']);
    fclose(fid);

    system(['ssh cluster.hnl.bcm.tmc.edu  PATH=$PATH:/usr/local/pbs/i686/bin bash ' pwd '/preprocessALL.pbs']);
    cd ..
elseif strcmp(singleOrCluster, 'single')
    handles = guidata(hObject);
    set(handles.infoTextBox, 'string', {'What will be done in preprocessing is:',...
        '1. Converting DICOM files',...
        '2. Re-align',...
        '3. Coregister',...
        '4. Slice timing',...
        '5. Normalize',...
        '6. Smooth','',...
        'This function creates a directory ''preprocessing'' under your homeDir and put all preprocessing files there. It also saves the head movement matrix in ''headmovement.mat''. The intermediate files are removed',''});
    for ii=1:size(P,1)
        tmp = get(handles.infoTextBox, 'string');
        set(handles.infoTextBox, 'string', [tmp; {['Preprocessing subject ' num2str(ii) ' / ' num2str(size(P,1)) ' ...']}]); 
        cuixuSmartPreprocess(P, ii);
    end
	tmp = get(handles.infoTextBox, 'string');
	set(handles.infoTextBox, 'string', [tmp; {'Done!'}]); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ROI: TimeSeries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_timeSeries(hObject, eventdata)

global TR_;
handles = guidata(hObject);
if ~isfield(handles,'currentDisplayMNI') | isempty(cell2mat(handles.currentDisplayMNI'))
    %errordlg('No voxels found.','error');
	set(handles.infoTextBox, 'string', 'No voxels found.'); 
    beep
    return;
end

dimensionmatch = questdlg('The dimension of your mask image should be identical to the dimensions of the images where the time series are extracted.  Are they identical?', ...
                         'Dimension Match Question', ...
                         'Yes', 'No', 'I Do Not Know', 'Yes');

switch dimensionmatch
    case 'No',
        msgbox('You should use spm reslice to change your mask image dimensions','warning');
        return;
    case 'I Do Not Know',
        msgbox('You should check the dimensions and use spm reslice to change your mask image dimensions','warning');
        return;
end % switch
   


answer = getVariable({'Image data', 'Event time data', 'Correlator', 'Effects to be removed (e.g. headmovement)'});
%answer=inputdlg('Input the name of the variable (in the base workspace) which contains the image file name list, or the 4-D data matrix. If you wish to select the image files directly, click cancel.', 'Variable Name?', 1, {'P'});
if isempty(answer)
    return;
end
if isempty(answer{1})
    return;
end

searchContent = get(handles.searchContentEdit,'string');
if isempty(searchContent)
    searchContent = num2str(coord(1,3));
end
searchContent = ['ROI_' searchContent];
[filename, pathname] = uiputfile('*.mat', 'Save result to', searchContent);
if isequal(filename,0) | isequal(pathname,0)
   return
else
   thisfilename = fullfile(pathname, filename);
end


P = evalin('base',answer{1}); % data
try
    E = evalin('base',answer{2}); % event
catch
    E = [];
end
try
    C = evalin('base',answer{3}); % correlator
catch
    C = [];
end
try
    headmovement = evalin('base',answer{4}); % headmovement
catch
    headmovement = [];
end

mni = cell2mat(handles.currentDisplayMNI');
coord = mni2cor(mni, handles.M{1});


if ischar(P)
    h = waitbar(0,'Please wait...');
    for ii=1:size(P,1)
        tmp = load(deblank(P(ii,:)));
        tmp1 = fields(tmp);
        eval(['M = tmp.' tmp1{1} ';']);
        if isnumeric(M)
            M = double(M);
        end
        if ~isempty(E)
            tmp = load(deblank(E(ii,:)));
            tmp1 = fields(tmp);
            eval(['event{ii} = tmp.' tmp1{1} ';']);        
            if ~isstruct(event{ii})
                error('I need you specify the event in a structure, not a cell array!');
                return
            end
        else
            event = {};
        end
     
        if ~isempty(headmovement) % if headmovement is present (or any other effects to be removed), then get rid of it first
            H = [];
            tmp = load(deblank(headmovement(ii,:)));  
            fd=fields(tmp); 
            eval(['thisheadmovement=tmp.' fd{1} ';']);
            
            tmphd = struct2cell(thisheadmovement);
            for kk=1:length(tmphd)
                [r,c] = size(tmphd{kk});
                if r==1; tmphd{kk} = tmphd{kk}'; end
                H = [H  tmphd{kk}];
            end    
            for kk=1:size(coord,1)
                b = linearregression(squeeze(M(coord(kk,1), coord(kk,2), coord(kk,3), :)),H);
                M(coord(kk,1), coord(kk,2), coord(kk,3), :) = squeeze(M(coord(kk,1), coord(kk,2), coord(kk,3), :)) - H*b(1:end-1);
            end
        end
        
        [eventResponse{ii}, time{ii}, wholeOriginal{ii}, wholeAbsolute{ii}, wholeBaseline{ii}, removedEvent{ii}, remainedEvent{ii}] = cuixuBOLDretrieve(M, coord, event{ii}, 50, -20, 100, 0, TR_);

        if ~isempty(C)
            tmp = load(deblank(C(ii,:)));
            tmp1 = fields(tmp);
            eval(['correlator{ii} = tmp.' tmp1{1} ';']);        
            if ~isstruct(correlator{ii})
                error('I need you specify the correlator in a structure, not a cell array!');
                return
            end

            names = fields(correlator{ii});
            for kk=1:length(names)
                eval(['names2 = fields(correlator{ii}.' names{kk} ');']);
                for mm=1:length(names2)
                    eval(['correlator{ii}.' names{kk} '.' names2{mm} ' = correlator{ii}.' names{kk} '.' names2{mm} '(remainedEvent{ii}.' names{kk} ');']);;
                end
            end
        end
        waitbar(ii/size(P,1),h,['finished ' num2str(ii) ' of ' num2str(size(P,1))]);
        disp(['finished ' num2str(ii) ' of ' num2str(size(P,1))]);
    end
    close(h)    
elseif isnumeric(P) % never used!
    P = double(P);
	if ~isempty(E)
        tmp = load(deblank(E(ii,:)));
        tmp1 = fields(tmp);
        eval(['event{ii} = tmp.' tmp1{1} ';']);        
	else
        event = {};
	end
    [eventResponse{ii}, time{ii}, wholeOriginal{ii}, wholeAbsolute{ii}, wholeBaseline{ii}, removedEvent{ii}, remainedEvent{ii}] = cuixuBOLDretrieve(P, coord, event{ii}, 50, -20, 100, 0, TR_);
end

if ~exist('correlator')
    correlator = {};
end

mat = handles.M{1};
save(thisfilename, 'eventResponse', 'time', 'wholeOriginal', 'wholeAbsolute', 'wholeBaseline', 'correlator', 'coord', 'mni', 'removedEvent','remainedEvent', 'event', 'mat', 'P', 'headmovement');
CallBack_plotROI(hObject, eventdata, thisfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% myerrorbar plot (change errorbar width)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hE = myerrorbar(varargin)
hE = errorbar(varargin{:});
xdata = get(hE, 'xdata');
errorbarwidth = (xdata(2) - xdata(1))/4;
% adjust error bar width
hE_c                   = ...
    get(hE     , 'Children'    );
errorbarXData          = ...
    get(hE_c(2), 'XData'       );
errorbarXData(4:9:end) = ...
    errorbarXData(1:9:end) - errorbarwidth;
errorbarXData(7:9:end) = ....
    errorbarXData(1:9:end) - errorbarwidth;
errorbarXData(5:9:end) = ...
    errorbarXData(1:9:end) + errorbarwidth;
errorbarXData(8:9:end) = ...
    errorbarXData(1:9:end) + errorbarwidth;
set(hE_c(2), 'XData', errorbarXData);        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_plotROI(hObject, eventdata, filename)
handles = guidata(hObject);

if ~exist('filename')
    if findstr('SPM2',spm('ver'))
        P = spm_get([0:1],'*.mat','Select ROI result file');
    else %if findstr('SPM5',spm('ver'))
        P = spm_select([0:1],'mat','Select ROI result file');
    end
else
    P = filename;
end

if isempty(P)
    return;
end
load(deblank(P));

if exist('mni')
    delete(gcf)
    xjview(mni);
end


% allow user to select which subject to plot
prompt=['I find you have ' num2str(length(eventResponse)) ' subjects.  Please let me know which subjects you want to include to generate the plots.'];
name='Which subject(s) to plot?';
numlines=1;
defaultanswer={['[1:'  num2str(length(eventResponse)) ']']};
 
whichtoplot=inputdlg(prompt,name,numlines,defaultanswer);
if isempty(whichtoplot)
    return;
end
whichtoplot = eval(whichtoplot{1});


response = eventResponse{1};

eventisstruct = 0;
if isstruct(response)
    eventisstruct = 1;
	f = fieldnames(response);
    response = struct2cell(response); 
else
    f = num2cell([1:length(response)]);
    for ii=1:length(response)
        f{ii} = num2str(f{ii});
    end
end
for jj=1:length(response)
    response{jj} = [];
end

if exist('correlator') & ~isempty(correlator) & eventisstruct==1
    C = correlator{1};
    cenames = fields(C);
    for jj=1:length(cenames)
        eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
        for kk=1:length(correlatornames{jj})
            eval(['C.' cenames{jj} '.' correlatornames{jj}{kk} '= [];']);
        end
    end
    for ii=whichtoplot%1:length(correlator)
        for jj=1:length(cenames)
            eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
            for kk=1:length(correlatornames{jj})
                tmp2 = eval(['correlator{ii}.' cenames{jj} '.' correlatornames{jj}{kk} ';']);
                if size(tmp2,1) == 1
                    tmp2 = tmp2';
                end
                eval(['C.' cenames{jj} '.' correlatornames{jj}{kk} '= [C.' cenames{jj} '.' correlatornames{jj}{kk} '; tmp2];']);
            end
        end
    end
end

if eventisstruct == 1
    for ii=1:length(eventResponse)
        for kk=1:length(f)
            eval(['tmptmp{kk} = eventResponse{ii}.' f{kk} ';']);
        end
        eventResponse{ii} = tmptmp;
    end
end

for ii=whichtoplot%1:length(eventResponse)
    for jj=1:length(eventResponse{1})
        response{jj} = [response{jj}; eventResponse{ii}{jj}];
    end
end

if isempty(response)
    msgbox('Nothing to plot');
    return
end




v = version;
figure;
%colors='rbkmcgyrbkmcgyrbkmcgyrbkmcgyrbkmcgyrbkmcgyrbkmcgy';
colors=[1 0 0;
    0 0 1;
    0 0 0;
    1 0 1;
    0 1 1;
    0 1 0;
    1 1 0;
    1 0.5 0.5;
    .5 .5 1;
    .5 .5 .5;
    1 .5 0;
    1 0 .5;
    0 1 .5;
    .5 1 0;
    0 .5 1;
    .5 0 1];
colors = repmat(colors, 10, 1);
thislabel = {};
for ii=1:length(response)
    if v(1)=='6'
        ff = myerrorbar(time{1}, meannan(response{ii}), stdnan(response{ii})/sqrt(size(response{ii},1)));
        set(ff,'color',colors(ii, :));
        thislabel = [thislabel {''} f(ii)];
    elseif v(1)=='7'
        hE = myerrorbar(time{1}, meannan(response{ii}), stdnan(response{ii})/sqrt(size(response{ii},1)), 'color', colors(ii, :));
%         errorbarwidth = (time{1}(2) - time{1}(1))/4;
%         % adjust error bar width
%         hE_c                   = ...
%             get(hE     , 'Children'    );
%         errorbarXData          = ...
%             get(hE_c(2), 'XData'       );
%         errorbarXData(4:9:end) = ...
%             errorbarXData(1:9:end) - errorbarwidth;
%         errorbarXData(7:9:end) = ....
%             errorbarXData(1:9:end) - errorbarwidth;
%         errorbarXData(5:9:end) = ...
%             errorbarXData(1:9:end) + errorbarwidth;
%         errorbarXData(8:9:end) = ...
%             errorbarXData(1:9:end) + errorbarwidth;
%         set(hE_c(2), 'XData', errorbarXData);        
        thislabel = [thislabel f(ii)];
    end
    hold on;
end
xlabel('time (s)');
ylabel('relative signal');
legend(thislabel)
hold off;

[row,col]  = size(response);
if row == 1; response = response'; end;
ResponseForPlot = cell2struct(response, f, 1);

if exist('correlator') & ~isempty(correlator)
    %list all plot names in command window for selection
    promt{1} = 'I will plot signal amplitude vs correlator. Here are all the correlators:';
    promt{2} = 'index   event    correlator';
    count = 1;
    for jj=1:length(cenames)
        eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
        for kk=1:length(correlatornames{jj})
            promt{count+2} = sprintf('%d         %s         %s', count,  cenames{jj}, correlatornames{jj}{kk});
            count = count + 1;
        end
    end
    promt{count+2} = 'Which plots do you want me to show?';

    promt = char(promt);
    name='Please select which plots to show';
    numlines=1;
    defaultanswer={['[1:'  num2str(count-1) ']']};

    whichplot=inputdlg(promt,name,numlines,defaultanswer);
    if isempty(whichplot)
        return;
    end    
    whichplot = eval(whichplot{1});
    
    count = 1;
    for jj=1:length(cenames)
        eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
        for kk=1:length(correlatornames{jj})
            if ~ismember(count, whichplot)
                count = count + 1;
                continue;
            end
            
            eval(['x = C.' cenames{jj} '.' correlatornames{jj}{kk} ';']);
            
            uniquex = unique(x);
            pos = find(isnan(uniquex));
            uniquex(pos) = [];
            
            if(length(uniquex)<10)  % only if the correlator is discrete, then plot BOLD vs time for each correlator
                figure;
                thislabel = {};
                for ii=1:length(uniquex)
                    pos = find(x==uniquex(ii));
                    eval(['z = meannan(ResponseForPlot.' cenames{jj} '(pos,:));']);
                    eval(['stdz = stdnan(ResponseForPlot.' cenames{jj} '(pos,:))/sqrt(length(pos));']);

                    if v(1)=='6'
                        ff = myerrorbar(time{1}, z, stdz);
                        set(ff,'color',colors(ii, :));
                        thislabel = [thislabel {''} {(num2str(uniquex(ii)))}];
                    elseif v(1)=='7'
                        myerrorbar(time{1}, z, stdz, 'color', colors(ii, :));
                        thislabel = [thislabel {(num2str(uniquex(ii)))}];
                    end
                    hold on;
                end
                legend(thislabel);
                xlabel('time in second');
                ylabel('BOLD');
                title([cenames{jj} ' ' correlatornames{jj}{kk}])
            end
            
            count = count + 1;
        end
    end

    peakpoint = inputdlg({'I will plot peak signal. What points in time do you considered as peak signal?'}, 'Peak point',1,{'4 6'});
    if isempty(peakpoint); return; end
    peakpoint = str2num(peakpoint{1});
    
    
	%list all plot names in command window for selection, effect(s) to be
	%removed
    clear promt;
    promt{1} = 'Which effect(s) do you want me to remove first? (If none to remove, leave it blank and click OK)';
    promt{2} = 'index   event    correlator';
    count = 1;
    for jj=1:length(cenames)
        eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
        for kk=1:length(correlatornames{jj})
            promt{count+2} = sprintf('%d         %s         %s', count,  cenames{jj}, correlatornames{jj}{kk});
            count = count + 1;
        end
    end
    %promt{count+2} = 'Which effect(s) do you want me to remove first?';

    promt = char(promt);
    name='Please select which effect to remove';
    numlines=1;
    defaultanswer={''};

	whichremovefirst=inputdlg(promt,name,numlines,defaultanswer);
    if isempty(whichremovefirst)
        whichremovefirst = [];
    elseif isempty(deblank(whichremovefirst{1}))
        whichremovefirst = [];
    else
        whichremovefirst = eval(whichremovefirst{1});
    end
    
   
%% remove effect first
    tmpResponseForPlot = ResponseForPlot;
    for jj=1:length(cenames)
        eval(['tmpResponseForPlot.' cenames{jj} ' = mean(tmpResponseForPlot.' cenames{jj} '(:,findind(time{1},peakpoint)), 2);']);
    end
    
    if ~isempty(whichremovefirst)
        count = 1;
        for jj=1:length(cenames)
            eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
            eval(['y = tmpResponseForPlot.' cenames{jj} ';']);
            x = [];
            if size(y,1)==1
                y = y';
            end
            for kk=1:length(correlatornames{jj})
                if ~ismember(count, whichremovefirst)
                    count = count + 1;
                    continue;
                end

                eval(['tmpx = C.' cenames{jj} '.' correlatornames{jj}{kk} ';']);

                % convert to colum vector
                if size(tmpx,1)==1
                    tmpx = tmpx';
                end
                x = [x tmpx];



                count = count + 1;
            end

            % remove NaN
            tmpposx = find(isnan(mean(x,2)) | isinf(mean(x,2)));
            tmpposy = find(isnan(mean(y,2)) | isinf(mean(y,2)));
            x([tmpposx; tmpposy],:) = [];
            y([tmpposx; tmpposy],:) = [];

            b = linearregression(y,x);
            y = y - x*b(1:end-1) - b(end);
            eval(['tmpResponseForPlot.' cenames{jj} ' = y;']);
        end
    end

%% now plot the residual amplitude against correlator
    count = 1;
    for jj=1:length(cenames)
        eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
        for kk=1:length(correlatornames{jj})
            if ~ismember(count, whichplot)
                count = count + 1;
                continue;
            end

            eval(['y = tmpResponseForPlot.' cenames{jj} ';']);
            eval(['x = C.' cenames{jj} '.' correlatornames{jj}{kk} ';']);
            
           
            figure
            plot(x,y, 'sb');
            hold on;            
            plot2(x,y,'r',1);
            
            xlabel(correlatornames{jj}{kk})
            ylabel('peak relative BOLD');
            legend(cenames{jj})
            try
                % convert to colum vector
                if size(x,1)==1
                    x = x';
                end
                if size(y,1)==1
                    y = y';
                end
                
                % remove NaN
                tmpposx = find(isnan(x) | isinf(x));
                tmpposy = find(isnan(y) | isinf(y));
                x([tmpposx; tmpposy]) = [];
                y([tmpposx; tmpposy]) = [];

                b = linearregression(y,x);
                
                totalvar = var(y);
                residvar = var(y - x*b(1) - b(2));
                modelvar = totalvar - residvar;
                dfm = 1;
                dft = length(x) - 1;
                dfr = dft - dfm;
                mmodelvar = modelvar/dfm;
                mresidvar = residvar/dfr;
                F = mmodelvar/mresidvar;
                pvalue = 1-spm_Fcdf(F, [dfm dfr]);
                
                title(sprintf('pValue=%s; slope=%g; constant=%g', pvalue, b(1), b(2)))
                xx=[min(x):(max(x)-min(x))/10:max(x)];
                yy=b(1) * xx + b(2);
                hold on;
                plot(xx,yy,'g');
            end
            count = count + 1;
        end
    end
end

function index = findind(vector1, vector2)
kk = 1;
for ii=1:length(vector2)
    tmp = find(vector1 == vector2(ii));
    if ~isempty(tmp)
        index(kk) = tmp;
        kk = kk+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot ROI, individual by individual, with behavior timing on it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_plotIndividualROIWithBehavior(hObject, eventdata)
global TR_;
handles = guidata(hObject);

if findstr('SPM2',spm('ver'))
    P = spm_get([0:1],'*.mat','Select ROI result file');
else%if findstr('SPM5',spm('ver'))
    P = spm_select([0:1],'mat','Select ROI result file');
end

if isempty(P)
    return;
end
load(P);
if ~exist('TR')
    TR = TR_;
end

colors=[1 0 0;
    0 0 1;
    0 0 0;
    1 0 1;
    0 1 1;
    0 1 0;
    1 1 0;
    1 0.5 0.5;
    .5 .5 1;
    .5 .5 .5;
    1 .5 0;
    1 0 .5;
    0 1 .5;
    .5 1 0;
    0 .5 1;
    .5 0 1];
colors = repmat(colors, 10, 1);


figure;
Maximize(gcf);
for ii=1:length(event)
    kk = mod(ii,8);
    if kk == 0
        kk = 8;
    end
    subplot(4,2,kk)

    m = mean(wholeOriginal{ii});
    s = std(wholeOriginal{ii});
    mx = max(wholeOriginal{ii});
    mn = min(wholeOriginal{ii});
    ylowerlimit = mn;
    
    hm = event{ii};
    hmcell = struct2cell(hm);
    hmname = fieldnames(hm);
    for jj=1:length(hmname)
        toplot1 = [];
        toplot2 = [];
        for mm = 1:length(hmcell{jj})
            toplot1 = [toplot1 hmcell{jj}(mm) hmcell{jj}(mm) nan];
            toplot2 = [toplot2 m-s m+s nan];
        end
        plot(toplot1, toplot2, 'color', colors(jj,:));
        hold on
    end
    legend(hmname)
    title(['subject ' num2str(ii)])
    xlabel('time in s');

    % plot imaging data
    plot(TR*[0:length(wholeOriginal{ii})-1], wholeOriginal{ii}, 'k');
    hold on;
    ylim([ylowerlimit mx])

    if mod(ii,8) == 0 && ii~=length(event)
        figure;
        Maximize(gcf);
    end     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot ROI, individual by individual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_plotIndividualROI(hObject, eventdata, filename)
handles = guidata(hObject);

if ~exist('filename')
    if findstr('SPM2',spm('ver'))
        P = spm_get([0:1],'*.mat','Select ROI result file');
    else%if findstr('SPM5',spm('ver'))
        P = spm_select([0:1],'mat','Select ROI result file');
    end
else
    P = filename;
end

if isempty(P)
    return;
end
load(deblank(P));

if exist('mni')
    delete(gcf)
    xjview(mni);
end

response = eventResponse{1};

eventisstruct = 0;
if isstruct(response)
    eventisstruct = 1;
	f = fieldnames(response);
    response = struct2cell(response); 
else
    f = num2cell([1:length(response)]);
    for ii=1:length(response)
        f{ii} = num2str(f{ii});
    end
end
for jj=1:length(response)
    response{jj} = [];
end

if exist('correlator') & ~isempty(correlator) & eventisstruct==1
    C = correlator{1};
    cenames = fields(C);
    for jj=1:length(cenames)
        eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
        for kk=1:length(correlatornames{jj})
            eval(['C.' cenames{jj} '.' correlatornames{jj}{kk} '= [];']);
        end
    end
    for ii=1:length(correlator)
        for jj=1:length(cenames)
            eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
            for kk=1:length(correlatornames{jj})
                tmp2 = eval(['correlator{ii}.' cenames{jj} '.' correlatornames{jj}{kk} ';']);
                if size(tmp2,1) == 1
                    tmp2 = tmp2';
                end
                eval(['C.' cenames{jj} '.' correlatornames{jj}{kk} '= [C.' cenames{jj} '.' correlatornames{jj}{kk} '; tmp2];']);
            end
        end
    end
end

if eventisstruct == 1
    for ii=1:length(eventResponse)
        for kk=1:length(f)
            eval(['tmptmp{kk} = eventResponse{ii}.' f{kk} ';']);
        end
        eventResponse{ii} = tmptmp;
    end
end

for ii=1:length(eventResponse)
    for jj=1:length(eventResponse{1})
        response{jj} = [response{jj}; eventResponse{ii}{jj}];
    end
end

if isempty(response)
    msgbox('Nothing to plot');
    return
end




v = version;
figure;
%colors='rbkmcgyrbkmcgyrbkmcgyrbkmcgyrbkmcgyrbkmcgyrbkmcgy';
colors=[1 0 0;
    0 0 1;
    0 0 0;
    1 0 1;
    0 1 1;
    0 1 0;
    1 1 0;
    1 0.5 0.5;
    .5 .5 1;
    .5 .5 .5;
    1 .5 0;
    1 0 .5;
    0 1 .5;
    .5 1 0;
    0 .5 1;
    .5 0 1];
colors = repmat(colors, 10, 1);
thislabel = {};
for ii=1:length(response)
    if v(1)=='6'
        ff = myerrorbar(time{1}, meannan(response{ii}), stdnan(response{ii})/sqrt(size(response{ii},1)));
        set(ff,'color',colors(ii, :));
        thislabel = [thislabel {''} f(ii)];
    elseif v(1)=='7'
        myerrorbar(time{1}, meannan(response{ii}), stdnan(response{ii})/sqrt(size(response{ii},1)), 'color', colors(ii, :));
        thislabel = [thislabel f(ii)];
    end
    hold on;
end
xlabel('time (s)');
ylabel('relative signal');
legend(thislabel)
hold off;

figure;
Maximize(gcf);
for jj=1:length(eventResponse)
    wheretoplot = mod(jj, 16);
    if wheretoplot == 0; 
        wheretoplot = 16; 
        subplot(4,4,wheretoplot)
    else
        subplot(4,4,wheretoplot)
    end
    for ii=1:length(response)
        if v(1)=='6'
            ff = myerrorbar(time{1}, meannan(eventResponse{jj}{ii}), stdnan(eventResponse{jj}{ii})/sqrt(size(eventResponse{jj}{ii},1)));
            set(ff,'color',colors(ii, :));
            %thislabel = [thislabel {''} f(ii)];
        elseif v(1)=='7'
            myerrorbar(time{1}, meannan(eventResponse{jj}{ii}), stdnan(eventResponse{jj}{ii})/sqrt(size(eventResponse{jj}{ii},1)), 'color', colors(ii, :));
            %myerrorbar(time{1}, meannan(response{ii}), stdnan(response{ii})/sqrt(size(response{ii},1)), 'color', colors(ii, :));
            %thislabel = [thislabel f(ii)];
        end
        hold on;
    end
    title(['subject ' num2str(jj)])
    
    if mod(jj, 16)==0 && jj~=length(eventResponse)
        figure;
        Maximize(gcf);
    end
end



[row,col]  = size(response);
if row == 1; response = response'; end;
for jj=1:length(eventResponse)
    ResponseForPlot{jj} = cell2struct(eventResponse{jj}, f, 2);
end

if exist('correlator') & ~isempty(correlator)
    %list all plot names in command window for selection
    promt{1} = 'I will plot signal amplitude vs correlator. Here are all the correlators:';
    promt{2} = 'index   event    correlator';
    count = 1;
    for jj=1:length(cenames)
        eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
        for kk=1:length(correlatornames{jj})
            promt{count+2} = sprintf('%d         %s         %s', count,  cenames{jj}, correlatornames{jj}{kk});
            count = count + 1;
        end
    end
    promt{count+2} = 'Which plots do you want me to show?';

    promt = char(promt);
    name='Please select which plots to show';
    numlines=1;
    defaultanswer={['[1:'  num2str(count-1) ']']};

    whichplot=inputdlg(promt,name,numlines,defaultanswer);
    if isempty(whichplot)
        return;
    end    
    whichplot = eval(whichplot{1});

    peakpoint = inputdlg({'What points in time do you considered as peak signal?'}, 'Peak point',1,{'4 6'});
    if isempty(peakpoint); return; end
    peakpoint = str2num(peakpoint{1});   
    
    count = 1;
    for jj=1:length(cenames)
        eval(['correlatornames{jj} = fields(C.' cenames{jj} ');']);
        for kk=1:length(correlatornames{jj})
            if ~ismember(count, whichplot)
                count = count + 1;
                continue;
            end
            figure;
            Maximize(gcf);
            for ii=1:length(eventResponse)
                wheretoplot = mod(ii, 16);
                if wheretoplot == 0;
                    wheretoplot = 16;
                    subplot(4,4,wheretoplot)
                else
                    subplot(4,4,wheretoplot)
                end


                eval(['y = mean(ResponseForPlot{ii}.' cenames{jj} '(:,findind(time{1},peakpoint)), 2);']);
                eval(['x = correlator{ii}.' cenames{jj} '.' correlatornames{jj}{kk} ';']);

                %plot(x,y, 'sb');
                %hold on;            
                pos = find(isnan(x));
                x(pos) = [];
                y(pos) = [];
                pos = find(isnan(y));
                x(pos) = [];
                y(pos) = [];                
                if length(x) < 1
                    continue;
                end
                plot2(x,y,'r',1);

                xlabel(correlatornames{jj}{kk})
                ylabel(['peak signal at ' cenames{jj}]);
                %legend(cenames{jj})
                try
                    % convert to colum vector
                    if size(x,1)==1
                        x = x';
                    end
                    if size(y,1)==1
                        y = y';
                    end

                    % remove NaN
                    tmpposx = find(isnan(x) | isinf(x));
                    tmpposy = find(isnan(y) | isinf(y));
                    x([tmpposx; tmpposy]) = [];
                    y([tmpposx; tmpposy]) = [];

                    b = linearregression(y,x);

                    totalvar = var(y);
                    residvar = var(y - x*b(1) - b(2));
                    modelvar = totalvar - residvar;
                    dfm = 1;
                    dft = length(x) - 1;
                    dfr = dft - dfm;
                    mmodelvar = modelvar/dfm;
                    mresidvar = residvar/dfr;
                    F = mmodelvar/mresidvar;
                    pvalue = 1-spm_Fcdf(F, [dfm dfr]);

                    title(sprintf('subject %d, pValue=%s', ii, pvalue))
                    xx=[min(x):(max(x)-min(x))/10:max(x)];
                    yy=b(1) * xx + b(2);
                    hold on;
                    plot(xx,yy,'g');
                end
                if mod(ii, 16)==0 && ii~=length(eventResponse)
                    figure;
                    Maximize(gcf);
                end
            end % end of loop subject
            count = count + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot ROI: correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_plotCorrelationROI(hObject, eventdata)

if findstr('SPM2',spm('ver'))
    P1 = spm_get([0:100],'*.mat','Select ROI result file(s)');
else%if findstr('SPM5',spm('ver'))
    P1 = spm_select([0:100],'mat','Select ROI result file(s)');
end

if isempty(P1)
    return;
end

N = size(P1,1);
t = [-20:20];

v = version;

if N == 1
    load(deblank(P1));
    tmp = [];
    for jj=1:length(wholeAbsolute)
        tmp = [tmp; corrlag(wholeAbsolute{jj},wholeAbsolute{jj},t)];
    end
    figure;
    myerrorbar(t, mean(tmp,1), std(tmp)/sqrt(length(wholeAbsolute)));
    xlabel('lag in scan number');
    ylabel('auto correlation');
    title(deblank(P1))
    return
else
    figure;
    colors='rbkmcgyrbkmcgyrbkmcgyrbkmcgyrbkmcgyrbkmcgyrbkmcgy';
    legendlabel = {};
    totalPermute = nchoosek([1:N],2);
    for ii=1:size(totalPermute,1)
        tmp = [];
        s1 = load(deblank(P1(totalPermute(ii,1),:)));
        s2 = load(deblank(P1(totalPermute(ii,2),:)));        
        for jj=1:length(s1.wholeAbsolute)
            tmp = [tmp; corrlag(s1.wholeAbsolute{jj},s2.wholeAbsolute{jj},t)];
        end
        myerrorbar(t, mean(tmp,1), std(tmp)/sqrt(length(s1.wholeAbsolute)), colors(ii));
        xlabel('lag in scan number');
        ylabel('cross correlation');
        hold on;
        if v(1) == '6'
            legendlabel = [legendlabel, {'', P(totalPermute(ii,:),:)}];
        elseif v(1) == '7'
            legendlabel = [legendlabel, {P1(totalPermute(ii,:),:)}];
        end
    end
    legend(legendlabel);
    hold off;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% whole brain correlation analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_wholeBrainCorrelation(hObject, eventdata)
% answer = getVariable({'Which data file?'});
% if isempty(answer)
%     return;
% end
% if isempty(answer{1})
%     return;
% end
% P = evalin('base',answer{1}); % subject directories
if findstr('SPM2',spm('ver'))
    P = spm_get([0:1],'*.mat','Whole Brain Signal file');
else%if findstr('SPM5',spm('ver'))
    P = spm_select([0:1],'mat','Whole Brain Signal file');
end

if isempty(P)
    return
end

load(deblank(P));
M = double(M);

N=zeros(size(M,4), size(M,1) * size(M,2) * size(M,3));

mm = 1;
for ii=1:2:size(M,1)
    for jj=1:3:size(M,2)
        for kk=1:2:size(M,3)
            N(:, mm) = squeeze(M(ii,jj,kk,:));
            mm = mm + 1;
            cor(mm,:) = [ii jj kk];
        end
    end
end

N(:, mm:end) = [];

C=corrcoef(N);
C(find(isnan(C))) = 0;

for ii=1:4
    CC{ii} = abs(C)>(0.5+ii/10);
end

K = {zeros(1, size(C,1)),zeros(1, size(C,1)),zeros(1, size(C,1)),zeros(1, size(C,1))};
for ii=1:size(C,1)
    for jj=1:4
        K{jj}(ii) = sum(CC{jj}(ii,:)) - 1;
    end
end

x=[0:size(C,1)];
for jj=1:4
    y{jj}=zeros(size(x));
end

for ii=1:length(x)
    for jj=1:4
        y{jj}(ii) = sum(K{jj} == x(ii));
    end
end

figure;
loglog(x,y{1}, x, y{2}, x, y{3}, x, y{4});
xlabel('degree');
ylabel('counts');
legend('0.6', '0.7', '0.8', '0.9');
axis equal


z = sum(abs(C))-1;
zx = 0:max(z);
for ii=1:length(zx)-1
    pos = find(z>zx(ii) & z<=zx(ii+1));
    zy(ii) = length(pos);
end
figure;
loglog(zx(1:length(zy)),zy);
xlabel('weighted degree');
ylabel('counts');
axis equal

xjview(cor2mni(coord, handles.M{1}), z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get/select variable names in base workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vars = getVariable(titles)
vars = evalin('base','who');
height = 0.07;
f = dialog('unit','normalized', 'menubar','none', 'position', [0.3 0.2 0.4 0.4], 'name', 'pick variables', 'NumberTitle','off');
uicontrol(f,'style','text',...
        'unit','normalized',...
        'String', 'Available variables', ...
        'position',[0 0.9 0.5 0.1/2]);    
variableListbox = uicontrol(f,'style','listbox','tag','variableListbox',...
        'unit','normalized',...
        'String', vars, ...
        'value',1,...
        'position',[0 0 0.5 0.9]);
okPush = uicontrol(f,'style','push',...
        'unit','normalized',...
        'String', 'OK', ...
        'position',[0.55 0.05 0.2 height],...
        'callback','uiresume');
cancelPush = uicontrol(f,'style','push',...
        'unit','normalized',...
        'String', 'Cancel', ...
        'position',[0.75 0.05 0.2 height],...
        'callback','delete(gcf)');  

uicontrol(f,'style','text',...
        'unit','normalized',...
        'String', titles{1}, ...
        'position',[0.5 0.9 0.5 0.1/2]);   
imageDataPush = uicontrol(f,'style','push',...
        'unit','normalized',...
        'String', '->', ...
        'position',[0.5 0.8 0.1 height],...
        'callback', [...
            'xyouneverfindme=findobj(get(gcf,''Children''),''flat'',''tag'',''imageDataEdit'');'...
            'yyouneverfindme=findobj(get(gcf,''Children''),''flat'',''tag'',''variableListbox'');'...
            'allyouneverfindme=get(yyouneverfindme,''string'');'...
            'selectedyouneverfindme=allyouneverfindme{get(yyouneverfindme,''value'')};'...
            'set(xyouneverfindme,''string'',selectedyouneverfindme);'...
            'clear allyouneverfindme selectedyouneverfindme xyouneverfindme yyouneverfindme;']);      
imageDataEdit = uicontrol(f,'style','edit','tag','imageDataEdit',...
        'unit','normalized',...
        'String', '', ...
        'BackgroundColor', 'w',...
        'position',[0.6 0.8 0.3 height]);

if length(titles)>=2
uicontrol(f,'style','text',...
        'unit','normalized',...
        'String', titles{2}, ...
        'position',[0.5 0.7 0.5 0.1/2]);   
eventDataPush = uicontrol(f,'style','push',...
        'unit','normalized',...
        'String', '->', ...
        'position',[0.5 0.6 0.1 height],...
        'callback', [...
            'xyouneverfindme=findobj(get(gcf,''Children''),''flat'',''tag'',''eventDataEdit'');'...
            'yyouneverfindme=findobj(get(gcf,''Children''),''flat'',''tag'',''variableListbox'');'...
            'allyouneverfindme=get(yyouneverfindme,''string'');'...
            'selectedyouneverfindme=allyouneverfindme{get(yyouneverfindme,''value'')};'...
            'set(xyouneverfindme,''string'',selectedyouneverfindme);'...
            'clear allyouneverfindme selectedyouneverfindme xyouneverfindme yyouneverfindme;']);           

eventDataEdit = uicontrol(f,'style','edit','tag','eventDataEdit',...
        'unit','normalized',...
        'String', '', ...
        'BackgroundColor', 'w',...
        'position',[0.6 0.6 0.3 height]);    
end
if length(titles)>=3
uicontrol(f,'style','text',...
        'unit','normalized',...
        'String', titles{3}, ...
        'position',[0.5 0.5 0.5 0.1/2]);   
correlatorDataPush = uicontrol(f,'style','push',...
        'unit','normalized',...
        'String', '->', ...
        'position',[0.5 0.4 0.1 height],...
        'callback', [...
            'xyouneverfindme=findobj(get(gcf,''Children''),''flat'',''tag'',''correlatorDataEdit'');'...
            'yyouneverfindme=findobj(get(gcf,''Children''),''flat'',''tag'',''variableListbox'');'...
            'allyouneverfindme=get(yyouneverfindme,''string'');'...
            'selectedyouneverfindme=allyouneverfindme{get(yyouneverfindme,''value'')};'...
            'set(xyouneverfindme,''string'',selectedyouneverfindme);'...
            'clear allyouneverfindme selectedyouneverfindme xyouneverfindme yyouneverfindme;']);           
    
correlatorDataEdit = uicontrol(f,'style','edit','tag','correlatorDataEdit',...
        'unit','normalized',...
        'String', '', ...
        'BackgroundColor', 'w',...
        'position',[0.6 0.4 0.3 height]); 
end
if length(titles)>=4
uicontrol(f,'style','text',...
        'unit','normalized',...
        'String', titles{4}, ...
        'position',[0.5 0.3 0.5 0.1/2]);   
otherDataPush = uicontrol(f,'style','push',...
        'unit','normalized',...
        'String', '->', ...
        'position',[0.5 0.2 0.1 height],...
        'callback', [...
            'xyouneverfindme=findobj(get(gcf,''Children''),''flat'',''tag'',''otherDataEdit'');'...
            'yyouneverfindme=findobj(get(gcf,''Children''),''flat'',''tag'',''variableListbox'');'...
            'allyouneverfindme=get(yyouneverfindme,''string'');'...
            'selectedyouneverfindme=allyouneverfindme{get(yyouneverfindme,''value'')};'...
            'set(xyouneverfindme,''string'',selectedyouneverfindme);'...
            'clear allyouneverfindme selectedyouneverfindme xyouneverfindme yyouneverfindme;']);           
   
otherDataEdit = uicontrol(f,'style','edit','tag','otherDataEdit',...
        'unit','normalized',...
        'String', '', ...
        'BackgroundColor', 'w',...
        'position',[0.6 0.2 0.3 height]); 
end

uiwait(f);
try
    var{1} = get(imageDataEdit,'string');
    if length(titles)>=2
        var{2} = get(eventDataEdit,'string');
    end
    if length(titles)>=3
        var{3} = get(correlatorDataEdit,'string');    
    end
    if length(titles)>=4
        var{4} = get(otherDataEdit,'string');    
    end
    vars = var;
    delete(f);
catch
    vars = {};
    try    
        delete(f);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% quit xjview
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_quit(hObject, eventdata, warnstate)
warning(warnstate)
% try
%     rmdir('xjviewtmp');
% end
delete(gcf);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mouse double click
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figureMouseUpFcn(hObject, eventdata)
status = get(hObject, 'SelectionType');

% double click
if strcmp(status, 'open')
    handles = guidata(hObject);
    CallBack_loadImagePush(handles.loadImagePush, eventdata);
else
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% change edit image file name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_imageFileEdit(hObject, eventdata)
handles = guidata(hObject);
filename = get(hObject, 'String');
filename = str2cell(filename);
CallBack_loadImagePush(handles.loadImagePush, eventdata, filename);

%%%%%%%%%%%%%%%%
function Callback_SetOverlayColor(hObject, eventdata)
%Change the overlay's color
%dawnsong, 20120415
nc=uisetcolor(hObject, get(get(hObject, 'UserData'), 'string'));
set(hObject,'BackGroundColor', nc);
set(get(hObject, 'UserData'),'ForeGroundColor', nc);
handles = guidata(hObject);
CallBack_slider(hObject, eventdata, -log10(handles.pValue));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% click load image file button
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_loadImagePush(hObject, eventdata, thisfilename)
handles = guidata(hObject);

set(handles.instructionText, 'visible', 'off');

handles.imageFileName=[]; handles.M=[]; handles.DIM=[]; handles.TF=[]; handles.df=[]; 
handles.mni=[]; handles.intensity=[]; handles.currentmni=[]; handles.currentintensity=[]; handles.currentDisplayMNI=[]; handles.currentDisplayIntensity=[];

if ~exist('thisfilename')
    thisfilename = ''; 
end

if isstruct(thisfilename)
    handles.imageFileName = {''};
    handles.mni = {thisfilename.mni};
	handles.intensity = {thisfilename.intensity};
    handles.M = {thisfilename.M};
    handles.DIM = {thisfilename.DIM};        
    handles.TF = {'S'};
    handles.df = {1};
    handles.pValue = 1;
    set(handles.pValueEdit, 'string', '1');
    handles.clusterSizeThreshold = 0;
    set(handles.clusterSizeThresholdEdit, 'String', '0');
else
    [handles.imageFileName, handles.M, handles.DIM, handles.TF, handles.df, handles.mni, handles.intensity] = getImageFile(thisfilename);
end

if isempty(handles.imageFileName)
    return
end

different = 0; % image files are same or different?
if length(handles.TF)>1
    for ii=1:length(handles.TF)
        if ( strcmp(handles.TF{ii},handles.TF{1}) || isempty(handles.TF{ii}) && isempty(handles.TF{1}) ) & isequal(handles.df{ii},handles.df{1}) & isequal(handles.M{ii},handles.M{1}) & isequal(handles.DIM{ii},handles.DIM{1})
            continue;
        else
            warndlg('Images are from different statistics or sources.', 'Warning');
            %set(handles.infoTextBox, 'string', 'Images are from different statistics or sources.'); 
            beep;
            different = 1;
            break;
        end
    end
end

% reset files with empty df/TF to df=1 and TF='S'. 'S'=='T' but has a tag
% meaning it is changed.
resetTF = 0;
for ii=1:length(handles.TF)
    if isempty(handles.TF{ii}) 
        handles.TF{ii} = 'S';
        resetTF = 1;
    end
    if isempty(handles.df{ii}) | isequal(handles.df{ii},0)
        handles.df{ii} = 1;
        resetTF = 1;
    end
end
    
handles.currentmni = handles.mni;
handles.currentintensity = handles.intensity;

set(handles.dfEdit, 'String', cell2str(handles.df));
set(handles.imageFileEdit, 'String', cell2str(handles.imageFileName)); % s=-log10(p)
maxs = maxcell(t2s(cellmax(handles.intensity,'abs'),handles.df, handles.TF),'abs'); 
if isinf(maxs); maxs = 20; end
set(handles.slider, 'Max', maxs, 'Min', 0, 'sliderstep',[min(maxs/100,0.05),min(maxs/100,0.05)]);
if handles.TF{1}=='T' & different == 0
    str = [blanks(length(' intensit')) 'T=']; 
elseif handles.TF{1}=='F' & different == 0
    str = [blanks(length(' intensit')) 'F=']; 
else
    str=' intensity='; 
end
set(handles.intensityThresholdText, 'String', str);
set(handles.figure,'Name',['xjView: ' cell2str(handles.imageFileName)]);
try
    for ii=1:length(handles.hLegend)
        delete(handles.hLegend{ii});
        delete(handles.hLegendBtn{ii});
    end
end
if length(handles.TF)>1
    colours = repmat([1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1],100,1); 
    for ii=1:min(length(handles.TF),10)
        [tmp,filename] = fileparts(handles.imageFileName{ii});
        if ii == 10; filename = '......'; end
        pos0 = handles.sectionViewPosition;
        pos(1) = pos0(1)+pos0(3)/2+0.03;
        pos(2) = pos0(2)+ii/50-0.03;
        pos(3) = 0.12;
        pos(4) = 0.02;
        handles.hLegend{ii}=uicontrol(handles.figure, 'style','text',...
        'unit','normalized','position',pos,...
        'string', filename,...
        'horizontal','left',...
        'fontweight','bold',...
        'ForeGroundColor',colours(ii,:));
        set(handles.hLegend{ii}, 'unit','pixel', ...
            'Tooltipstring', sprintf('Click left button to change Overlay color: %s', filename));
        pos=get(handles.hLegend{ii}, 'position');
        pos=[pos(1)-12,pos(2),12,12];
        handles.hLegendBtn{ii}=uicontrol(handles.figure, 'style','pushbutton',...
            'unit','pixel','position',pos,...
            'callback',@Callback_SetOverlayColor, ...
            'Tooltipstring', sprintf('Change Overlay color: %s', filename), ...
            'UserData', handles.hLegend{ii},...
            'BackGroundColor',colours(ii,:));
    end
	c_names = {'red';'yellow';'green';'cyan';'blue';'magenta'};
end
guidata(hObject, handles);

global M_;
global DIM_;
global XJVIEWURL_;
M_ = handles.M{1};
DIM_ = handles.DIM{1};

if resetTF==1
    set(handles.pValueEdit,'string',1);
end
CallBack_pValueEdit(handles.pValueEdit, eventdata);

% display info
% try
%     %urlread([XJVIEWURL_ '/guestbook/stat.php']);    
% 	s = urlread([XJVIEWURL_ '/toUser.txt']);
%     report{1} = s;
% catch
%     report{1} = 'Welcome to xjView 8';
% end
urlread([XJVIEWURL_ '/stat.php']);    

for jj=1:length(handles.imageFileName)
    report{2+6*(jj-1)} = cell2str(handles.imageFileName(jj));
    if handles.TF{jj} == 'T' | handles.TF{jj} == 'F'
        report{3+6*(jj-1)} = ['This is a ' handles.TF{jj} ' test image.'];
    else
        report{3+6*(jj-1)} = '';%['I don''t know what test this image came from.'];
    end
    report{4+6*(jj-1)} = 'mat = ';
    report{5+6*(jj-1)} = num2str(handles.M{jj});
    report{6+6*(jj-1)} = 'dimension = ';
    report{7+6*(jj-1)} = num2str(handles.DIM{jj});    
end
set(handles.infoTextBox, 'string', report);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pValueEdit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_pValueEdit(hObject, eventdata)
handles = guidata(hObject);

theVer =version;
theVer =str2num( theVer(1:3) );
if(theVer<=7.1),
  tmp=str2num(get(hObject, 'String'));
else
  %multi-p values allowed, dawnsong, 20120418
  mp=regexp(get(hObject, 'String'),';','split');
  tmp = str2double(mp);
end

%tmp = str2double(get(hObject, 'String'));
tmps = -log10(tmp);
if sum(isnan(tmp))>0 | tmp < 0 | tmp > 1
    errordlg('I don''t understand the input.','error');
    set(hObject, 'String', handles.pValue);
    return
end

if tmps>get(handles.slider,'max') | tmps<get(handles.slider,'min')
    errordlg('pValue is too small. No suprathreshold voxels.','error');
    set(hObject, 'String', handles.pValue);
    return
end

if isempty(handles.df) | handles.df{1}==0
    set(handles.pValueEdit,'String', 'NaN');
end

handles.pValue = tmp;
handles.intensityThreshold = p2t(num2cell(handles.pValue.*ones(1,length(handles.TF))), handles.df, handles.TF);

set(handles.pValueEdit,'string', cell2str(num2cell(handles.pValue)));
%set(handles.pValueEdit,'string', num2str(handles.pValue));
guidata(hObject, handles);
CallBack_slider(hObject, eventdata, -log10(handles.pValue));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% intensity threshold edit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_intensityThresholdEdit(hObject, eventdata)
handles = guidata(hObject);

%multi-p values allowed, dawnsong, 20120418
mp=regexp(get(hObject, 'String'),';','split');

tmp = str2double(mp);
if sum(isnan(tmp))>0 | tmp<0
    errordlg('Please input a single valid number bigger than 0.','error');
    return
end

if tmp > maxcell(cellmax(handles.intensity,'abs'))
    tmp = maxcell(cellmax(handles.intensity,'abs'));
end
for i=1:length(tmp),
    handles.pValue = t2p(tmp, handles.df{i}, handles.TF{i});
end
%handles.pValue = t2p(tmp, handles.df{1}, handles.TF{1});
handles.intensityThreshold = p2t(num2cell(handles.pValue.*ones(1,length(handles.TF))), handles.df, handles.TF);
set(handles.slider,'Value', -log10(handles.pValue));
set(handles.intensityThresholdEdit, 'String', cell2str(handles.intensityThreshold));
guidata(hObject, handles);
CallBack_slider(hObject, eventdata, -log10(handles.pValue));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% slider bar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_slider(hObject, eventdata, value)

handles = guidata(hObject);

if exist('value')
    s = value; %only display 1st p-value, dawnsong
    if s(1) > get(handles.slider,'max')
        s(1) = get(handles.slider,'max')*0.99;
    end
    if s(1) < get(handles.slider,'min')
        s(1) = get(handles.slider,'min');
    end
else
    s = get(hObject,'Value');
end

set(handles.slider, 'value', s(1));

%multi-thrd, dawnsong, 20120418
for i=1:length(s),
    pvalue(i)=10^(-s(i));
end
%pvalue = 10^(-s(1));

set(handles.pValueEdit,'string', cell2str(num2cell(pvalue)));
%set(handles.pValueEdit,'String',num2str(pvalue));
t = p2t(num2cell(pvalue.*ones(1,length(handles.TF))), handles.df, handles.TF);
handles.intensityThreshold = t;
set(handles.intensityThresholdEdit, 'String', cell2str(handles.intensityThreshold));
handles.pValue = pvalue;
for ii=1:length(handles.TF)
    pos{ii} = find(abs(handles.intensity{ii})>=t{ii});
    handles.currentintensity{ii} = handles.intensity{ii}(pos{ii});
    handles.currentmni{ii} = handles.mni{ii}(pos{ii},:);
end

guidata(hObject,handles);

if get(handles.allIntensityRadio, 'Value')
    CallBack_allIntensityRadio(handles.allIntensityRadio, eventdata);
elseif get(handles.positiveIntensityRadio, 'Value')
    CallBack_allIntensityRadio(handles.positiveIntensityRadio, eventdata, '+');
elseif get(handles.negativeIntensityRadio, 'Value')
    CallBack_allIntensityRadio(handles.negativeIntensityRadio, eventdata, '-');
end

% set(handles.infoTextBox, 'string', {'Don''t drag the slider bar too fast. Release your mouse button at least 1 second later.', ...
%     'This sounds stupid. But there is a bug (probably MatLab bug) which I can'' fix right now.', ...
%     'I suggest you confirm the correctness of the current display by press Enter in the pValue edit box.'}); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% display intensity all+- radios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_allIntensityRadio(hObject, eventdata, pnall)
% pnall = '+', '-', or 'c'. c means current (simply update drawing)
%
handles = guidata(hObject);

currentselect = [];
if get(handles.allIntensityRadio, 'Value'); currentselect = handles.allIntensityRadio; thispnall = 'a'; end
if get(handles.positiveIntensityRadio, 'Value'); currentselect = handles.positiveIntensityRadio; thispnall = '+';  end
if get(handles.negativeIntensityRadio, 'Value'); currentselect = handles.negativeIntensityRadio; thispnall = '-';  end

set(handles.allIntensityRadio, 'Value', 0);
set(handles.positiveIntensityRadio, 'Value', 0);
set(handles.negativeIntensityRadio, 'Value', 0);

if exist('pnall')
    if pnall=='c'
        hObject = currentselect;
        pnall = thispnall;
    end
end

set(hObject, 'Value', 1);

if ~isfield(handles,'currentintensity')
    return
end
for ii=1:length(handles.TF)
	if exist('pnall')
        if pnall == '-'
            pos{ii} = find(handles.currentintensity{ii} < 0);
        elseif pnall == '+'
            pos{ii} = find(handles.currentintensity{ii} > 0);
        elseif pnall == 'a'
            pos{ii} = 1:length(handles.currentintensity{ii});
        end
	else
        pos{ii} = 1:length(handles.currentintensity{ii});
	end
    intensity{ii} = handles.currentintensity{ii}(pos{ii});
    mni{ii} = handles.currentmni{ii}(pos{ii},:);
    cor{ii} = mni2cor(mni{ii}, handles.M{ii});

	if ~isempty(cor{ii})
		A = spm_clusters(cor{ii}');
		pos0 = [];
		for kk = 1:max(A)
            jj = find(A == kk);
            if length(jj) >= handles.clusterSizeThreshold; pos0 = [pos0 jj]; end
		end
		handles.currentDisplayMNI{ii} = mni{ii}(pos0,:);
		handles.currentDisplayIntensity{ii} = intensity{ii}(pos0);
	else
		handles.currentDisplayMNI{ii} = mni{ii}([],:);
		handles.currentDisplayIntensity{ii} = intensity{ii}([]);    
	end
end
[handles.hReg, handles.hSection, handles.hcolorbar] = Draw(handles.currentDisplayMNI, handles.currentDisplayIntensity, hObject, handles);
if get(handles.allIntensityRadio, 'Value') & max(handles.currentDisplayIntensity{1}) < 0
    warndlg('No supra-threshold positive intensity. Only negative intensity is displayed.');
    set(handles.negativeIntensityRadio, 'Value', 1);
    set(handles.allIntensityRadio, 'Value', 0);    
    guidata(hObject, handles);
    return
end
if get(handles.allIntensityRadio, 'Value') & min(handles.currentDisplayIntensity{1}) > 0
    warndlg('No supra-threshold negative intensity. Only positive intensity is displayed.');
    set(handles.positiveIntensityRadio, 'Value', 1);
    set(handles.allIntensityRadio, 'Value', 0);    
    guidata(hObject, handles);
    return
end

try
	set(handles.figure,'currentaxes', handles.glassViewAxes);
	xrange = xlim;
	yrange = ylim;
	try
        delete(handles.hGlassText)
	end
    if ~isempty(handles.selectedCluster,1)
    	%handles.hGlassText = text(xrange(1)+diff(xrange)*0.6, yrange(1)+diff(yrange)*0.9, [num2str(size(handles.selectedCluster,1)) ' clusters selected']);
        set(handles.infoTextBox, 'string', [num2str(size(handles.selectedCluster,1)) ' clusters selected']);
    end
end
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% render view check?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_renderViewCheck(hObject, eventdata)
handles = guidata(hObject);
check = get(hObject, 'Value');
%if check
    CallBack_allIntensityRadio(hObject, eventdata, 'c');
%end

try
if check; % if display render, will make feed invisible
    set(handles.feed,'visible','off');
    set(handles.feedclick,'visible','off');
    set(handles.postfeedclick,'visible','off');
    set(handles.refreshclick,'visible','off');
else
    set(handles.feed,'visible','on');
    set(handles.feedclick,'visible','on');
    set(handles.postfeedclick,'visible','on');
    set(handles.refreshclick,'visible','on');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% change rend style
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_renderStylePop(hObject, eventdata)    
CallBack_renderViewCheck(hObject, eventdata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% which section view target file? list box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_sectionViewListbox(hObject, eventdata)
handles = guidata(hObject);
contents = get(handles.sectionViewListbox,'String');
currentsel = contents{get(handles.sectionViewListbox,'Value')};
handles.sectionViewTargetFile = getSectionViewTargetFile(handles.spmdir, currentsel);
guidata(hObject, handles);
CallBack_allIntensityRadio(hObject, eventdata, 'c');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get sectionviewtargetfile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sectionViewTargetFile = getSectionViewTargetFile(spmdir, selectedcontent)
if findstr('SPM2',spm('ver'))
    fileext = 'mnc';
else%if findstr('SPM5',spm('ver'))
    fileext = 'nii';
end
currentsel = selectedcontent;
if ~isempty(strfind(currentsel, 'single'))
    sectionViewTargetFile = fullfile(spmdir, 'canonical', ['single_subj_T1.' fileext]);
elseif ~isempty(strfind(currentsel, '152PD'))
    sectionViewTargetFile = fullfile(spmdir, 'canonical', ['avg152PD.' fileext]);
elseif ~isempty(strfind(currentsel, '152T1'))
    sectionViewTargetFile = fullfile(spmdir, 'canonical', ['avg152T1.' fileext]);    
elseif ~isempty(strfind(currentsel, '152T2'))
    sectionViewTargetFile = fullfile(spmdir, 'canonical', ['avg152T2.' fileext]);
elseif ~isempty(strfind(currentsel, '305T1'))
    sectionViewTargetFile = fullfile(spmdir, 'canonical', ['avg305T1.' fileext]);
elseif strcmp(currentsel, 'ch2')
    sectionViewTargetFile = fullfile(spmdir, 'canonical', 'ch2.img');    
elseif strcmp(currentsel, 'ch2bet')
    sectionViewTargetFile = fullfile(spmdir, 'canonical', 'ch2bet.img');      
elseif strcmp(currentsel, 'aal')
    sectionViewTargetFile = fullfile(spmdir, 'canonical', 'aal.img');      
elseif strcmp(currentsel, 'brodmann')
    sectionViewTargetFile = fullfile(spmdir, 'canonical', 'brodmann.img');      
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% xhairs in section view?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_xHairCheck(hObject, eventdata)
check = get(hObject, 'Value');
if check
    spm_orthviews('Xhairs','off');
else
    spm_orthviews('Xhairs','on');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% other target file for section view, push button
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_sectionViewMoreTargetPush(hObject, eventdata)    
handles = guidata(hObject);
[filename, pathname, filterindex] = uigetfile('*', 'Pick an target file');
if isequal(filename,0) | isequal(pathname,0)
   return;
end
handles.sectionViewTargetFile = fullfile(pathname, filename);
guidata(hObject, handles);
CallBack_allIntensityRadio(hObject, eventdata, 'c');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set colorbar range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_setTRangeEdit(hObject, eventdata)    
global TMAX_;
handles = guidata(hObject);
TMAX_ = get(hObject, 'String');
if isempty(str2num(TMAX_)) && ~strcmp(TMAX_, 'auto')
    return;
end
guidata(hObject, handles);
CallBack_allIntensityRadio(hObject, eventdata, 'c');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set colorbar range (min)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_setTRangeEdit2(hObject, eventdata)    
global TMIN_;
handles = guidata(hObject);
TMIN_ = get(hObject, 'String');
if isempty(str2num(TMIN_)) && ~strcmp(TMIN_, 'auto')
    return;
end
guidata(hObject, handles);
CallBack_allIntensityRadio(hObject, eventdata, 'c');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% change degree of freedome 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_dfEdit(hObject, eventdata) 
handles = guidata(hObject);
tmp = str2cell(get(hObject, 'String'));
if iscellstr(tmp)
    errordlg('Please input a valid number.','error');
    try
        set(hObject, 'String', handles.df);
    catch
        set(hObject, 'String', '');
    end
    return
end
handles.df=tmp;

if(length(handles.df)==2) % if F test
    handles.df = {[handles.df{1} handles.df{2}]};
end

if isfield(handles,'TF')
    t = p2t(mat2cell(handles.pValue*ones(1,length(handles.TF))), handles.df, handles.TF); 
    handles.intensityThreshold = t;
    set(handles.intensityThresholdEdit, 'String', cell2str(t));    
end

guidata(hObject, handles);
CallBack_pValueEdit(handles.pValueEdit, eventdata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get structure push
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_getStructurePush(hObject, eventdata) 
handles = guidata(hObject);
xyz = spm_XYZreg('GetCoords',handles.hReg);
tmp_coor = cuixuFindStructure(xyz', handles.DB);
set(handles.structureEdit,'String', tmp_coor{1});
handles.currentxyz = xyz';
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% structure edit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_structureEdit(action, xyz, thisfun, hReg)

try
    handles=guidata(hReg);
catch
    return;
end

handles.currentxyz = xyz';
handles.hReg = hReg;
guidata(hReg, handles);

try
    [tmp_coor, cellstructure] = cuixuFindStructure(xyz', handles.DB);
catch
    return;
end

set(handles.structureEdit,'String', tmp_coor{1});

for ii=[5 3 2 6 1 4]
    if strfind('undefined', cellstructure{ii})
        continue;
    else
        set(handles.searchContentEdit,'string', trimStructureStr(cellstructure{ii}));
        return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% trim str
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = trimStructureStr(str)
pos = findstr('(', str);
if ~isempty(pos)
    str(pos-1:end)=[];
end
out = str;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cluster size threshold edit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_clusterSizeThresholdEdit(hObject, eventdata) 
handles = guidata(hObject);
tmp = str2double(get(hObject, 'String'));
if isnan(tmp) | tmp<0
    errordlg('Please input a valid number.','error');
    try
        set(hObject, 'String', handles.clusterSizeThreshold);
    catch
        set(hObject, 'String', '5');
    end
    return
end
handles.clusterSizeThreshold=tmp;

guidata(hObject, handles);
CallBack_allIntensityRadio(hObject, eventdata, 'c');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% save image push button
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_saveImagePush(hObject, eventdata, thisfilename, isMask) 
handles = guidata(hObject);

if exist('thisfilename')
    if ~strcmp(deblank(thisfilename), '')
        if isfield(handles,'imageFileName')
            if ~isempty(handles.imageFileName)
                mni2mask(cell2mat(handles.currentDisplayMNI'), thisfilename, cell2mat(handles.currentDisplayIntensity), handles.M{1}, handles.DIM{1}, handles.imageFileName{1});
                return;
            end
        end
        mni2mask(cell2mat(handles.currentDisplayMNI'), thisfilename, cell2mat(handles.currentDisplayIntensity), handles.M{1}, handles.DIM{1});
        return;
    end
end

[filename, pathname] = uiputfile('*.img', 'Save image file as', get(handles.saveImageFileEdit, 'string'));
if isequal(filename,0) | isequal(pathname,0)
   return
else
   thisfilename = fullfile(pathname, filename);
end

if isfield(handles,'imageFileName')
    if ~isempty(handles.imageFileName{1})
        mni2mask(cell2mat(handles.currentDisplayMNI'), thisfilename, cell2mat(handles.currentDisplayIntensity'), handles.M{1}, handles.DIM{1}, handles.imageFileName{1}, isMask);
        return;
    end
end

mni2mask(cell2mat(handles.currentDisplayMNI'), thisfilename, cell2mat(handles.currentDisplayIntensity), handles.M{1}, handles.DIM{1}, '', isMask);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% save image edit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_saveImageFileEdit(hObject, eventdata) 
handles = guidata(hObject);
CallBack_saveImagePush(handles.saveImagePush, eventdata, get(hObject,'string'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% save result push
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_saveResultPSPush(hObject, eventdata, thisfilename) 

handles = guidata(hObject);

if exist('thisfilename')
    if ~strcmp(deblank(thisfilename), '')
        spm_print(handles.figure);
        if strcmp(handles.system, 'linux')
            system(['ps2pdf spm2.ps']);
            system(['mv spm2.ps ' thisfilename]);
            [p,f,ext]=fileparts(thisfilename);
            system(['mv spm2.pdf ' fullfile(p,f) '.pdf']);            
        elseif strcmp(handles.system, 'windows')
            system(['move spm2.ps ' '"' thisfilename '"']);
		end        
        return;
    end
end

[filename, pathname] = uiputfile('*.ps', 'Save result as', get(handles.saveResultPSEdit, 'string'));
if isequal(filename,0) | isequal(pathname,0)
   return
else
   thisfilename = fullfile(pathname, filename);
end


% print
H  = findobj(get(handles.figure,'Children'),'flat','Type','axes');
un = cellstr(get(H,'Units'));
pos = get(H,'position');
index = [];

for ii=1:length(H)
    if findstr('pixels', un{ii})
        continue;
    end
    if pos{ii}(1)>0.4 & pos{ii}(2) > 0.5
        set(H(ii),'position',[pos{ii}(1), pos{ii}(2), pos{ii}(3)*3/4, pos{ii}(4)]);
        index = [index ii];
    end
end


spm_print(handles.figure);
if strcmp(handles.system, 'linux')
    system(['ps2pdf spm2.ps']);
    system(['mv spm2.ps ' thisfilename]);
    [p,f,ext]=fileparts(thisfilename);
    system(['mv spm2.pdf ' fullfile(p,f) '.pdf']);            
elseif strcmp(handles.system, 'windows')
    system(['move spm2.ps ' '"' thisfilename '"']);
end 

% set the position back
set(H(index), {'position'}, pos(index)); 


% printstr = ['print -dpsc2 -painters -noui ' '''' thisfilename ''''];
% try
%     orient portrait
%     eval(printstr);
%     printsuccess = 1;
% catch
%     errordlg('Print to ps file failed', 'print error');
%     printsuccess = 0;
% end
%set(H,{'Units'},un);
% if strcmp(handles.system, 'linux') & printsuccess == 1
%     system(['ps2pdf ' '''' thisfilename '''']);
% end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% save result edit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_saveResultPSEdit(hObject, eventdata) 
handles = guidata(hObject);
CallBack_saveResultPSPush(hObject, eventdata, get(handles.saveResultPSEdit,'string'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% select a cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_selectThisClusterPush(hObject, eventdata) 
handles = guidata(hObject);

try
    xyz = handles.currentxyz';
catch
    xyz = spm_XYZreg('GetCoords',handles.hReg);
end
try
    handles.selectedCluster = [handles.selectedCluster; xyz'];
catch
    handles.selectedCluster = xyz';
end

set(handles.figure,'currentaxes', handles.glassViewAxes);
xrange = xlim;
yrange = ylim;
try
    delete(handles.hGlassText)
end
%handles.hGlassText = text(xrange(1)+diff(xrange)*0.6, yrange(1)+diff(yrange)*0.9, [num2str(size(handles.selectedCluster,1)) ' clusters selected']);
set(handles.infoTextBox, 'string', [num2str(size(handles.selectedCluster,1)) ' clusters selected']);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% unselect a cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_clearSelectedClusterPush(hObject, eventdata) 
handles = guidata(hObject);
try
    handles = rmfield(handles,'selectedCluster');
end
set(handles.figure,'currentaxes', handles.glassViewAxes);
xrange = xlim;
yrange = ylim;
try
    %delete(handles.hGlassText)
    set(handles.infoTextBox, 'string', ['No clusters selected']);
end

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pick a cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_pickThisClusterPush(hObject, eventdata) 
handles = guidata(hObject);
mni = cell2mat(handles.currentDisplayMNI');
if isempty(mni)
    %errordlg('No cluster is picked up.','oops');
    set(handles.infoTextBox, 'string', 'No cluster is picked up.'); 
    beep
    return;
end

if ~isfield(handles, 'selectedCluster') | isempty(handles.selectedCluster)
    try
        xyz = handles.currentxyz';
    catch
        xyz = spm_XYZreg('GetCoords',handles.hReg);
    end
    handles.selectedCluster = xyz';
end

intensity = cell2mat(handles.currentDisplayIntensity');
cor = mni2cor(mni, handles.M{1});
A = spm_clusters(cor');
xyzcor = mni2cor(handles.selectedCluster, handles.M{1});

pos = [];
for ii = 1:size(xyzcor,1)
    pos0 = find(cor(:,1)==xyzcor(ii,1) & cor(:,2)==xyzcor(ii,2) & cor(:,3)==xyzcor(ii,3));
    if isempty(pos0)
        continue;
    end
    pos = [pos find(A==A(pos0(1)))];
end
if isempty(pos)
    %errordlg('No cluster is picked up.','oops');
    set(handles.infoTextBox, 'string', 'No cluster is picked up.'); 
    beep    
    return
end

pos = unique(pos);

tmpmni = mni(pos,:);
tmpintensity = intensity(pos);

[B,I,J] = unique(tmpmni, 'rows');
handles.currentDisplayMNI = {B};
handles.currentDisplayIntensity = {tmpintensity(I,:)};

handles.currentxyz = handles.selectedCluster(end,:);

[handles.hReg, handles.hSection, handles.hcolorbar] = Draw(handles.currentDisplayMNI, handles.currentDisplayIntensity, hObject, handles);
set(handles.figure,'currentaxes', handles.glassViewAxes);
xrange = xlim;
yrange = ylim;

handles.selectedCluster = [];
try
    delete(handles.hGlassText)
end

guidata(hObject, handles);

set(handles.thisClusterSizeEdit,'string', num2str(size(handles.currentDisplayMNI{1},1)));
str = get(handles.searchContentEdit, 'string');
pos = findstr(' ', str);
str(pos) = [];
files = [];
for ii=1:length(handles.imageFileName)
    [a,b,c] = fileparts(handles.imageFileName{ii});
    files = [files b];
end
if ~isempty(files)
    files = ['_from_' files];
end

set(handles.saveImageFileEdit, 'string', [str files '.img']);
        
% list structure of voxels in this cluster
[a, b] = cuixuFindStructure(cell2mat(handles.currentDisplayMNI'), handles.DB);
names = unique(b(:));
index = NaN*zeros(length(b(:)),1);
for ii=1:length(names)
    pos = find(strcmp(b(:),names{ii}));
    index(pos) = ii;
end

for ii=1:max(index)
    report{ii,1} = names{ii};
    report{ii,2} = length(find(index==ii));
end
for ii=1:size(report,1)
    for jj=ii+1:size(report,1)
        if report{ii,2} < report{jj,2}
            tmp = report(ii,:);
            report(ii,:) = report(jj,:);
            report(jj,:) = tmp;
        end
    end
end
report = [{'structure','# voxels'}; {'--TOTAL # VOXELS--', length(a)}; report];
% format long
% disp(b)
% disp(report)
% format
report2 = {sprintf('%s\t%s',report{1,2}, report{1,1}),''};
for ii=2:size(report,1)
    if strcmp('undefined', report{ii,1}); continue; end
    report2 = [report2, {sprintf('%5d\t%s',report{ii,2}, report{ii,1})}];
end

report2 = [report2, {'','select and Ctrl-C to copy'}];
% f = figure('unit','normalized', 'menubar','none', 'position', [0.3 0.2 0.2 min(0.7,0.016*length(report2))], 'name', 'cluster information', 'NumberTitle','off');
% hEdit = uicontrol(f,'style','edit',...
%         'unit','normalized','position',[0 0 1 1],...
%         'horizontal','left',...
%         'BackgroundColor', 'w',...
%         'String', report2,...
%         'max',2,'min',0);

set(handles.infoTextBox, 'string', report2);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_reportPush(hObject, eventdata) 
handles = guidata(hObject);

set(handles.infoTextBox, 'string', 'This function will generate a report of this functional image. Check out your MatLab command window. You may copy and past to a file for later use.');
if(length(handles.imageFileName)~=1)
    set(handles.infoTextBox, 'string', 'This function only applies to 1 functional image.');
    return;
end

mni = cell2mat(handles.currentDisplayMNI');
if isempty(mni)
    %errordlg('No cluster is picked up.','oops');
    set(handles.infoTextBox, 'string', 'No cluster is found. So no report will be generated.'); 
    beep
    return;
end

disp(handles.imageFileName{1});
disp(['Type: ' num2str(handles.TF{1})]);
disp(['df: ' num2str(handles.df{1})]);
disp('Threshold')
disp(['-- p value = ' num2str(handles.pValue)]);
disp(['-- intensity = ' num2str(handles.intensityThreshold{1})]);
disp(['-- cluster size = ' num2str(handles.clusterSizeThreshold)]);

intensity = cell2mat(handles.currentDisplayIntensity');
cor = mni2cor(mni, handles.M{1});
A = spm_clusters(cor');

clusterID = unique(A);
numClusters = length(clusterID);
disp(['Number of clusters found: ' num2str(numClusters)]);

for mm = clusterID
    pos = find(A == clusterID(mm));
    numVoxels = length(pos);
    tmpmni = mni(pos,:);
    tmpintensity = intensity(pos);
    
    peakpos = find(abs(tmpintensity) == max(abs(tmpintensity)));
    peakpos = peakpos(1);
    peakcoord = tmpmni(peakpos,:);
    peakintensity = tmpintensity(peakpos);

    % list structure of voxels in this cluster
    [a, b] = cuixuFindStructure(tmpmni, handles.DB);
    names = unique(b(:));
    index = NaN*zeros(length(b(:)),1);
    for ii=1:length(names)
        pos = find(strcmp(b(:),names{ii}));
        index(pos) = ii;
    end

    report = {};
    
    for ii=1:max(index)
        report{ii,1} = names{ii};
        report{ii,2} = length(find(index==ii));
    end
    for ii=1:size(report,1)
        for jj=ii+1:size(report,1)
            if report{ii,2} < report{jj,2}
                tmp = report(ii,:);
                report(ii,:) = report(jj,:);
                report(jj,:) = tmp;
            end
        end
    end
    report = [{'structure','# voxels'}; {'--TOTAL # VOXELS--', length(a)}; report];

    report2 = {sprintf('%s\t%s',report{1,2}, report{1,1}),''};
    for ii=2:size(report,1)
        if strcmp('undefined', report{ii,1}); continue; end
        report2 = [report2, {sprintf('%5d\t%s',report{ii,2}, report{ii,1})}];
    end

    disp(['----------------------'])
    disp(['Cluster ' num2str(mm)])
    disp(['Number of voxels: ' num2str(numVoxels)])
    disp(['Peak MNI coordinate: ' num2str(peakcoord)])
    [a,b] = cuixuFindStructure(peakcoord, handles.DB);
    disp(['Peak MNI coordinate region: ' a{1}]);
    disp(['Peak intensity: ' num2str(peakintensity)])
    for kk=1:length(report2)
        disp(report2{kk});
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sliceView push
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_sliceViewPush(hObject, eventdata) 
if(nargin == 1)
    eventdata = [];
end
handles = guidata(hObject);


global sliceview

if(~isfield(sliceview, 'fig') || ~ishandle(sliceview.fig))
    sliceview.viewtype = 't';    
    sliceview.row = 8;
    sliceview.col = 8;
    sliceview.spacing = 4;
    sliceview.page = 1;
    sliceview.data = {{},{},{}}; % t,s,c
    sliceview.slices = {[],[],[]};% t,s,c
    sliceview.colormap = '';
    sliceview.fig = figure('color','k', 'unit','normalized','position',[0.1 0.1 .6 .8],'toolbar','none', 'name', 'xjView slice view', 'NumberTitle','off');
    sliceview.ax = axes('Visible','on','DrawMode','fast','Parent',sliceview.fig,...
    'YDir','normal','Ydir','normal','XTick',[],'YTick',[], 'position', [0.15 0.05 .8 .9]);
    %handles.sliceview.d  = image([],'Tag','Transverse','Parent',handles.sliceview.ax);
    set(sliceview.ax,'XTick',[],'YTick',[]);
    axis equal
    set(sliceview.ax,'color','k');
    %setcolormap(colormp)
    width = 0.05;
    height = 0.025;
    step = 0.025;
    labeloffset = step/2;
    
     uicontrol(sliceview.fig,'style','text',...
        'unit','normalized','position',[0 step-labeloffset width height],...
        'horizontal','left',...
        'BackGroundColor', 'k', ...
        'ForeGroundColor', [1 1 1], ...
        'String', 'row #');
    uicontrol(sliceview.fig,'style','text',...
        'unit','normalized','position',[0 step*3-labeloffset width height],...
        'horizontal','left',...
        'BackGroundColor', 'k', ...
        'ForeGroundColor', [1 1 1], ...
        'String', 'column #');
    uicontrol(sliceview.fig,'style','text',...
        'unit','normalized','position',[0 step*5-labeloffset width*2 height],...
        'horizontal','left',...
        'BackGroundColor', 'k', ...
        'ForeGroundColor', [1 1 1], ...
        'String', 'spacing (mm)');
    uicontrol(sliceview.fig,'style','text',...
        'unit','normalized','position',[0 step*7-labeloffset width height],...
        'horizontal','left',...
        'BackGroundColor', 'k', ...
        'ForeGroundColor', [1 1 1], ...
        'String', 'scroll');
    
    uicontrol(sliceview.fig,'style','edit',...
        'unit','normalized','position',[0 0 width height],...
        'horizontal','left',...
        'BackgroundColor', 'w',...
        'String', sliceview.row,...
        'callback',@sliceViewRowChanged);
   
    uicontrol(sliceview.fig,'style','edit',...
        'unit','normalized','position',[0 step*2 width height],...
        'horizontal','left',...
        'BackgroundColor', 'w',...
        'String', sliceview.col,...
        'callback',@sliceViewColChanged);
    
    uicontrol(sliceview.fig,'style','edit',...
        'unit','normalized','position',[0 step*4 width height],...
        'horizontal','left',...
        'BackgroundColor', 'w',...
        'String', sliceview.spacing,...
        'callback',@sliceViewSpacingChanged);
    
%     uicontrol(sliceview.fig,'style','edit',...
%         'unit','normalized','position',[0 step*6 width height],...
%         'horizontal','left',...
%         'BackgroundColor', 'w',...
%         'String', sliceview.page,...
%         'callback',@sliceViewPageChanged);
    uicontrol(sliceview.fig,'style','push',...
        'unit','normalized','position',[0 step*6 width height],...
        'horizontal','left',...
        'BackgroundColor', 'k',...
        'ForegroundColor', 'w', ...
        'String', 'up',...
        'callback',@sliceViewPageChangedP);
    uicontrol(sliceview.fig,'style','push',...
        'unit','normalized','position',[0.05 step*6 width height],...
        'horizontal','left',...
        'BackgroundColor', 'k',...
        'ForegroundColor', 'w', ...
        'String', 'down',...
        'callback',@sliceViewPageChangedN);
    uicontrol(...
		'Units','normalized', ...
		'ListboxTop',0, ...
		'Position',[0 step*8 width*2 height],...
		'String',{'transverse';'coronal';'sagittal'}, ...
		'Style','popupmenu', ...
		'value',1,...
        'callback',@sliceViewTypeChanged);
    
end

figure(sliceview.fig)

viewtype =    sliceview.viewtype;   
row =         sliceview.row;
col =         sliceview.col;
spacing =     sliceview.spacing;
page =        sliceview.page;
slice_fig = sliceview.fig;
ax = sliceview.ax;
%d = handles.sliceview.d;

[slicedata, colormp, slices] = cuixu_getSliceViewData(viewtype,row,col, spacing, page);

if isempty(slices)
    return;
end

for ii=1:length(slices)
    if(viewtype == 's')
        postmp = find(slices(ii) - sliceview.slices{2} == 0);
        if(isempty(postmp))
            sliceview.data{2}{end+1} = slicedata{ii};
            sliceview.slices{2}(end+1) = slices(ii);
        end   
    elseif(viewtype == 't')
        postmp = find(slices(ii) - sliceview.slices{1} == 0);
        if(isempty(postmp))
            sliceview.data{1}{end+1} = slicedata{ii};
            sliceview.slices{1}(end+1) = slices(ii);
        end    
    elseif(viewtype == 'c')
        postmp = find(slices(ii) - sliceview.slices{3} == 0);
        if(isempty(postmp))
            sliceview.data{3}{end+1} = slicedata{ii};
            sliceview.slices{3}(end+1) = slices(ii);
        end   
    end  
    
end



%slice_fig = figure('color','k', 'unit','normalized','position',[0.1 0.1 .6 .8],'toolbar','none');


if(length(size(slicedata{1})) == 3)
    [nx, ny, nz] = size(slicedata{1});
    slicedatafinal = zeros(nx*row, ny*col, nz );
    for ii=1:length(slicedata)
        slicedatafinal(nx*(floor((ii-1)/col))+1:nx*(1+floor((ii-1)/col)), ny*(mod(ii-1,col))+1:ny*(mod(ii-1,col)+1), :) = slicedata{ii};
    end
else
    [nx, ny] = size(slicedata{1});
    slicedatafinal = zeros(nx*row, ny*col );
    for ii=1:length(slicedata)
        slicedatafinal(nx*(floor((ii-1)/col))+1:nx*(1+floor((ii-1)/col)), ny*(mod(ii-1,col))+1:ny*(mod(ii-1,col)+1)) = slicedata{ii};
    end
end


try
    delete(handles.sliceview.d)
catch
    [];
end

handles.sliceview.d  = image(slicedatafinal,'Tag','Transverse','Parent',sliceview.ax);

% put slice positions
for ii=1:length(slicedata)
    %text(nx*(floor((ii-1)/col))+1:nx*(1+floor((ii-1)/col)), ny*(mod(ii-1,col))+1:ny*(mod(ii-1,col)+1), num2str(slices(ii)), 'color', 'w');
    text(ny*(mod(ii-1,col))+1, nx*(floor((ii-1)/col))+1+20, num2str(slices(ii)), 'color', 'w');
end

set(sliceview.ax,'XTick',[],'YTick',[]);
axis(sliceview.ax, 'equal');
set(sliceview.ax,'color','k');
%setcolormap(colormp)
% for ii=1:length(slicedata)
%     
%     ax = subplot(row,col,ii);
% %	ax = axes('Visible','on','DrawMode','fast','Parent',slice_fig,...
% %		'YDir','normal','Ydir','normal','XTick',[],'YTick',[]);
%     set(ax,'Visible','on','DrawMode','fast','Parent',slice_fig,...
% 		'YDir','normal','Ydir','normal','XTick',[],'YTick',[], 'color','k');
%     d  = image(slicedata{ii},'Tag','Transverse','Parent',ax);
%     set(ax,'XTick',[],'YTick',[]);
%     title(pos{ii},'color','w')
%     axis equal
%     set(ax,'color','k');
% 
%     %set(d,'Cdata',imgt);
%     if mn*mx < 0
%         setcolormap('gray-hot-cold')                    
%     elseif mx > 0
%         setcolormap('gray-hot');
%     else
%         setcolormap('gray-cold')
%     end  
% end

guidata(hObject, handles);

function sliceViewRowChanged(hObject, eventdata) 
global sliceview
sliceview.row = str2num(get(hObject,'string'));
CallBack_sliceViewPush(hObject, eventdata); 

function sliceViewColChanged(hObject, eventdata) 
global sliceview
sliceview.col = str2num(get(hObject,'string'));
CallBack_sliceViewPush(hObject, eventdata); 

function sliceViewSpacingChanged(hObject, eventdata) 
global sliceview
sliceview.spacing = str2num(get(hObject,'string'));
CallBack_sliceViewPush(hObject, eventdata); 

function sliceViewPageChangedN(hObject, eventdata) 
global sliceview
sliceview.page = sliceview.page + 1;
CallBack_sliceViewPush(hObject, eventdata); 

function sliceViewPageChangedP(hObject, eventdata) 
global sliceview
sliceview.page = sliceview.page - 1;
CallBack_sliceViewPush(hObject, eventdata); 

function sliceViewTypeChanged(hObject, eventdata) 
global sliceview
type = get(hObject,'value');
if(type == 1)
    sliceview.viewtype = 't';
elseif(type == 2)    
    sliceview.viewtype = 'c';
else
    sliceview.viewtype = 's';
end
%sliceview.data = {{},{},{}};
%sliceview.slices = {[],[],[]};
sliceview.page = 0;
CallBack_sliceViewPush(hObject, eventdata); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get slice view data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [slicedata, colormp, slices] = cuixu_getSliceViewData(viewtype,row,col, spacing,page)

global sliceview


slicedata = {};
pos = {};
slices = -80+(page-1)*spacing*col:spacing:100+(page)*spacing;

if(isempty(slices))
    colormp = '';
    return
end

global st
bb   = st.bb;
Dims = round(diff(bb)'+1);
is   = inv(st.Space);
cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);

kk = 1;
for slice = slices
    if(kk > row*col)
        slices = slices(1:kk-1);
        break;
    end
    
if(viewtype == 's')
    cent(2) = -slice;
    postmp = find(slice - sliceview.slices{2} == 0);
    if(~isempty(postmp))
        slicedata{kk} = sliceview.data{2}{postmp(1)};
        kk = kk+1;
        continue;
    end
elseif(viewtype == 't')
    cent(3) = slice;
    postmp = find(slice - sliceview.slices{1} == 0);
    if(~isempty(postmp))
        slicedata{kk} = sliceview.data{1}{postmp(1)};
        kk = kk+1;
        continue;
    end
elseif(viewtype == 'c')
    cent(1) = slice;
    postmp = find(slice - sliceview.slices{3} == 0);
    if(~isempty(postmp))
        slicedata{kk} = sliceview.data{3}{postmp(1)};
        kk = kk+1;
        continue;
    end
end       

if(viewtype == 's')

for i = 1
	M = st.vols{i}.premul*st.vols{i}.mat;


	CM0 = [	1 0 0 -bb(1,1)+1
		0 0 1 -bb(1,3)+1
		0 1 0 -cent(2)
		0 0 0 1];
	CM = inv(CM0*(st.Space\M));
	CD = Dims([1 3]);



	ok=1;
   
    eval('imgc  = (spm_slice_vol(st.vols{i},CM,CD,st.hld))'';','ok=0;');

	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.vols{i}.fname);
	else,
                % get min/max threshold
                if strcmp(st.vols{i}.window,'auto')
                        mn = -Inf;
                        mx = Inf;
                else
                        mn = min(st.vols{i}.window);
                        mx = max(st.vols{i}.window);
                end;
                % threshold images
                imgc = max(imgc,mn); imgc = min(imgc,mx);
                % compute intensity mapping, if histeq is available
                if license('test','image_toolbox') == 0
                    st.vols{i}.mapping = 'linear';
                end;
                switch st.vols{i}.mapping,
                 case 'linear',
                 case 'histeq',
                  % scale images to a range between 0 and 1
                  imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                  img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                  imgc = reshape(img(numel(imgt1)+[1:numel(imgc1)]),size(imgc1));
                  mn = 0;
                  mx = 1;
                 case 'quadhisteq',
                  % scale images to a range between 0 and 1
                  imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                  imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                  imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                  img  = histeq([imgt1(:).^2; imgc1(:).^2; imgs1(:).^2],1024);
                  imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                  imgc = reshape(img(numel(imgt1)+[1:numel(imgc1)]),size(imgc1));
                  imgs = reshape(img(numel(imgt1)+numel(imgc1)+[1:numel(imgs1)]),size(imgs1));
                  mn = 0;
                  mx = 1;
                 case 'loghisteq',
                  warning off % messy - but it may avoid extra queries
                  imgt = log(imgt-min(imgt(:)));
                  imgc = log(imgc-min(imgc(:)));
                  imgs = log(imgs-min(imgs(:)));
                  warning on
                  imgt(~isfinite(imgt)) = 0;
                  imgc(~isfinite(imgc)) = 0;
                  imgs(~isfinite(imgs)) = 0;
                  % scale log images to a range between 0 and 1
                  imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                  imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                  imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                  img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                  imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                  imgc = reshape(img(numel(imgt1)+[1:numel(imgc1)]),size(imgc1));
                  imgs = reshape(img(numel(imgt1)+numel(imgc1)+[1:numel(imgs1)]),size(imgs1));
                  mn = 0;
                  mx = 1;
                end;
                % recompute min/max for display
                if strcmp(st.vols{i}.window,'auto')
                    mx = -inf; mn = inf;
                end;
                
                if ~isempty(imgc),
			tmp = imgc(isfinite(imgc));
                        mx = max([mx max(max(tmp))]);
                        mn = min([mn min(min(tmp))]);
                end;
                
                if mx==mn, mx=mn+eps; end;

		if isfield(st.vols{i},'blobs'),
			if ~isfield(st.vols{i}.blobs{1},'colour'),
				% Add blobs for display using the split colourmap
				scal = 64/(mx-mn);
				dcoff = -mn*scal;
				imgc = imgc*scal+dcoff;

				if isfield(st.vols{i}.blobs{1},'max'),
					mx = st.vols{i}.blobs{1}.max;
				else,
					mx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
					st.vols{i}.blobs{1}.max = mx;
				end;
				if isfield(st.vols{i}.blobs{1},'min'),
					mn = st.vols{i}.blobs{1}.min;
				else,
					mn = min([0 minval(st.vols{i}.blobs{1}.vol)]);
					st.vols{i}.blobs{1}.min = mn;
				end;

				vol  = st.vols{i}.blobs{1}.vol;
				M    = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;

                tmpc = spm_slice_vol(vol,inv(CM0*(st.Space\M)),CD,[0 NaN])';

				

				%tmpt_z = find(tmpt==0);tmpt(tmpt_z) = NaN;
				%tmpc_z = find(tmpc==0);tmpc(tmpc_z) = NaN;
				%tmps_z = find(tmps==0);tmps(tmps_z) = NaN;

				sc   = 64/(mx-mn);
				off  = 65.51-mn*sc;
				msk  = find(isfinite(tmpc)); imgc(msk) = off+tmpc(msk)*sc;

				cmap = get(st.fig,'Colormap');

                %figure(st.fig)
                if mn*mx < 0
                    setcolormap('gray-hot-cold')                    
                elseif mx > 0
                    setcolormap('gray-hot');
                else
                    setcolormap('gray-cold')
                end                
                %                redraw_colourbar(i,1,[mn mx],[1:64]'+64); 
			elseif isstruct(st.vols{i}.blobs{1}.colour),
				% Add blobs for display using a defined
                                % colourmap

				% colourmaps
				gryc = [0:63]'*ones(1,3)/63;
				actc = ...
				    st.vols{1}.blobs{1}.colour.cmap;
				actp = ...
				    st.vols{1}.blobs{1}.colour.prop;
				
				% scale grayscale image, not finite -> black
				imgc = scaletocmap(imgc,mn,mx,gryc,65);
				gryc = [gryc; 0 0 0];
				
				% get max for blob image
				vol = st.vols{i}.blobs{1}.vol;
				mat = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				if isfield(st.vols{i}.blobs{1},'max'),
					cmx = st.vols{i}.blobs{1}.max;
				else,
					cmx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
				end;
				if isfield(st.vols{i}.blobs{1},'min'),
					cmn = st.vols{i}.blobs{1}.min;
				else,
					cmn = -cmx;
				end;

				% get blob data
				vol  = st.vols{i}.blobs{1}.vol;
				M    = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				tmpc = spm_slice_vol(vol,inv(CM0*(st.Space\M)),CD,[0 NaN])';
				
				% actimg scaled round 0, black NaNs
				topc = size(actc,1)+1;
				tmpc = scaletocmap(tmpc,cmn,cmx,actc,topc);
				actc = [actc; 0 0 0];
				
				% combine gray and blob data to
				% truecolour
				imgc = reshape(actc(tmpc(:),:)*actp+ ...
					       gryc(imgc(:),:)*(1-actp), ...
					       [size(imgc) 3]);
				
				
			else,
				% Add full colour blobs - several sets at once
				scal  = 1/(mx-mn);
				dcoff = -mn*scal;

				wc = zeros(size(imgc));

				imgc  = repmat(imgc*scal+dcoff,[1,1,3]);

				cimgc = zeros(size(imgc));

				for j=1:length(st.vols{i}.blobs), % get colours of all images first
					if isfield(st.vols{i}.blobs{j},'colour'),
						colour(j,:) = reshape(st.vols{i}.blobs{j}.colour, [1 3]);
					else,
						colour(j,:) = [1 0 0];
					end;
				end;
				%colour = colour/max(sum(colour));

				for j=1:length(st.vols{i}.blobs),
					if isfield(st.vols{i}.blobs{j},'max'),
						mx = st.vols{i}.blobs{j}.max;
					else,
						mx = max([eps max(st.vols{i}.blobs{j}.vol(:))]);
						st.vols{i}.blobs{j}.max = mx;
					end;
					if isfield(st.vols{i}.blobs{j},'min'),
						mn = st.vols{i}.blobs{j}.min;
					else,
						mn = min([0 min(st.vols{i}.blobs{j}.vol(:))]);
						st.vols{i}.blobs{j}.min = mn;
					end;

					vol  = st.vols{i}.blobs{j}.vol;
					M    = st.Space\st.vols{i}.premul*st.vols{i}.blobs{j}.mat;
                    tmpc = spm_slice_vol(vol,inv(CM0*M),CD,[0 NaN])';
                    % check min/max of sampled image
                    % against mn/mx as given in st
                    tmpc(tmpc(:)<mn) = mn;
                    tmpc(tmpc(:)>mx) = mx;
					tmpc = (tmpc-mn)/(mx-mn);
					tmpc(~isfinite(tmpc)) = 0;

					cimgc = cimgc + cat(3,tmpc*colour(j,1),tmpc*colour(j,2),tmpc*colour(j,3));

					wc = wc + tmpc;
                                        cdata=permute(shiftdim([1/64:1/64:1]'* ...
                                                               colour(j,:),-1),[2 1 3]);
                                        redraw_colourbar(i,j,[mn mx],cdata);
				end;

				imgc = repmat(1-wc,[1 1 3]).*imgc+cimgc;

				imgc(imgc<0)=0; imgc(imgc>1)=1;
			end;
		else,
			scal = 64/(mx-mn);
			dcoff = -mn*scal;
			imgc = imgc*scal+dcoff;
		end;

% 		set(st.vols{i}.ax{1}.d,'HitTest','off', 'Cdata',imgt);
% 		set(st.vols{i}.ax{1}.lx,'HitTest','off',...
% 			'Xdata',[0 TD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
% 		set(st.vols{i}.ax{1}.ly,'HitTest','off',...
% 			'Ydata',[0 TD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));
% 
% 		set(st.vols{i}.ax{2}.d,'HitTest','off', 'Cdata',imgc);
% 		set(st.vols{i}.ax{2}.lx,'HitTest','off',...
% 			'Xdata',[0 CD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
% 		set(st.vols{i}.ax{2}.ly,'HitTest','off',...
% 			'Ydata',[0 CD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));
% 
% 		set(st.vols{i}.ax{3}.d,'HitTest','off','Cdata',imgs);
% 		if st.mode ==0,
% 			set(st.vols{i}.ax{3}.lx,'HitTest','off',...
% 				'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
% 			set(st.vols{i}.ax{3}.ly,'HitTest','off',...
% 				'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(cent(3)-bb(1,3)+1));
% 		else,
% 			set(st.vols{i}.ax{3}.lx,'HitTest','off',...
% 				'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
% 			set(st.vols{i}.ax{3}.ly,'HitTest','off',...
% 				'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(bb(2,2)+1-cent(2)));
% 		end;


	end;
end;

elseif(viewtype == 't')

for i = 1
	M = st.vols{i}.premul*st.vols{i}.mat;
	TM0 = [	1 0 0 -bb(1,1)+1
		0 1 0 -bb(1,2)+1
		0 0 1 -cent(3)
		0 0 0 1];
	TM = inv(TM0*(st.Space\M));
	TD = Dims([1 2]);



	ok=1;
    eval('imgt  = (spm_slice_vol(st.vols{i},TM,TD,st.hld))'';','ok=0;');

	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.vols{i}.fname);
	else,
                % get min/max threshold
                if strcmp(st.vols{i}.window,'auto')
                        mn = -Inf;
                        mx = Inf;
                else
                        mn = min(st.vols{i}.window);
                        mx = max(st.vols{i}.window);
                end;
                % threshold images
                imgt = max(imgt,mn); imgt = min(imgt,mx);
                % compute intensity mapping, if histeq is available
                if license('test','image_toolbox') == 0
                    st.vols{i}.mapping = 'linear';
                end;
                switch st.vols{i}.mapping,
                 case 'linear',
                 case 'histeq',
                  % scale images to a range between 0 and 1
                  imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                  img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                  imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                  mn = 0;
                  mx = 1;
                 case 'quadhisteq',
                  % scale images to a range between 0 and 1
                  imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                  imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                  imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                  img  = histeq([imgt1(:).^2; imgc1(:).^2; imgs1(:).^2],1024);
                  imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                  imgc = reshape(img(numel(imgt1)+[1:numel(imgc1)]),size(imgc1));
                  imgs = reshape(img(numel(imgt1)+numel(imgc1)+[1:numel(imgs1)]),size(imgs1));
                  mn = 0;
                  mx = 1;
                 case 'loghisteq',
                  warning off % messy - but it may avoid extra queries
                  imgt = log(imgt-min(imgt(:)));
                  imgc = log(imgc-min(imgc(:)));
                  imgs = log(imgs-min(imgs(:)));
                  warning on
                  imgt(~isfinite(imgt)) = 0;
                  imgc(~isfinite(imgc)) = 0;
                  imgs(~isfinite(imgs)) = 0;
                  % scale log images to a range between 0 and 1
                  imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                  imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                  imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                  img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                  imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                  imgc = reshape(img(numel(imgt1)+[1:numel(imgc1)]),size(imgc1));
                  imgs = reshape(img(numel(imgt1)+numel(imgc1)+[1:numel(imgs1)]),size(imgs1));
                  mn = 0;
                  mx = 1;
                end;
                % recompute min/max for display
                if strcmp(st.vols{i}.window,'auto')
                    mx = -inf; mn = inf;
                end;
                if ~isempty(imgt),
			tmp = imgt(isfinite(imgt));
                        mx = max([mx max(max(tmp))]);
                        mn = min([mn min(min(tmp))]);
                end;
                
                if mx==mn, mx=mn+eps; end;

		if isfield(st.vols{i},'blobs'),
			if ~isfield(st.vols{i}.blobs{1},'colour'),
				% Add blobs for display using the split colourmap
				scal = 64/(mx-mn);
				dcoff = -mn*scal;
				imgt = imgt*scal+dcoff;

				if isfield(st.vols{i}.blobs{1},'max'),
					mx = st.vols{i}.blobs{1}.max;
				else,
					mx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
					st.vols{i}.blobs{1}.max = mx;
				end;
				if isfield(st.vols{i}.blobs{1},'min'),
					mn = st.vols{i}.blobs{1}.min;
				else,
					mn = min([0 minval(st.vols{i}.blobs{1}.vol)]);
					st.vols{i}.blobs{1}.min = mn;
				end;

				vol  = st.vols{i}.blobs{1}.vol;
				M    = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
                
                tmpt = spm_slice_vol(vol,inv(TM0*(st.Space\M)),TD,[0 NaN])';


				%tmpt_z = find(tmpt==0);tmpt(tmpt_z) = NaN;
				%tmpc_z = find(tmpc==0);tmpc(tmpc_z) = NaN;
				%tmps_z = find(tmps==0);tmps(tmps_z) = NaN;

				sc   = 64/(mx-mn);
				off  = 65.51-mn*sc;
				msk  = find(isfinite(tmpt)); imgt(msk) = off+tmpt(msk)*sc;

				cmap = get(st.fig,'Colormap');

                %figure(st.fig)
                if mn*mx < 0
                    setcolormap('gray-hot-cold')                    
                elseif mx > 0
                    setcolormap('gray-hot');
                else
                    setcolormap('gray-cold')
                end                
                %                redraw_colourbar(i,1,[mn mx],[1:64]'+64); 
			elseif isstruct(st.vols{i}.blobs{1}.colour),
				% Add blobs for display using a defined
                                % colourmap

				% colourmaps
				gryc = [0:63]'*ones(1,3)/63;
				actc = ...
				    st.vols{1}.blobs{1}.colour.cmap;
				actp = ...
				    st.vols{1}.blobs{1}.colour.prop;
				
				% scale grayscale image, not finite -> black
				imgt = scaletocmap(imgt,mn,mx,gryc,65);
				gryc = [gryc; 0 0 0];
				
				% get max for blob image
				vol = st.vols{i}.blobs{1}.vol;
				mat = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				if isfield(st.vols{i}.blobs{1},'max'),
					cmx = st.vols{i}.blobs{1}.max;
				else,
					cmx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
				end;
				if isfield(st.vols{i}.blobs{1},'min'),
					cmn = st.vols{i}.blobs{1}.min;
				else,
					cmn = -cmx;
				end;

				% get blob data
				vol  = st.vols{i}.blobs{1}.vol;
				M    = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				tmpt = spm_slice_vol(vol,inv(TM0*(st.Space\M)),TD,[0 NaN])';
				
				% actimg scaled round 0, black NaNs
				topc = size(actc,1)+1;
				tmpt = scaletocmap(tmpt,cmn,cmx,actc,topc);
				actc = [actc; 0 0 0];
				
				% combine gray and blob data to
				% truecolour
				imgt = reshape(actc(tmpt(:),:)*actp+ ...
					       gryc(imgt(:),:)*(1-actp), ...
					       [size(imgt) 3]);
% 				
%                                 redraw_colourbar(i,1,[cmn cmx],[1:64]'+64); 
				
			else,
				% Add full colour blobs - several sets at once
				scal  = 1/(mx-mn);
				dcoff = -mn*scal;

				wt = zeros(size(imgt));

				imgt  = repmat(imgt*scal+dcoff,[1,1,3]);

				cimgt = zeros(size(imgt));

				for j=1:length(st.vols{i}.blobs), % get colours of all images first
					if isfield(st.vols{i}.blobs{j},'colour'),
						colour(j,:) = reshape(st.vols{i}.blobs{j}.colour, [1 3]);
					else,
						colour(j,:) = [1 0 0];
					end;
				end;
				%colour = colour/max(sum(colour));

				for j=1:length(st.vols{i}.blobs),
					if isfield(st.vols{i}.blobs{j},'max'),
						mx = st.vols{i}.blobs{j}.max;
					else,
						mx = max([eps max(st.vols{i}.blobs{j}.vol(:))]);
						st.vols{i}.blobs{j}.max = mx;
					end;
					if isfield(st.vols{i}.blobs{j},'min'),
						mn = st.vols{i}.blobs{j}.min;
					else,
						mn = min([0 min(st.vols{i}.blobs{j}.vol(:))]);
						st.vols{i}.blobs{j}.min = mn;
					end;

					vol  = st.vols{i}.blobs{j}.vol;
					M    = st.Space\st.vols{i}.premul*st.vols{i}.blobs{j}.mat;
                    tmpt = spm_slice_vol(vol,inv(TM0*M),TD,[0 NaN])';
                    % check min/max of sampled image
                    % against mn/mx as given in st
                    tmpt(tmpt(:)<mn) = mn;
                    tmpt(tmpt(:)>mx) = mx;
                    tmpt = (tmpt-mn)/(mx-mn);
					tmpt(~isfinite(tmpt)) = 0;

					cimgt = cimgt + cat(3,tmpt*colour(j,1),tmpt*colour(j,2),tmpt*colour(j,3));

					wt = wt + tmpt;
                                        cdata=permute(shiftdim([1/64:1/64:1]'* ...
                                                               colour(j,:),-1),[2 1 3]);
                                        %redraw_colourbar(i,j,[mn mx],cdata);
				end;

				imgt = repmat(1-wt,[1 1 3]).*imgt+cimgt;

				imgt(imgt<0)=0; imgt(imgt>1)=1;
			end;
		else,
			scal = 64/(mx-mn);
			dcoff = -mn*scal;
			imgt = imgt*scal+dcoff;
		end;

% 		set(st.vols{i}.ax{1}.d,'HitTest','off', 'Cdata',imgt);
% 		set(st.vols{i}.ax{1}.lx,'HitTest','off',...
% 			'Xdata',[0 TD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
% 		set(st.vols{i}.ax{1}.ly,'HitTest','off',...
% 			'Ydata',[0 TD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));
% 
% 		set(st.vols{i}.ax{2}.d,'HitTest','off', 'Cdata',imgc);
% 		set(st.vols{i}.ax{2}.lx,'HitTest','off',...
% 			'Xdata',[0 CD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
% 		set(st.vols{i}.ax{2}.ly,'HitTest','off',...
% 			'Ydata',[0 CD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));
% 
% 		set(st.vols{i}.ax{3}.d,'HitTest','off','Cdata',imgs);
% 		if st.mode ==0,
% 			set(st.vols{i}.ax{3}.lx,'HitTest','off',...
% 				'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
% 			set(st.vols{i}.ax{3}.ly,'HitTest','off',...
% 				'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(cent(3)-bb(1,3)+1));
% 		else,
% 			set(st.vols{i}.ax{3}.lx,'HitTest','off',...
% 				'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
% 			set(st.vols{i}.ax{3}.ly,'HitTest','off',...
% 				'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(bb(2,2)+1-cent(2)));
% 		end;


	end;
end;

elseif(viewtype == 'c')
    
for i = 1
	M = st.vols{i}.premul*st.vols{i}.mat;
	

	if st.mode ==0,
		SM0 = [	0 0 1 -bb(1,3)+1
			0 1 0 -bb(1,2)+1
			1 0 0 -cent(1)
			0 0 0 1];
		SM = inv(SM0*(st.Space\M)); SD = Dims([3 2]);
	else,
		SM0 = [	0  1 0 -bb(1,2)+1
			0  0 1 -bb(1,3)+1
			1  0 0 -cent(1)
			0  0 0 1];
		SM0 = [	0 -1 0 +bb(2,2)+1
			0  0 1 -bb(1,3)+1
			1  0 0 -cent(1)
			0  0 0 1];
		SM = inv(SM0*(st.Space\M));
		SD = Dims([2 3]);
	end;

	ok=1;

        eval('imgs  = (spm_slice_vol(st.vols{i},SM,SD,st.hld))'';','ok=0;');
    
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.vols{i}.fname);
	else,
                % get min/max threshold
                if strcmp(st.vols{i}.window,'auto')
                        mn = -Inf;
                        mx = Inf;
                else
                        mn = min(st.vols{i}.window);
                        mx = max(st.vols{i}.window);
                end;
                % threshold images
                imgs = max(imgs,mn); imgs = min(imgs,mx);
                % compute intensity mapping, if histeq is available
                if license('test','image_toolbox') == 0
                    st.vols{i}.mapping = 'linear';
                end;
                switch st.vols{i}.mapping,
                 case 'linear',
                 case 'histeq',
                  % scale images to a range between 0 and 1
                  imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                  img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                  imgs = reshape(img(numel(imgt1)+numel(imgc1)+[1:numel(imgs1)]),size(imgs1));
                  mn = 0;
                  mx = 1;
                 case 'quadhisteq',
                  % scale images to a range between 0 and 1
                  imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                  img  = histeq([imgt1(:).^2; imgc1(:).^2; imgs1(:).^2],1024);
                  imgs = reshape(img(numel(imgt1)+numel(imgc1)+[1:numel(imgs1)]),size(imgs1));
                  mn = 0;
                  mx = 1;
                 case 'loghisteq',
                  warning off % messy - but it may avoid extra queries
                  imgs = log(imgs-min(imgs(:)));
                  warning on
                  imgs(~isfinite(imgs)) = 0;
                  % scale log images to a range between 0 and 1
                  imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                  img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                  imgs = reshape(img(numel(imgt1)+numel(imgc1)+[1:numel(imgs1)]),size(imgs1));
                  mn = 0;
                  mx = 1;
                end;
                % recompute min/max for display
                if strcmp(st.vols{i}.window,'auto')
                    mx = -inf; mn = inf;
                end;
                
                if ~isempty(imgs),
			tmp = imgs(isfinite(imgs));
                        mx = max([mx max(max(tmp))]);
                        mn = min([mn min(min(tmp))]);
                end;
                if mx==mn, mx=mn+eps; end;

		if isfield(st.vols{i},'blobs'),
			if ~isfield(st.vols{i}.blobs{1},'colour'),
				% Add blobs for display using the split colourmap
				scal = 64/(mx-mn);
				dcoff = -mn*scal;
				imgs = imgs*scal+dcoff;

				if isfield(st.vols{i}.blobs{1},'max'),
					mx = st.vols{i}.blobs{1}.max;
				else,
					mx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
					st.vols{i}.blobs{1}.max = mx;
				end;
				if isfield(st.vols{i}.blobs{1},'min'),
					mn = st.vols{i}.blobs{1}.min;
				else,
					mn = min([0 minval(st.vols{i}.blobs{1}.vol)]);
					st.vols{i}.blobs{1}.min = mn;
				end;

				vol  = st.vols{i}.blobs{1}.vol;
				M    = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;

                tmps = spm_slice_vol(vol,inv(SM0*(st.Space\M)),SD,[0 NaN])';

				

				%tmpt_z = find(tmpt==0);tmpt(tmpt_z) = NaN;
				%tmpc_z = find(tmpc==0);tmpc(tmpc_z) = NaN;
				%tmps_z = find(tmps==0);tmps(tmps_z) = NaN;

				sc   = 64/(mx-mn);
				off  = 65.51-mn*sc;
				msk  = find(isfinite(tmps)); imgs(msk) = off+tmps(msk)*sc;

				cmap = get(st.fig,'Colormap');

                %figure(st.fig)
                if mn*mx < 0
                    setcolormap('gray-hot-cold')                    
                elseif mx > 0
                    setcolormap('gray-hot');
                else
                    setcolormap('gray-cold')
                end                
                %                redraw_colourbar(i,1,[mn mx],[1:64]'+64); 
			elseif isstruct(st.vols{i}.blobs{1}.colour),
				% Add blobs for display using a defined
                                % colourmap

				% colourmaps
				gryc = [0:63]'*ones(1,3)/63;
				actc = ...
				    st.vols{1}.blobs{1}.colour.cmap;
				actp = ...
				    st.vols{1}.blobs{1}.colour.prop;
				
				% scale grayscale image, not finite -> black
				imgs = scaletocmap(imgs,mn,mx,gryc,65);
				gryc = [gryc; 0 0 0];
				
				% get max for blob image
				vol = st.vols{i}.blobs{1}.vol;
				mat = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				if isfield(st.vols{i}.blobs{1},'max'),
					cmx = st.vols{i}.blobs{1}.max;
				else,
					cmx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
				end;
				if isfield(st.vols{i}.blobs{1},'min'),
					cmn = st.vols{i}.blobs{1}.min;
				else,
					cmn = -cmx;
				end;

				% get blob data
				vol  = st.vols{i}.blobs{1}.vol;
				M    = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				tmps = spm_slice_vol(vol,inv(SM0*(st.Space\M)),SD,[0 NaN])';
				
				% actimg scaled round 0, black NaNs
				topc = size(actc,1)+1;
				tmps = scaletocmap(tmps,cmn,cmx,actc,topc);
				actc = [actc; 0 0 0];
				
				% combine gray and blob data to
				% truecolour

				imgs = reshape(actc(tmps(:),:)*actp+ ...
					       gryc(imgs(:),:)*(1-actp), ...
					       [size(imgs) 3]);
				
                                redraw_colourbar(i,1,[cmn cmx],[1:64]'+64); 
				
			else,
				% Add full colour blobs - several sets at once
				scal  = 1/(mx-mn);
				dcoff = -mn*scal;

				ws = zeros(size(imgs));

				imgs  = repmat(imgs*scal+dcoff,[1,1,3]);

				cimgs = zeros(size(imgs));

				for j=1:length(st.vols{i}.blobs), % get colours of all images first
					if isfield(st.vols{i}.blobs{j},'colour'),
						colour(j,:) = reshape(st.vols{i}.blobs{j}.colour, [1 3]);
					else,
						colour(j,:) = [1 0 0];
					end;
				end;
				%colour = colour/max(sum(colour));

				for j=1:length(st.vols{i}.blobs),
					if isfield(st.vols{i}.blobs{j},'max'),
						mx = st.vols{i}.blobs{j}.max;
					else,
						mx = max([eps max(st.vols{i}.blobs{j}.vol(:))]);
						st.vols{i}.blobs{j}.max = mx;
					end;
					if isfield(st.vols{i}.blobs{j},'min'),
						mn = st.vols{i}.blobs{j}.min;
					else,
						mn = min([0 min(st.vols{i}.blobs{j}.vol(:))]);
						st.vols{i}.blobs{j}.min = mn;
					end;

					vol  = st.vols{i}.blobs{j}.vol;
					M    = st.Space\st.vols{i}.premul*st.vols{i}.blobs{j}.mat;

                    tmps = spm_slice_vol(vol,inv(SM0*M),SD,[0 NaN])';
                    % check min/max of sampled image
                    % against mn/mx as given in st

                    tmps(tmps(:)<mn) = mn;

                    tmps(tmps(:)>mx) = mx;

					tmps = (tmps-mn)/(mx-mn);

					tmps(~isfinite(tmps)) = 0;

					cimgs = cimgs + cat(3,tmps*colour(j,1),tmps*colour(j,2),tmps*colour(j,3));

					ws = ws + tmps;
                                        cdata=permute(shiftdim([1/64:1/64:1]'* ...
                                                               colour(j,:),-1),[2 1 3]);
                                        redraw_colourbar(i,j,[mn mx],cdata);
				end;

				imgs = repmat(1-ws,[1 1 3]).*imgs+cimgs;

				imgs(imgs<0)=0; imgs(imgs>1)=1;
			end;
		else,
			scal = 64/(mx-mn);
			dcoff = -mn*scal;
			imgs = imgs*scal+dcoff;
		end;

% 		set(st.vols{i}.ax{1}.d,'HitTest','off', 'Cdata',imgt);
% 		set(st.vols{i}.ax{1}.lx,'HitTest','off',...
% 			'Xdata',[0 TD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
% 		set(st.vols{i}.ax{1}.ly,'HitTest','off',...
% 			'Ydata',[0 TD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));
% 
% 		set(st.vols{i}.ax{2}.d,'HitTest','off', 'Cdata',imgc);
% 		set(st.vols{i}.ax{2}.lx,'HitTest','off',...
% 			'Xdata',[0 CD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
% 		set(st.vols{i}.ax{2}.ly,'HitTest','off',...
% 			'Ydata',[0 CD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));
% 
% 		set(st.vols{i}.ax{3}.d,'HitTest','off','Cdata',imgs);
% 		if st.mode ==0,
% 			set(st.vols{i}.ax{3}.lx,'HitTest','off',...
% 				'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
% 			set(st.vols{i}.ax{3}.ly,'HitTest','off',...
% 				'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(cent(3)-bb(1,3)+1));
% 		else,
% 			set(st.vols{i}.ax{3}.lx,'HitTest','off',...
% 				'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
% 			set(st.vols{i}.ax{3}.ly,'HitTest','off',...
% 				'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(bb(2,2)+1-cent(2)));
% 		end;

	end;
end;

end



if(viewtype == 's')
    slicedata{kk} = flipdim(imgc,1);
elseif(viewtype == 't')
    slicedata{kk} = flipdim(flipdim(permute(imgt,[2,1,3]),1),2);
elseif(viewtype == 'c')
    slicedata{kk} = flipdim(imgs,1);
end            
pos{kk} = slice;
kk = kk+1;
end
if(isempty(sliceview.colormap))
    if mn*mx < 0
        colormp = 'gray-hot-cold';                    
    elseif mx > 0
        colormp = 'gray-hot';
    else
        colormp = 'gray-cold';
    end  
    sliceview.colormap = colormp;
else
    colormp = sliceview.colormap;
end

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reslice push
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_reslicePush(hObject, eventdata) 

handles = guidata(hObject);

if findstr('SPM2',spm('ver'))
    p = spm_get([],'*.img','Select the image(s) you want to reslice');
else%if findstr('SPM5',spm('ver'))
    p = spm_select([0:200],'IMAGE','Select the image(s) you want to reslice');
end
if findstr('SPM2',spm('ver'))
    target = spm_get([1],'*.img','Select the space-defining image');
else%if findstr('SPM5',spm('ver'))
    target = spm_select([1],'IMAGE','Select the space-defining image');
end

m = size(p,2);
n = length(target);
target1 = target;
template = target1;
if m>n
    target1 = [template blanks(m-n)];
else
    p = [p repmat(blanks(n-m), size(p,1),1)];
end
P = [target1; p];
flags = struct('interp',1,'mask',0,'mean',0,'which',1,'wrap',[0 0 0]');
tic;spm_reslice(P, flags);toc;

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% commonRegion push
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_commonRegionPush(hObject, eventdata) 
handles = guidata(hObject);
if ~isfield(handles,'TF') | length(handles.TF)<=1
    %msgbox('Two or more images need to be loaded to find the common region.', 'info');
    set(handles.infoTextBox, 'string', 'Two or more images need to be loaded to find the common region.'); 
    beep
    return;
end

common = handles.currentDisplayMNI{1};
for ii=2:length(handles.TF)
    common = intersect(common, handles.currentDisplayMNI{ii}, 'rows');
end

if isempty(common)
    %msgbox('No common region found.', 'info');
    set(handles.infoTextBox, 'string', 'No common region found.'); 
    beep
    return;
end

tmpMNI = cell2mat(handles.currentDisplayMNI');
tmpIntensity = cell2mat(handles.currentDisplayIntensity');
intensity = zeros(size(common,1),1);

for ii=1:size(common,1)
    pos = find(abs(tmpMNI(:,1)-common(ii,1))<0.1 & abs(tmpMNI(:,2)-common(ii,2))<0.1 & abs(tmpMNI(:,3)-common(ii,3))<0.1);
    intensity(ii) = 1;%max(abs(tmpIntensity(pos))); %2009/11/18
    %intensity(ii) = prod(tmpIntensity(pos)); %2009/11/18 comment
end

common = [common; cor2mni([1 1 1], handles.M{1}); cor2mni([1 1 2], handles.M{1})];
intensity = [intensity; 0; 1.5];
%intensity = max(intensity) * ones(size(intensity)) + randn(size(intensity))/10;

handles.currentDisplayMNI = {common};
handles.currentDisplayIntensity = {intensity};
[handles.hReg, handles.hSection, handles.hcolorbar] = Draw(handles.currentDisplayMNI, handles.currentDisplayIntensity, hObject, handles);

% if isfield(handles,'hLegend')
%     try
%         set(cell2mat(handles.hLegend'),'visible',{'off'});
%     end
% end
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% small volume push
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_smallVolume(hObject, eventdata) 
handles = guidata(hObject);

mni = cell2mat(handles.currentDisplayMNI');
intensity = cell2mat(handles.currentDisplayIntensity');

xSPM.XYZ = mni2cor(mni, handles.M{1});
xSPM.XYZ = xSPM.XYZ';
xSPM.XYZmm = mni';
xSPM.Z = (intensity');
xSPM.M = handles.M{1};
xSPM.DIM = handles.DIM{1};
xSPM.STAT = handles.TF{1};
if xSPM.STAT == 'T'
    xSPM.STATstr = [xSPM.STAT '_{' num2str(handles.df{1}) '}'];
    xSPM.df = [1 handles.df{1}];
elseif xSPM.STAT == 'F'
    xSPM.STATstr = [xSPM.STAT '_{' num2str(handles.df{1}(1)) ',' num2str(handles.df{1}(2)) '}'];
    xSPM.df = [handles.df{1}];
else
    xSPM.STAT = 'T';
    xSPM.STATstr = [xSPM.STAT '_{' num2str(handles.df{1}) '}'];
    xSPM.df = [1 handles.df{1}];
end
xSPM.k = str2num(get(handles.clusterSizeThresholdEdit, 'string'));
xSPM.u = str2num(get(handles.intensityThresholdEdit, 'string'));
xSPM.u = xSPM.u(1);
xSPM.VOX = abs([xSPM.M(1,1) xSPM.M(2,2) xSPM.M(3,3)]);
xSPM.n = 1;
xSPM.uc = [4.7070 4.2941 338 338];
xSPM.Pc = [0.0000    0.0000    0.0000    0.0668    0.1401    0.1723    0.1723    0.2412    0.2723    0.3092    0.5687    0.5687    0.7008    0.7008];
xSPM.Pp = [0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000 0.0000    0.0000    0.0001    0.0001    0.0001    0.0005    0.0005    0.0015    0.0107    0.0121    0.0123    0.0383    0.0435    0.1019    0.1955    0.2257    0.2966    0.5612    0.6060    0.7725    0.9944];


xSPM.S = length(handles.intensity{1});
xSPM.R = [3 27.0931 276.3307 498.1985];
xSPM.FWHM = [2.9746 3.1923 2.8600];

if get(handles.positiveIntensityRadio, 'Value')
    xSPM.Z = abs(xSPM.Z);
    if xSPM.STAT == 'T';     xSPM.Ps = (1-spm_Tcdf(handles.intensity{1}, xSPM.df(2))); end
    if xSPM.STAT == 'F';     xSPM.Ps = (1-spm_Fcdf(handles.intensity{1}, xSPM.df)); end  
end
    
if get(handles.negativeIntensityRadio, 'Value');
    xSPM.Z = abs(xSPM.Z);    
    if xSPM.STAT == 'T';     xSPM.Ps = (1-spm_Tcdf(-handles.intensity{1}, xSPM.df(2))); end
    if xSPM.STAT == 'F';     xSPM.Ps = (1-spm_Fcdf(-handles.intensity{1}, xSPM.df)); end  
end
   
if get(handles.allIntensityRadio, 'Value');
    if xSPM.STAT == 'T';     xSPM.Ps = (1-spm_Tcdf((handles.intensity{1}), xSPM.df(2))); end
    if xSPM.STAT == 'F';     xSPM.Ps = (1-spm_Fcdf((handles.intensity{1}), xSPM.df)); end  
end

if findstr('SPM2',spm('ver'))
    P = spm_get([0 1],'SPM.mat','locate the corresponding SPM.mat');
else%if findstr('SPM5',spm('ver'))
    P = spm_select([0:1],'SPM.mat','locate the corresponding SPM.mat');
end

if ~isempty(P)
    load(P);
    xSPM.FWHM = SPM.xVol.FWHM;
    xSPM.R =  SPM.xVol.R;
    xSPM.S =  SPM.xVol.S; 
else
    %warndlg(['You did not input SPM.mat. The listed result may not be correct.'], 'SPM.mat missing');
    warndlg('You must input SPM.mat.', 'SPM.mat missing');
    return;
end

if findstr('SPM2',spm('ver'))
    warndlg('Small volume correction does not work in spm 2', 'spm 2');
    return;
elseif findstr('SPM5',spm('ver'))
    spm_VOI(SPM,xSPM,handles.hReg);
elseif findstr('SPM8',spm('ver'))
    spm_VOI8(SPM,xSPM,handles.hReg);        
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% spm_VOI, for spm5, modified by Xu Cui, 2011/03/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TabDat = spm_VOI(SPM,xSPM,hReg)
% List of local maxima and adjusted p-values for a small Volume of Interest
% FORMAT TabDat = spm_VOI(SPM,xSPM,hReg)
%
% SPM   - structure containing analysis details (see spm_spm)
%
% xSPM  - structure containing SPM, distribution & filtering details
%        - required fields are:
% .swd   - SPM working directory - directory containing current SPM.mat
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests
% .STAT  - distribution {Z, T, X or F}
% .df    - degrees of freedom [df{interest}, df{residual}]
% .u     - height threshold
% .k     - extent threshold {resels}
% .XYZ   - location of voxels {voxel coords}
% .XYZmm - location of voxels {mm}
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}
% .M     - voxels - > mm matrix
% .VOX   - voxel dimensions {mm}
% .DIM   - image dimensions {voxels} - column vector
% .Vspm  - Mapped statistic image(s)
% .Ps    - P vlues in searched voxels (for FDR)
%
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% TabDat - Structure containing table data
%        - see spm_list for definition
%
%_______________________________________________________________________
%
% spm_VOI is  called by the SPM results section and takes variables in
% SPM to compute p-values corrected for a specified volume of interest.
%
% The volume of interest may be defined as a box or sphere centred on
% the current voxel or by a mask image.
%
% If the VOI is defined by a mask this mask must have been defined
% independently of the SPM (e.g.using a mask based on an orthogonal
% contrast)
%
% External mask images should be in the same orientation as the SPM
% (i.e. as the input used in stats estimation). The VOI is defined by
% voxels with values greater than 0.
%
% FDR computations are similarly resticted by the small search volume
%
% See also: spm_list
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_VOI.m 112 2005-05-04 18:20:52Z john $


%-Parse arguments
%-----------------------------------------------------------------------
if nargin < 2,   error('insufficient arguments'), end
if nargin < 3,	 hReg = []; end

Num      = 16;			% maxima per cluster
Dis      = 04;			% distance among maxima (mm)

%-Title
%-----------------------------------------------------------------------
spm('FigName',['SPM{',xSPM.STAT,'}: Small Volume Correction']);

%-Get current location {mm}
%-----------------------------------------------------------------------
%xyzmm    = spm_results_ui('GetCoords');
xyzmm = spm_XYZreg('GetCoords',hReg);% added by xu cui

% added by xu cui, open a figure so users can input parameters
f = figure('unit','normalized','position',[0.1 0.1 0.3 0.3]);

%-Specify search volume
%-----------------------------------------------------------------------
str      = sprintf(' at [%.0f,%.0f,%.0f]',xyzmm(1),xyzmm(2),xyzmm(3));
SPACE    = spm_input('Search volume...',-1,'m',...
		{['Sphere',str],['Box',str],'Image'},['S','B','I']);

% voxels in entire search volume {mm}
%-----------------------------------------------------------------------
XYZmm    = SPM.xVol.M(1:3,:)*[SPM.xVol.XYZ; ones(1, SPM.xVol.S)];
Q        = ones(1,size(xSPM.XYZmm,2));
O        = ones(1,size(     XYZmm,2));
FWHM     = xSPM.FWHM;


switch SPACE

	case 'S' %-Sphere
	%---------------------------------------------------------------
	D          = spm_input('radius of VOI {mm}',-2);
	str        = sprintf('%0.1fmm sphere',D);
	j          = find(sum((xSPM.XYZmm - xyzmm*Q).^2) <= D^2);
	k          = find(sum((     XYZmm - xyzmm*O).^2) <= D^2);
	D          = D./xSPM.VOX;


	case 'B' %-Box
	%---------------------------------------------------------------
	D          = spm_input('box dimensions [k l m] {mm}',-2);
	str        = sprintf('%0.1f x %0.1f x %0.1f mm box',D(1),D(2),D(3));
	j          = find(all(abs(xSPM.XYZmm - xyzmm*Q) <= D(:)*Q/2));
	k          = find(all(abs(     XYZmm - xyzmm*O) <= D(:)*O/2));
	D          = D./xSPM.VOX;


	case 'I' %-Mask Image
	%---------------------------------------------------------------
	Msk   = spm_select(1,'image','Image defining search volume');
	D     = spm_vol(Msk);
	str   = sprintf('image mask: %s',spm_str_manip(Msk,'a30'));
	VOX   = sqrt(sum(D.mat(1:3,1:3).^2));
	FWHM  = FWHM.*(xSPM.VOX./VOX);
	XYZ   = D.mat \ [xSPM.XYZmm; ones(1, size(xSPM.XYZmm, 2))];
	j     = find(spm_sample_vol(D, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0);
	XYZ   = D.mat \ [     XYZmm; ones(1, size(    XYZmm, 2))];
	k     = find(spm_sample_vol(D, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0);

end

xSPM.S     = length(k);
xSPM.R     = spm_resels(FWHM,D,SPACE);
xSPM.Z     = xSPM.Z(j);
xSPM.XYZ   = xSPM.XYZ(:,j);
xSPM.XYZmm = xSPM.XYZmm(:,j);
xSPM.Ps    = xSPM.Ps(k);

%-Tabulate p values
%-----------------------------------------------------------------------
str       = sprintf('search volume: %s',str);
if any(strcmp(SPACE,{'S','B'}))
	str = sprintf('%s at [%.0f,%.0f,%.0f]',str,xyzmm(1),xyzmm(2),xyzmm(3));
end

TabDat    = spm_list('List',xSPM,hReg,Num,Dis,str);

%-Reset title
%-----------------------------------------------------------------------
spm('FigName',['SPM{',xSPM.STAT,'}: Results']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% spm_VOI, for spm8, modified by Xu Cui, 2011/02/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TabDat = spm_VOI8(SPM,xSPM,hReg)
% List of local maxima and adjusted p-values for a small Volume of Interest
% FORMAT TabDat = spm_VOI(SPM,xSPM,hReg)
%
% SPM   - structure containing analysis details (see spm_spm)
%
% xSPM  - structure containing SPM, distribution & filtering details
%        - required fields are:
% .swd   - SPM working directory - directory containing current SPM.mat
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests
% .STAT  - distribution {Z, T, X or F}
% .df    - degrees of freedom [df{interest}, df{residual}]
% .u     - height threshold
% .k     - extent threshold {resels}
% .XYZ   - location of voxels {voxel coords}
% .XYZmm - location of voxels {mm}
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}
% .M     - voxels - > mm matrix
% .VOX   - voxel dimensions {mm}
% .DIM   - image dimensions {voxels} - column vector
% .Vspm  - Mapped statistic image(s)
% .Ps    - uncorrected P values in searched volume (for voxel FDR)
% .Pp    - uncorrected P values of peaks (for peak FDR)
% .Pc    - uncorrected P values of cluster extents (for cluster FDR)
%
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% TabDat - Structure containing table data
%        - see spm_list for definition
%
%__________________________________________________________________________
%
% spm_VOI is  called by the SPM results section and takes variables in
% SPM to compute p-values corrected for a specified volume of interest.
%
% The volume of interest may be defined as a box or sphere centred on
% the current voxel or by a mask image.
%
% If the VOI is defined by a mask this mask must have been defined
% independently of the SPM (e.g.using a mask based on an orthogonal
% contrast)
%
% External mask images should be in the same orientation as the SPM
% (i.e. as the input used in stats estimation). The VOI is defined by
% voxels with values greater than 0.
%
% FDR computations are similarly resticted by the small search volume
%
% See also: spm_list
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_VOI.m 2782 2009-02-24 18:53:13Z guillaume $


%-Parse arguments
%--------------------------------------------------------------------------
if nargin < 2,   error('insufficient arguments'), end
if nargin < 3,   hReg = []; end

Num      = 16;          % maxima per cluster
Dis      = 04;          % distance among maxima (mm)

%-Title
%--------------------------------------------------------------------------
spm('FigName',['SPM{',xSPM.STAT,'}: Small Volume Correction']);

%-Get current location {mm}
%--------------------------------------------------------------------------
%xyzmm    = spm_results_ui('GetCoords');
xyzmm = spm_XYZreg('GetCoords',hReg);% added by xu cui

% added by xu cui, open a figure so users can input parameters
f = figure('unit','normalized','position',[0.1 0.1 0.3 0.3]);

%-Specify search volume
%--------------------------------------------------------------------------
str      = sprintf(' at [%.0f,%.0f,%.0f]',xyzmm(1),xyzmm(2),xyzmm(3));
SPACE    = spm_input('Search volume...',-1,'m',...
        {['Sphere',str],['Box',str],'Image'},['S','B','I']);

% voxels in entire search volume {mm}
%--------------------------------------------------------------------------
XYZmm    = SPM.xVol.M(1:3,:)*[SPM.xVol.XYZ; ones(1, SPM.xVol.S)];
Q        = ones(1,size(xSPM.XYZmm,2));
O        = ones(1,size(     XYZmm,2));
FWHM     = xSPM.FWHM;

switch SPACE

    case 'S' %-Sphere
    %----------------------------------------------------------------------
    D          = spm_input('radius of VOI {mm}',-2);
    str        = sprintf('%0.1fmm sphere',D);
    j          = find(sum((xSPM.XYZmm - xyzmm*Q).^2) <= D^2);
    k          = find(sum((     XYZmm - xyzmm*O).^2) <= D^2);
    D          = D./xSPM.VOX;


    case 'B' %-Box
    %----------------------------------------------------------------------
    D          = spm_input('box dimensions [k l m] {mm}',-2);
    if length(D)~=3, D = ones(1,3)*D(1); end
    str        = sprintf('%0.1f x %0.1f x %0.1f mm box',D(1),D(2),D(3));
    j          = find(all(abs(xSPM.XYZmm - xyzmm*Q) <= D(:)*Q/2));
    k          = find(all(abs(     XYZmm - xyzmm*O) <= D(:)*O/2));
    D          = D./xSPM.VOX;


    case 'I' %-Mask Image
    %----------------------------------------------------------------------
    Msk   = spm_select(1,'image','Image defining search volume');
    D     = spm_vol(Msk);
    str   = strrep(spm_str_manip(Msk,'a30'),'\','\\');
    str   = strrep(str,'^','\^'); str   = strrep(str,'_','\_');
    str   = strrep(str,'{','\{'); str   = strrep(str,'}','\}');
    str   = sprintf('image mask: %s',str); 
    VOX   = sqrt(sum(D.mat(1:3,1:3).^2));
    FWHM  = FWHM.*(xSPM.VOX./VOX);
    XYZ   = D.mat \ [xSPM.XYZmm; ones(1, size(xSPM.XYZmm, 2))];
    j     = find(spm_sample_vol(D, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0);
    XYZ   = D.mat \ [     XYZmm; ones(1, size(    XYZmm, 2))];
    k     = find(spm_sample_vol(D, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0);

end

xSPM.S     = length(k);
xSPM.R     = spm_resels(FWHM,D,SPACE);
xSPM.Z     = xSPM.Z(j);
xSPM.XYZ   = xSPM.XYZ(:,j);
xSPM.XYZmm = xSPM.XYZmm(:,j);

%-Restrict FDR to the search volume
%--------------------------------------------------------------------------
df         = xSPM.df;
STAT       = xSPM.STAT;
DIM        = xSPM.DIM;
R          = xSPM.R;
n          = xSPM.n;
Z          = xSPM.Z;
u          = xSPM.u;
S          = xSPM.S;

try, xSPM.Ps  = xSPM.Ps(k); end
[up, xSPM.Pp] = spm_uc_peakFDR(0.05,df,STAT,R,n,Z,SPM.xVol.XYZ(:,k),u);
try % if STAT == 'T'
    V2R               = 1/prod(xSPM.FWHM(DIM>1));
    [uc, xSPM.Pc, ue] = spm_uc_clusterFDR(0.05,df,STAT,R,n,Z,SPM.xVol.XYZ(:,k),V2R,u);
catch
    uc                = NaN;
    ue                = NaN;
    xSPM.Pc           = [];
end
uu            = spm_uc(0.05,df,STAT,R,n,S);
xSPM.uc       = [uu up ue uc];

%-Tabulate p values
%--------------------------------------------------------------------------
str       = sprintf('search volume: %s',str);
if any(strcmp(SPACE,{'S','B'}))
    str = sprintf('%s at [%.0f,%.0f,%.0f]',str,xyzmm(1),xyzmm(2),xyzmm(3));
end

TabDat    = spm_list8('List',xSPM,hReg,Num,Dis,str);

%-Reset title
%--------------------------------------------------------------------------
spm('FigName',['SPM{',xSPM.STAT,'}: Results']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% volume push
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_volumePush(hObject, eventdata) 
handles = guidata(hObject);
mni = cell2mat(handles.currentDisplayMNI');
intensity = cell2mat(handles.currentDisplayIntensity');

xSPM.XYZ = mni2cor(mni, handles.M{1});
xSPM.XYZ = xSPM.XYZ';
xSPM.XYZmm = mni';
xSPM.Z = (intensity');
xSPM.M = handles.M{1};
xSPM.DIM = handles.DIM{1};
xSPM.STAT = handles.TF{1};
if xSPM.STAT == 'T'
    xSPM.STATstr = [xSPM.STAT '_{' num2str(handles.df{1}) '}'];
    xSPM.df = [1 handles.df{1}];
elseif xSPM.STAT == 'F'
    xSPM.STATstr = [xSPM.STAT '_{' num2str(handles.df{1}(1)) ',' num2str(handles.df{1}(2)) '}'];
    xSPM.df = [handles.df{1}];
else
    xSPM.STAT = 'T';
    xSPM.STATstr = [xSPM.STAT '_{' num2str(handles.df{1}) '}'];
    xSPM.df = [1 handles.df{1}];
end
xSPM.k = str2num(get(handles.clusterSizeThresholdEdit, 'string'));
xSPM.u = str2num(get(handles.intensityThresholdEdit, 'string'));
xSPM.u = xSPM.u(1);
xSPM.VOX = abs([xSPM.M(1,1) xSPM.M(2,2) xSPM.M(3,3)]);
xSPM.n = 1;
xSPM.uc = [4.7070 4.2941 338 338];
xSPM.Pc = [0.0000    0.0000    0.0000    0.0668    0.1401    0.1723    0.1723    0.2412    0.2723    0.3092    0.5687    0.5687    0.7008    0.7008];
xSPM.Pp = [0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000 0.0000    0.0000    0.0001    0.0001    0.0001    0.0005    0.0005    0.0015    0.0107    0.0121    0.0123    0.0383    0.0435    0.1019    0.1955    0.2257    0.2966    0.5612    0.6060    0.7725    0.9944];

% SPM    - structure containing SPM, distribution & filtering details
%        - required fields are:
% .swd   - SPM working directory - directory containing current SPM.mat
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests        
% .STAT  - distribution {Z, T, X or F}     
% .df    - degrees of freedom [df{interest}, df{residual}]
% .u     - height threshold
% .k     - extent threshold {voxels}
% .XYZ   - location of voxels {voxel coords}
% .XYZmm - location of voxels {mm}
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}     
% .M     - voxels - > mm matrix
% .VOX   - voxel dimensions {mm}
% .Vspm  - mapped statistic image(s)
% .Ps    - uncorrected P values in searched volume (for FDR)

xSPM.S = length(handles.intensity{1});
xSPM.R = [3 27.0931 276.3307 498.1985];
xSPM.FWHM = [2.9746 3.1923 2.8600];

if get(handles.positiveIntensityRadio, 'Value')
    xSPM.Z = abs(xSPM.Z);
    if xSPM.STAT == 'T';     xSPM.Ps = (1-spm_Tcdf(handles.intensity{1}, xSPM.df(2))); end
    if xSPM.STAT == 'F';     xSPM.Ps = (1-spm_Fcdf(handles.intensity{1}, xSPM.df)); end  
end
    
if get(handles.negativeIntensityRadio, 'Value');
    xSPM.Z = abs(xSPM.Z);    
    if xSPM.STAT == 'T';     xSPM.Ps = (1-spm_Tcdf(-handles.intensity{1}, xSPM.df(2))); end
    if xSPM.STAT == 'F';     xSPM.Ps = (1-spm_Fcdf(-handles.intensity{1}, xSPM.df)); end  
end
   
if get(handles.allIntensityRadio, 'Value');
    if xSPM.STAT == 'T';     xSPM.Ps = (1-spm_Tcdf((handles.intensity{1}), xSPM.df(2))); end
    if xSPM.STAT == 'F';     xSPM.Ps = (1-spm_Fcdf((handles.intensity{1}), xSPM.df)); end  
end

if findstr('SPM2',spm('ver'))
    P = spm_get([0 1],'SPM.mat','locate the corresponding SPM.mat');
else%if findstr('SPM5',spm('ver'))
    P = spm_select([0:1],'SPM.mat','locate the corresponding SPM.mat');
end

if ~isempty(P)
    load(P);
    xSPM.FWHM = SPM.xVol.FWHM;
    xSPM.R =  SPM.xVol.R;
    xSPM.S =  SPM.xVol.S; 
else
    warndlg(['You did not input SPM.mat. The listed result may not be correct.'], 'SPM.mat missing');
end

if findstr('SPM8',spm('ver'))
    spm_list8('List',xSPM,handles.hReg);
else%if findstr('SPM5',spm('ver'))
    spm_list('List',xSPM,handles.hReg);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_displayPush(hObject, eventdata) 
handles = guidata(hObject);
spm_image('init', handles.imageFileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% all in one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_allinonePush(hObject, eventdata) 
try
    if findstr('SPM2',spm('ver'))
        P = spm_get([0:100],'*IMAGE','Select image files');
    else%if findstr('SPM5',spm('ver'))
        P = spm_select(Inf,'image','Select image files');
    end
    
catch
    return;
end
if isempty(P)
    return
end

handles = guidata(hObject);
%if ~isfield(handles,'hReg') | ~isfield(handles,'hSection') | ~isfield(handles,'hcolorbar')
    tmp = [0 0 0];
    [handles.hReg, handles.hSection, handles.hcolorbar] = Draw(tmp([],:), [], hObject, handles);
    %end

colors = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
for ii=1:size(P,1)
    thisfilename = deblank(P(ii,:));
    [handles.imageFileName, handles.M, handles.DIM, handles.TF, handles.df, handles.mni, handles.intensity] = getImageFile(thisfilename);
    cor = mni2cor(handles.mni, handles.M);
    spm_orthviews('addcolouredblobs',1,cor',handles.intensity',handles.M,colors(mod(ii,6)+1,:));
end
spm_orthviews('Redraw');

guidata(hObject, handles);

% contents = get(handles.sectionViewListbox,'String');
% currentsel = contents{get(handles.sectionViewListbox,'Value')};
% sectionViewTargetFile = getSectionViewTargetFile(handles.spmdir, currentsel);
% 
% % spm_image('init', sectionViewTargetFile);
% % spm_image('addblobs');
% % return

% guidata(hObject, handles);
% load xSPM;
%         VOL.XYZ = xSPM.XYZ;
%         VOL.Z = xSPM.Z;
%         VOL.M = handles.M;

%addcolouredimage(handles.hSection, '333.img',[1 0 1]);
% uigetfiles
% nblobs = 4;
% 	for i=1:nblobs,
% 		%[SPM,VOL] = spm_getSPM;
% 		%c = spm_input('Colour','+1','m','Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
%         c = i;
%         VOL.XYZ = ceil(rand(3,20)*20);
%         VOL.Z = randn(1,20);
%         VOL.M = handles.M;
% 		colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
% 		spm_orthviews('addcolouredblobs',1,VOL.XYZ,VOL.Z,VOL.M,colours(c,:));
%     end;
	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% overlay a brain region Push
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_overlayPush(hObject, eventdata)    
handles = guidata(hObject);

tobeoverlay = deblank(get(handles.overlayEdit, 'string'));

if isempty(tobeoverlay)
    return;
end

tobeoverlay = str2cell(tobeoverlay, ' ');
if isnumeric(tobeoverlay{1})
    if length(tobeoverlay{1})>1
        errordlg(['Your input, ' deblank(get(handles.overlayEdit, 'string')) ', seems coincide with a matlab constant.'], 'error');
        return
    end
    for ii=1:length(tobeoverlay)
        tobeoverlay{ii} = num2str(tobeoverlay{ii});
    end
end
fn = fieldnames(handles.wholeMaskMNIAll);
for ii=1:length(tobeoverlay)
    pos{ii} = [];
    for jj=1:length(fn)
%         x = [];
%         if ~isempty(str2num(tobeoverlay{ii}))
%             y = str2cell(fn{jj},'_');
%             for kk=1:length(y)
%                 x = isequal(tobeoverlay{ii}, y{kk});
%                 if x; break;end;
%             end
%             if x==0; x = []; end;
%         else
%             x = findstr(lower(tobeoverlay{ii}), lower(fn{jj}));
%         end
        if strcmp(tobeoverlay{ii}, fn{jj})
            pos{ii} = [pos{ii} jj];
        end
    end
end
common = pos{1};
for ii=2:length(pos)
    common = intersect(common, pos{ii});
end

if isempty(common)
    %warndlg(['I don'' find ' deblank(get(handles.overlayEdit, 'string')) '.'], 'oops');
    set(handles.infoTextBox, 'string', ['I don'' find ' deblank(get(handles.overlayEdit, 'string')) '.']); 
    return
end

mask = [];
for ii=1:length(common)
    eval(['mask = [mask; handles.wholeMaskMNIAll.' fn{common(ii)} '];']);
end

if isempty(mask)
    %warndlg(['I don'' find ' deblank(get(handles.overlayEdit, 'string')) '.'], 'oops');
    set(handles.infoTextBox, 'string', ['I don'' find ' deblank(get(handles.overlayEdit, 'string')) '.']); 
    return;
end

try
    handles.mni;
catch
    delete(gcf);
    xjview(mask);
    return;
end

handles.imageFileName = [handles.imageFileName, {deblank(get(handles.overlayEdit, 'string'))}];
handles.M = [handles.M, {handles.M{1}}];
handles.DIM = [handles.DIM, {handles.DIM{1}}];
handles.currentDisplayMNI = [handles.currentDisplayMNI, {mask}];
m = max(abs(cell2mat(handles.currentDisplayIntensity')));
if ~isempty(m)
    handles.currentDisplayIntensity = [handles.currentDisplayIntensity, {-2*m*ones(size(mask,1),1)}];
else
    handles.currentDisplayIntensity = [handles.currentDisplayIntensity, {ones(size(mask,1),1)}];
end
[handles.hReg, handles.hSection, handles.hcolorbar] = Draw(handles.currentDisplayMNI, handles.currentDisplayIntensity, hObject, handles);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% overlay a brain region Edit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_overlayEdit(hObject, eventdata)    
handles = guidata(hObject);
CallBack_overlayPush(handles.overlayPush, eventdata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% overlay a brain region popup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_overlayPop(hObject, eventdata)    
handles = guidata(hObject);
names = get(hObject, 'string');
value = get(hObject, 'value');
set(handles.overlayEdit, 'string', names{value});
CallBack_overlayEdit(handles.overlayEdit, eventdata);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% search xBrain.org and other databases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CallBack_searchPush(hObject, eventdata)    

try
	handles = guidata(hObject);
	searchContent = get(handles.searchContentEdit,'string');
	searchEngine = get(handles.searchEnginePop,'value');
    xbrainSearchField = 'region';
	if isempty(deblank(searchContent))
        xyz = spm_XYZreg('GetCoords',handles.hReg);
        handles.currentxyz = xyz';        
        searchContent = num2str(xyz');
        set(handles.searchContentEdit,'string',searchContent);
        guidata(hObject, handles);
        xbrainSearchField = 'mni or tal&mnidistance=20';
    end
    
	if searchEngine == 1
        urlstr = ['http://people.hnl.bcm.tmc.edu/cuixu/cgi-bin/bmd/paper.pl?search_content=' searchContent '&search_field=' xbrainSearchField];
	elseif searchEngine == 2
        urlstr = ['http://scholar.google.com/scholar?q=%22' searchContent '%22'];
	elseif searchEngine == 3
        urlstr = ['http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=PureSearch&db=pubmed&details_term=%22' searchContent '%22'];
	elseif searchEngine == 4
        urlstr = ['http://en.wikipedia.org/w/index.php?search=%22' searchContent '%22'];
	end
catch
    urlstr = 'http://www.xbrain.org';
end

web(urlstr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% control panel on or off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function controlHide(handles, status)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get multiple values from a string deliminated delim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = str2cell(str, delim)
if ~exist('delim')
    delim = '; ';
end

[out{1},b] = strtok(str, delim);
ii = 2;
while ~isempty(b)
    [out{ii},b] = strtok(b, delim);
    ii = ii+1;
end

for ii=1:length(out)
    out2{ii} = str2num(out{ii});
    if isempty(out2{ii})
        return;
    end
end

for ii=1:length(out)
    out{ii} = str2num(out{ii});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cell2str
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = cell2str(acell, delim)
if ~exist('delim')
    delim = ';';
end

str = [];

if length(acell)==1
    if isstr(acell{1})
        str = acell{1};
    else
        str = num2str(acell{1});
    end
else
    for ii=1:length(acell)
        if isstr(acell{ii})
            str = [str acell{ii}, delim ' '];
        else
            str = [str num2str(acell{ii}), delim ' '];
        end
    end
    str(end-1:end)=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cellmax, find max in each element, return a cell 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MAXX = cellmax(acell, absolute)

if ~exist('absolute')
    absolute = '';
end

MAXX = [];
for ii=1:length(acell)
    if strfind('abs',absolute)
        MAXX = [MAXX max(abs(acell{ii}))];
    else
        MAXX = [MAXX max(acell{ii})];
    end
end
MAXX = num2cell(MAXX);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% maxcell
%%% find max of all numbers in the whole cell, return a single value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MAXX = maxcell(acell, absolute)

if ~exist('absolute')
    absolute = '';
end
MAXX = cellmax(acell, absolute);
MAXX = max(cell2mat(MAXX));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% p2t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = p2t(p, df, TF)
if ~iscell(p)
	if upper(TF)=='T' | upper(TF)=='S'
        t = spm_invTcdf(1-p,df);
	elseif upper(TF) == 'F'
        t = spm_invFcdf(1-p,df);
	end
else
    for ii=1:length(p)
        t{ii} = p2t(p{ii},df{ii},TF{ii});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% s2t, s is defined as -log10(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = s2t(s, df, TF)

if ~iscell(s)
    p = 10^(-s);
	if upper(TF)=='T' | upper(TF)=='S'
        t = spm_invTcdf(1-p,df);
	elseif upper(TF) == 'F'
        t = spm_invFcdf(1-p,df);
	end
else
    for ii=1:length(p)
        t{ii} = s2t(s{ii},df{ii},TF{ii});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% t2p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = t2p(t, df, TF)
if ~iscell(t)
	if upper(TF)=='T' | upper(TF)=='S'
        p = 1-spm_Tcdf(t,df);
	elseif upper(TF) == 'F'
        p = 1-spm_Fcdf(t,df);
	end
else
    for ii=1:length(t)
        p{ii} = t2p(t{ii},df{ii},TF{ii});
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% t2s, s is defined as -log10(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = t2s(t, df, TF)
if ~iscell(t)
	if upper(TF)=='T' | upper(TF)=='S'
        p = 1-spm_Tcdf(t,df);
	elseif upper(TF) == 'F'
        p = 1-spm_Fcdf(t,df);
    else
        p = 0.1;
	end
    s = -log10(p); 
else
    s = {10};
    for ii=1:length(t)
        s{ii} = t2s(t{ii},df{ii},TF{ii});
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mni2mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mask = mni2mask(coords, targetfilename, intensity, M, DIM, templateFile, isMask)
% function mask = mni2mask(coords, targetfilename, intensity, templateFile)
% make mask from coordinates
%
% coords: a Nx3 column of 3-d coordinates in MNI
% space
% targetfilename: (optional) the image files to be written
% intensity: (optional) Nx1, the values of each coordinate.
% M: rotation matrix
% DIM: dimension
% templateFile: (optional) the templateFile from which we can get the right
% dimensions.
% isMask: if this variable exist and equal to 1, then all non-zero
% intensities will be set to 1
%
% Xu Cui
% 2004/11/18

if ~exist('intensity')
    intensity = ones(size(coords,1),1);
end
thistemplateFile = '';
if exist('templateFile')
    if ~isempty(templateFile)
        thistemplateFile = templateFile;
    end
end

if isempty(thistemplateFile)
    V.mat = [...
        -4     0     0    84; ...
         0     4     0  -116; ...
         0     0     4   -56; ...
         0     0     0     1];
    V.dim = [41 48 35 16];
    V.dt = [4,0];
    if exist('M')
        V.mat = M;
    end
    if exist('DIM')
        V.dim = DIM;
        if findstr('SPM2',spm('ver'))
            V.dim(4) = 16;
        end   
    end
    V.fname = targetfilename;
    V.descrip = 'our own image';
else
    V = spm_vol(templateFile);
    V.fname = targetfilename;
    if isfield(V, 'descrip')
        V.descrip = ['my image from ' V.descrip];
    else
        V.descrip = 'my own image';
    end
end

thisismask = 0;
if exist('isMask')
    if isMask == 1
        thisismask = 1;
    end
end
if thisismask
    V.descrip = 'my mask';
    intensity = ones(size(intensity));
end

O = zeros(V.dim(1),V.dim(2),V.dim(3));

coords = mni2cor(coords,V.mat);


for ii=1:size(coords,1)
    O(coords(ii,1),coords(ii,2),coords(ii,3)) = intensity(ii);
end

if exist('targetfilename')
    V = spm_write_vol(V,O);
end

mask = O;

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get image file information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [imageFile,M,DIM,TF,df,mni,intensity] = getImageFile(thisfilename)
% function imageFile = getImageFile(filename)
% get the information of this/these image file(s)
%
% thisfilename: (optional) if doesn't give file name, then show a open file
% dialog
% imageFile: the full name of the selected file (if user pressed cancel,
% imageFile == 0
% M: M matrix (rotation matrix)
% DIM: dimension
% TF: t test or f test? 'T' or 'F'
% df: degree of freedome
% mni: mni coord
% intensity: intensity of each mni coord
%
% Note: The returned variables are cellarrays.
%
% Xu Cui
% last revised: 2005-05-03

if nargin < 1 | isempty(thisfilename)
    if findstr('SPM2',spm('ver'))
        P0 = spm_get([0:100],'*IMAGE','Select image files');
    else%if findstr('SPM5',spm('ver'))
        P0 = spm_select(Inf,'image','Select image files');
    end    

%     try
%         P0 = spm_get([0:100],'*IMAGE','Select image files');
%     catch
%         P0 = [];
%         [FileName,PathName] = uigetfile({'*.img';'*.IMG';'*.*'},'Select image files','MultiSelect','on');
%         if isstr(FileName)
%             P0 = {fullfile(PathName, FileName)};
%         elseif iscellstr(FileName)
%             for ii=1:length(FileName)
%                 P0{ii} = fullfile(PathName, FileName{ii});
%             end
%         else
%             P0 = [];
%         end
%     end
    try
        if isempty(P0)
            imageFile = '';M=[];DIM=[];TF=[];df=[];mni=[];intensity=[];
            return
        end
    end
    for ii=1:size(P0,1)
        P{ii} = deblank(P0(ii,:));
    end
else
    if isstr(thisfilename)
        P = {thisfilename};
    elseif iscellstr(thisfilename)
        P = thisfilename;
    else
        disp('Error: In getImageFile: I don''t understand the input.');
        imageFile = '';M=[];DIM=[];TF=[];df=[];mni=[];intensity=[];
        return
    end
end

global LEFTRIGHTFLIP_;

for ii=1:length(P)
    imageFile{ii} = P{ii};
	[cor{ii}, intensity{ii}, tmp{ii}, M{ii}, DIM{ii}, TF{ii}, df{ii}] = mask2coord(imageFile{ii}, 0);
    if LEFTRIGHTFLIP_ == 1
        M{ii}(1,:) = -M{ii}(1,:);
    end
	mni{ii} = cor2mni(cor{ii}, M{ii});
end














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cuixuFindStructure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [onelinestructure, cellarraystructure] = cuixuFindStructure(mni, DB)
% function [onelinestructure, cellarraystructure] = cuixuFindStructure(mni, DB)
%
% this function converts MNI coordinate to a description of brain structure
% in aal
%
%   mni: the coordinates (MNI) of some points, in mm.  It is Nx3 matrix
%   where each row is the coordinate for one point
%   LDB: the database.  This variable is available if you load
%   TDdatabase.mat
%
%   onelinestructure: description of the position, one line for each point
%   cellarraystructure: description of the position, a cell array for each point
%
%   Example:
%   cuixuFindStructure([72 -34 -2; 50 22 0], DB)
%
% Xu Cui
% 2007-11-20
%

N = size(mni, 1);

% round the coordinates
mni = round(mni/2) * 2;

T = [...
     2     0     0   -92
     0     2     0  -128
     0     0     2   -74
     0     0     0     1];

index = mni2cor(mni, T);

cellarraystructure = cell(N, length(DB));
onelinestructure = cell(N, 1);

for ii=1:N
    for jj=1:length(DB)
        graylevel = DB{jj}.mnilist(index(ii, 1), index(ii, 2),index(ii, 3));
        if graylevel == 0
            thelabel = 'undefined';
        else
            if jj==length(DB); tmp = ' (aal)'; else tmp = ''; end
            thelabel = [DB{jj}.anatomy{graylevel} tmp];
        end
        cellarraystructure{ii, jj} = thelabel;
        onelinestructure{ii} = [ onelinestructure{ii} ' // ' thelabel ];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mask2coord
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cor, intensity, cor_singlecolumn,M,DIM,TF,df] = mask2coord(mask, checkdimension)
% [cor, intensity, cor_singlecolumn] = mask2coord(mask, checkdimension)
%
% This is to retrieve the coordinate of a mask file, or a matrix of 3-D
%
% mask: an image file or a matrix (3-D), with some of the elements are
% non-zeros
% checkdimension: check if the dimension is checkdimension, if not, return empty
% matrix
% cor: a N*3 matrix, which each row a coordinate in matrix
% intensity: a N*1 matrix, which encodes the intensity of each voxel.
% cor_singlecolumn: a N*1 matrix, each row is the index in the matrix
% M: rotation matrix
% DIM: dimension
% TF: t test or f test? 'T','F' or []
% df: degree of freedome for T/F test
%
% Example:
%   mask = zeros(4,3,2);
%   mask(1,2,1) = 1;
%   mask(3,2,2) = 1;
%   mask2coord(mask)
%
%   mask2coord('spmT_0002.img')
%   mask2coord('spmT_0002.img',[41 48 35])
%
% Xu Cui
% 2004-9-20
% last revised: 2005-04-30

if nargin < 2
    checkdimension = 0;
end

if ischar(mask)
    V = spm_vol(mask);
    mask = spm_read_vols(V);
    M = V.mat;
    DIM = V.dim;
    TF = 'T';
	T_start = strfind(V.descrip,'SPM{T_[')+length('SPM{T_[');
    if isempty(T_start); T_start = strfind(V.descrip,'SPM{F_[')+length('SPM{F_['); TF='F'; end
    if isempty(T_start)
        TF=[]; df=[];
    else
    	T_end = strfind(V.descrip,']}')-1;
    	df = str2num(V.descrip(T_start:T_end));    
    end
else
    M = [];
    TF = [];
    df = [];
end

dim = size(mask);
if length(checkdimension)==3
    if dim(1)~= checkdimension(1) | dim(2)~= checkdimension(2) | dim(3)~= checkdimension(3)
        y = [];
        disp('dimension doesn''t match')
        return
    end
end

pos = find(mask~=0);
intensity = mask(pos);

y = zeros(length(pos),3);

y(:,3) = ceil(pos/(dim(1)*dim(2)));
pos = pos - (y(:,3)-1)*(dim(1)*dim(2));
y(:,2) = ceil(pos/dim(1));
pos = pos - (y(:,2)-1)*(dim(1));
y(:,1) = pos;

cor = y;
cor_singlecolumn = pos;
DIM = dim;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mni2cor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coordinate = mni2cor(mni, T)
% function coordinate = mni2cor(mni, T)
% convert mni coordinate to matrix coordinate
%
% mni: a Nx3 matrix of mni coordinate
% T: (optional) transform matrix
% coordinate is the returned coordinate in matrix
%
% caution: if T is not specified, we use:
% T = ...
%     [-4     0     0    84;...
%      0     4     0  -116;...
%      0     0     4   -56;...
%      0     0     0     1];
%
% xu cui
% 2004-8-18
%

if isempty(mni)
    coordinate = [];
    return;
end

if nargin == 1
	T = ...
        [-4     0     0    84;...
         0     4     0  -116;...
         0     0     4   -56;...
         0     0     0     1];
end

coordinate = [mni(:,1) mni(:,2) mni(:,3) ones(size(mni,1),1)]*(inv(T))';
coordinate(:,4) = [];
coordinate = round(coordinate);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cor2mni
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mni = cor2mni(cor, T)
% function mni = cor2mni(cor, T)
% convert matrix coordinate to mni coordinate
%
% cor: an Nx3 matrix
% T: (optional) rotation matrix
% mni is the returned coordinate in mni space
%
% caution: if T is not given, the default T is
% T = ...
%     [-4     0     0    84;...
%      0     4     0  -116;...
%      0     0     4   -56;...
%      0     0     0     1];
%
% xu cui
% 2004-8-18
% last revised: 2005-04-30

if nargin == 1
    T = ...
        [-4     0     0    84;...
         0     4     0  -116;...
         0     0     4   -56;...
         0     0     0     1];
end

cor = round(cor);
mni = T*[cor(:,1) cor(:,2) cor(:,3) ones(size(cor,1),1)]';
mni = mni';
mni(:,4) = [];
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% draw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hReg, hSection, hcolorbar] = Draw(mniCoord, intensity, hObject, handles)

try
    delete(handles.hcolorbar);
end
try
    delete(handles.hReg);    
end
try
    delete(handles.hSection);
end
try
    cla(handles.glassViewAxes);
end

try
	H  = findobj(get(handles.figure,'Children'),'flat','Type','axes');
	un = cellstr(get(H,'Units'));
	pos = get(H,'position');
	
	for ii=1:length(H)
        if findstr('pixels', un{ii})
            continue;
        end
        if pos{ii}(1)>0.4
            delete(H(ii));
        end
	end
end

sectionViewTargetFile = handles.sectionViewTargetFile;
rendStyle = get(handles.renderStylePop,'value');
if(rendStyle == 1)
    rendStyle = 'new';
else    
    rendStyle = 'old';
end


if ~iscell(mniCoord) | (iscell(mniCoord) & length(mniCoord)==1)% multiple input? no
    if (iscell(mniCoord) & length(mniCoord)==1)
        mniCoord = mniCoord{1};
        intensity = intensity{1};
        handles.M = handles.M{1};
        handles.DIM = handles.DIM{1};
    end
	if max(intensity)*min(intensity) < 0
        [hReg, hSection, hcolorbar] = cuixuSectionView(mniCoord,intensity,sectionViewTargetFile,hObject,handles);    
    else
        [hReg, hSection, hcolorbar] = cuixuSectionView(mniCoord,abs(intensity),sectionViewTargetFile,hObject,handles); 
	end
	
	if size(mniCoord,1)>1
        pos1 = find(intensity>=0);
        pos2 = find(intensity<0);
        if ~isempty(pos1) & ~isempty(pos2)
            if get(handles.renderViewCheck, 'Value'); cuixuRenderView(rendStyle, mniCoord(pos1,:),intensity(pos1,:),mniCoord(pos2,:),-intensity(pos2,:)  ); end
        else
            if get(handles.renderViewCheck, 'Value'); cuixuRenderView(rendStyle, mniCoord, abs(intensity)); end        
        end
	end   
else % multiple input? yes
    [hReg, hSection, hcolorbar] = cuixuSectionView(mniCoord,intensity,sectionViewTargetFile,hObject,handles);    

    mniCoordtmp=[];
    intensitytmp=[];
    for ii=1:length(mniCoord)
        mniCoordtmp = [mniCoordtmp; mniCoord{ii}];
        intensitytmp = [intensitytmp; intensity{ii}];
    end
    mniCoord = mniCoordtmp;
    intensity = intensitytmp;
	if size(mniCoord,1)>1
        pos1 = find(intensity>=0);
        pos2 = find(intensity<0);
        if ~isempty(pos1) & ~isempty(pos2)
            if get(handles.renderViewCheck, 'Value'); cuixuRenderView(rendStyle, mniCoord(pos1,:),intensity(pos1,:),mniCoord(pos2,:),-intensity(pos2,:)); end
        else
            if get(handles.renderViewCheck, 'Value'); cuixuRenderView(rendStyle, mniCoord, abs(intensity)); end        
        end
	end    
    
end
    
try
    spm_XYZreg('SetCoords',handles.currentxyz',hReg);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cuixuSectionView
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hReg, hSection, hcolorbar] = cuixuSectionView(mniCoord, intensity, targetFile, hObject,handles)
% function h = cuixuSectionView(mniCoord, intensity)
% This is to project your coordinate to section view
%
% mniCoord: the mni coordinate, N*3 matrix
% intensity: the plot intensity of the spots, 1*N matrix. The intensity
% could be t value, for example.
%
% a special feature of this function is: it automatically seperate the
% negative and positive intensity and use hot and cold color to represent
% them.
%
% h: the returned handle for the axes
%
% SEE ALSO: cuixuView cuixuGlassView cuixuRenderView
%
% Xu Cui
% 2004/11/11

if ~iscell(mniCoord) | (iscell(mniCoord) & length(mniCoord)==1)% multiple input? no
    multiple = 0;
    if (iscell(mniCoord) & length(mniCoord)==1)
        mniCoord = mniCoord{1};
        intensity = intensity{1};
        handles.M = handles.M{1};
        handles.DIM = handles.DIM{1};
    end
else
    multiple = 1;
    mniCoordtmp = [];
    intensitytmp = [];
    for ii=1:length(mniCoord)
        mniCoordtmp = [mniCoordtmp; mniCoord{ii}];
        intensitytmp = [intensitytmp; intensity{ii}];
        % for multiple input
        cSPM{ii}.XYZ =  mni2cor(mniCoord{ii}, handles.M{ii});
        cSPM{ii}.XYZ = cSPM{ii}.XYZ';
        cSPM{ii}.Z = abs(intensity{ii});
        cSPM{ii}.M = handles.M{ii};
        cSPM{ii}.DIM = handles.DIM{ii}';
    end
    mniCoord = mniCoordtmp;
    intensity = intensitytmp;
	handles.M = handles.M{1};
	handles.DIM = handles.DIM{1}; 
end

SPM.XYZ = mni2cor(mniCoord, handles.M);
SPM.XYZ = SPM.XYZ';
SPM.Z = intensity;
SPM.M = handles.M;
SPM.DIM = handles.DIM';

%%%%%%%%%%%%%%% Reg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
handles.colormap = colormap;

coor = mniCoord;
xSPM = SPM;
xSPM.XYZmm = mniCoord';
axes(handles.glassViewAxes)
WS     = spm('WinScale');
FS     = spm('FontSizes');

Finter = gcf;
hReg    = uicontrol(Finter,'Style','Frame','Position',[60 100 300 300].*WS,...
		'Visible','off');
[hReg,xyz] = spm_XYZreg('InitReg',hReg,xSPM.M,xSPM.DIM,[0;0;0]);

hFxyz      = spm_results_ui('DrawXYZgui',xSPM.M,xSPM.DIM,xSPM,xyz,Finter);
spm_XYZreg('XReg',hReg,hFxyz,'spm_results_ui');

hMIPax =gca ;
setcolormap('gray');

hMIPax = spm_mip_ui(xSPM.Z,coor',xSPM.M,xSPM.DIM,hMIPax);
spm_XYZreg('XReg',hReg,hMIPax,'spm_mip_ui');

colormap(handles.colormap);
setcolormap('gray-hot');

if isempty(handles.M)
    spm_XYZreg('SetReg',handles.structureEdit,hReg);
else
    set(handles.structureEdit, 'UserData', struct('hReg',hReg,'M',handles.M,'D', handles.DIM,'xyx',[0 0 0]));
end

spm_XYZreg('Add2Reg',hReg,handles.structureEdit,@CallBack_structureEdit);

%%%%%%%%%%%%%%% Reg end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%dawnsong, revise to allow color setting
if multiple == 0
    [hgraph, hSection, hcolorbar] = spm_sections(SPM,hReg,targetFile,handles);
else 
    [hgraph, hSection, hcolorbar] = spm_sections(cSPM,hReg,targetFile,handles);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% spm_sections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fgraph, hMe, hcolorbar] = spm_sections(SPM,hReg,targetFile,handles)
% rendering of regional effects [SPM{Z}] on orthogonal sections
% FORMAT spm_sections(SPM,hReg)
%
% SPM  - xSPM structure containing details of excursion set
% hReg - handle of MIP register
%
% see spm_getSPM for details
%_______________________________________________________________________
%
% spm_sections is called by spm_results and uses variables in SPM and
% VOL to create three orthogonal sections though a background image.
% Regional foci from the selected SPM are rendered on this image.
%
%_______________________________________________________________________
% @(#)spm_sections.m	2.14	John Ashburner 02/09/05

Fgraph = gcf;

spms = fullfile(spm('dir'),'canonical', 'single_subj_T1.mnc');
if exist('targetFile')
    spms = targetFile;
end

spm_orthviews('Reset');
global st
st.fig = Fgraph;
st.Space = spm_matrix([0 0 0  0 0 -pi/2])*st.Space;

spm_orthviews('Image',spms,handles.sectionViewPosition); % position
spm_orthviews('MaxBB');
spm_orthviews('register',hReg);
if ~iscell(SPM)0
    spm_orthviews('addblobs',1,SPM.XYZ,SPM.Z,SPM.M);
elseif length(SPM) == 1
    SPM = SPM{1};
    spm_orthviews('addblobs',1,SPM.XYZ,SPM.Z,SPM.M);
else
    %%colors = repmat([1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1],20,1);
    colors=flipud(pink(3));
    for ii=1:length(SPM)
      try
        if ii<=length(handles.hLegend),
          spm_orthviews('addcolouredblobs',1,SPM{ii}.XYZ, SPM{ii}.Z, SPM{ii}.M, get(handles.hLegend{ii},'ForeGroundColor')); %%colors(ii,:));
        else
          spm_orthviews('addcolouredblobs',1,SPM{ii}.XYZ, SPM{ii}.Z, SPM{ii}.M, colors(1+mod(ii-length(handles.hLegend), length(colors)),:)); 
        end
      catch
        spm_orthviews('addcolouredblobs',1,SPM{ii}.XYZ, SPM{ii}.Z, SPM{ii}.M, colors(ii,:)); 
      end
    end
end
spm_orthviews('Redraw');

hMe = st.registry.hMe;
try
    hcolorbar = st.vols{1}.blobs{1}.cbar;
catch
    hcolorbar = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cuixuRenderView
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = cuixuRenderView(style, mniCoord, intensity, varargin)
% function h = cuixuRenderView(style, mniCoord1, intensity1, mniCoord2, intensity2, mniCoord3, intensity3)
% This is to project your coordinate to render view
%
% style: either 'new' or 'old'
% mniCoord: the mni coordinate, N*3 matrix
% intensity: the plot intensity of the spots, 1*N matrix. The intensity
% could be t value, for example.
%
% You can input 1, 2, or 3 pairs of coordinates and intensity.
%
% h: the returned handle for the figure
%
% Xu Cui
% 2004/11/11

global M_;
global DIM_;

if nargin < 4
	dat.XYZ = mni2cor(mniCoord, M_);
	dat.XYZ = dat.XYZ';
	if nargin < 2
        dat.t = ones(size(mniCoord,1),1);
	else    
        dat.t = intensity;
	end
	dat.mat = M_;
	dat.dim = DIM_;
else
    if mod(nargin+1,2) ~=0
        disp('You should put even number of parameters.')
        return;
    end
    for ii=1:(2+length(varargin))/2
        if ii==1
    		dat(ii).XYZ = mni2cor(mniCoord, M_);
            dat(ii).t = intensity;            
        else
    		dat(ii).XYZ = mni2cor(varargin{2*(ii-1)-1}, M_);
            dat(ii).t = varargin{2*(ii-1)};
        end
    	dat(ii).XYZ = dat(ii).XYZ';        
		dat(ii).mat = M_;
		dat(ii).dim = DIM_;        
    end 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% spm_render
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(style, 'old')
    brt = nan;
else
    brt = 1;
end
h = spm_render(dat,brt,fullfile(spm('dir'), 'rend', 'render_single_subj.mat'));

function Fgraph = spm_render(dat,brt,rendfile)
% Render blobs on surface of a 'standard' brain
% FORMAT spm_render(dat,brt,rendfile)
%
% dat - a vertical cell array of length 1 to 3
%       - each element is a structure containing:
%         - XYZ - the x, y & z coordinates of the transformed t values.
%                 in units of voxels.
%         - t   - the SPM{.} values
%         - mat - affine matrix mapping from XYZ voxels to Talairach.
%         - dim - dimensions of volume from which XYZ is drawn.
% brt - brightness control:
%            If NaN, then displays using the old style with hot
%            metal for the blobs, and grey for the brain.
%            Otherwise, it is used as a ``gamma correction'' to
%            optionally brighten the blobs up a little.
% rendfile - the file containing the images to render on to. See also
%            spm_xbrain.m.
%
% Without arguments, spm_render acts as its own UI.
%_______________________________________________________________________
% 
% spm_render prompts for details of up to three SPM{Z}s or SPM{t}s that
% are then displayed superimposed on the surface of a standard brain.
%
% The first is shown in red, then green then blue.
%
% The blobs which are displayed are the integral of all transformed t
% values, exponentially decayed according to their depth. Voxels that
% are 10mm behind the surface have half the intensity of ones at the
% surface.
%_______________________________________________________________________
% @(#)spm_render.m	2.19 John Ashburner FIL 02/10/29

%-Parse arguments, get data if not passed as parameters
%=======================================================================
if nargin < 1
	SPMid = spm('FnBanner',mfilename,'2.19');
	[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Results: render',0);

	num   = spm_input('Number of sets',1,'1 set|2 sets|3 sets',[1 2 3]);

	for i = 1:num,
		[SPM,VOL] = spm_getSPM;
		dat(i)    = struct(	'XYZ',	VOL.XYZ,...
					't',	VOL.Z',...
					'mat',	VOL.M,...
					'dim',	VOL.DIM);
	end;
	showbar = 1;
else,
	num     = length(dat);
	showbar = 0;
end;

% get surface
%-----------------------------------------------------------------------
if nargin < 3,
	rendfile = spm_get(1,'render*.mat','Render file',fullfile(spm('Dir'),'rend'));
end;

% get brightness
%-----------------------------------------------------------------------
if nargin < 2,
	brt = 1;
	if num==1,
		brt = spm_input('Style',1,'new|old',[1 NaN], 1);
    end;
	if isfinite(brt),
		brt = spm_input('Brighten blobs',1,'none|slightly|more|lots',[1 0.75 0.5 0.25], 1);
    end;
end;



% Perform the rendering
%=======================================================================
spm('Pointer','Watch')

load(rendfile);

if (exist('rend') ~= 1), % Assume old format...
	rend = cell(size(Matrixes,1),1);
	for i=1:size(Matrixes,1),
		rend{i}=struct('M',eval(Matrixes(i,:)),...
			'ren',eval(Rens(i,:)),...
			'dep',eval(Depths(i,:)));
		rend{i}.ren = rend{i}.ren/max(max(rend{i}.ren));
	end;
end;

if showbar, spm_progress_bar('Init', size(dat,1)*length(rend),...
			'Formatting Renderings', 'Number completed'); end;
for i=1:length(rend),
	rend{i}.max=0;
	rend{i}.data = cell(size(dat,1),1);
	if issparse(rend{i}.ren),
		% Assume that images have been DCT compressed
		% - the SPM99 distribution was originally too big.
		d = size(rend{i}.ren);
		B1 = spm_dctmtx(d(1),d(1));
		B2 = spm_dctmtx(d(2),d(2));
		rend{i}.ren = B1*rend{i}.ren*B2';
		% the depths did not compress so well with
		% a straight DCT - therefore it was modified slightly
		rend{i}.dep = exp(B1*rend{i}.dep*B2')-1;
	end;
	msk = find(rend{i}.ren>1);rend{i}.ren(msk)=1;
	msk = find(rend{i}.ren<0);rend{i}.ren(msk)=0;
	if showbar, spm_progress_bar('Set', i); end;
end;
if showbar, spm_progress_bar('Clear'); end;

if showbar, spm_progress_bar('Init', length(dat)*length(rend),...
			'Making pictures', 'Number completed'); end;

mx = zeros(length(rend),1)+eps;
mn = zeros(length(rend),1);

for j=1:length(dat),
	XYZ = dat(j).XYZ;
	t   = dat(j).t;
	dim = dat(j).dim;
	mat = dat(j).mat;

	for i=1:length(rend),

		% transform from Taliarach space to space of the rendered image
		%-------------------------------------------------------
		M1  = rend{i}.M*dat(j).mat;
		zm  = sum(M1(1:2,1:3).^2,2).^(-1/2);
		M2  = diag([zm' 1 1]);
		M  = M2*M1;
		cor = [1 1 1 ; dim(1) 1 1 ; 1 dim(2) 1; dim(1) dim(2) 1 ;
		       1 1 dim(3) ; dim(1) 1 dim(3) ; 1 dim(2) dim(3); dim(1) dim(2) dim(3)]';
		tcor= M(1:3,1:3)*cor + M(1:3,4)*ones(1,8);
		off = min(tcor(1:2,:)');
		M2  = spm_matrix(-off+1)*M2;
		M  = M2*M1;
		xyz = (M(1:3,1:3)*XYZ + M(1:3,4)*ones(1,size(XYZ,2)));
		d2  = ceil(max(xyz(1:2,:)'));

		% calculate 'depth' of values
		%-------------------------------------------------------
		dep = spm_slice_vol(rend{i}.dep,spm_matrix([0 0 1])*inv(M2),d2,1);
		z1  = dep(round(xyz(1,:))+round(xyz(2,:)-1)*size(dep,1));

		if ~isfinite(brt), msk = find(xyz(3,:) < (z1+20) & xyz(3,:) > (z1-5));
		else,      msk = find(xyz(3,:) < (z1+60) & xyz(3,:) > (z1-5)); end;

		if ~isempty(msk),

			% generate an image of the integral of the blob values.
			%-----------------------------------------------
			xyz = xyz(:,msk);
			if ~isfinite(brt), t0  = t(msk);
			else,	dst = xyz(3,:) - z1(msk);
				dst = max(dst,0);
				t0  = t(msk).*exp((log(0.5)/10)*dst)';
			end;
			X0  = full(sparse(round(xyz(1,:)), round(xyz(2,:)), t0, d2(1), d2(2)));
			hld = 1; if ~isfinite(brt), hld = 0; end;
			X   = spm_slice_vol(X0,spm_matrix([0 0 1])*M2,size(rend{i}.dep),hld);
			msk = find(X<0);
			X(msk) = 0;
		else,
			X = zeros(size(rend{i}.dep));
		end;

		% Brighten the blobs
		if isfinite(brt), X = X.^brt; end;

		mx(j) = max([mx(j) max(max(X))]);
		mn(j) = min([mn(j) min(min(X))]);

		rend{i}.data{j} = X;

		if showbar, spm_progress_bar('Set', i+(j-1)*length(rend)); end;
	end;
end;

mxmx = max(mx);
mnmn = min(mn);

if showbar, spm_progress_bar('Clear'); end;
Fgraph = gcf;%spm_figure('GetWin','Graphics');
%spm_results_ui('Clear',Fgraph);

nrow = ceil(length(rend)/2);
hght = 0.25;
width = 0.25;
x0 = 0.5;
y0 = 0.01;
% subplot('Position',[0, 0, 1, hght]);
ax=axes('Parent',Fgraph,'units','normalized','Position',[0, 0, 0.5, hght],'Visible','off');
%ax=axes;
%image(0,'Parent',ax);
set(ax,'YTick',[],'XTick',[]);

if ~isfinite(brt),
	% Old style split colourmap display.
	%---------------------------------------------------------------
	load Split;
	colormap(split);
	for i=1:length(rend),
		ren = rend{i}.ren;
		X   = (rend{i}.data{1}-mnmn)/(mxmx-mnmn);
		msk = find(X);
		ren(msk) = X(msk)+(1+1.51/64);
		ax=axes('Parent',Fgraph,'units','normalized',...
			'Position',[x0+rem(i-1,2)*width, y0+floor((i-1)/2)*hght/nrow*2, width, hght/nrow*2],...
			'Visible','off');
		image(ren*64,'Parent',ax);
		set(ax,'DataAspectRatio',[1 1 1], ...
			'PlotBoxAspectRatioMode','auto',...
			'YTick',[],'XTick',[],'XDir','normal','YDir','normal');
	end;
else,
	% Combine the brain surface renderings with the blobs, and display using
	% 24 bit colour.
	%---------------------------------------------------------------
	for i=1:length(rend),
		ren = rend{i}.ren;
		X = cell(3,1);
		for j=1:length(rend{i}.data),
			X{j} = rend{i}.data{j}/(mxmx-mnmn)-mnmn;
		end
		for j=(length(rend{i}.data)+1):3
			X{j}=zeros(size(X{1}));
		end

		rgb = zeros([size(ren) 3]);
		tmp = ren.*max(1-X{1}-X{2}-X{3},0);
		rgb(:,:,1) = tmp + X{1};
		rgb(:,:,2) = tmp + X{2};
		rgb(:,:,3) = tmp + X{3};

		ax=axes('Parent',Fgraph,'units','normalized',...
			'Position',[x0+rem(i-1,2)*width, y0+floor((i-1)/2)*hght/nrow*2, width, hght/nrow*2],...
			'Visible','off');
		image(rgb,'Parent',ax);
		set(ax,'DataAspectRatio',[1 1 1], ...
			'PlotBoxAspectRatioMode','auto',...
			'YTick',[],'XTick',[],...
			'XDir','normal','YDir','normal');
	end;
end;

spm('Pointer')
return;

function [P,p,Em,En,EN] = spm_P(c,k,Z,df,STAT,R,n,S)
% Returns the [un]corrected P value using unifed EC theory
% FORMAT [P p Em En EN] = spm_P(c,k,Z,df,STAT,R,n,S)
%
% c     - cluster number 
% k     - extent {RESELS}
% Z     - height {minimum over n values}
% df    - [df{interest} df{error}]
% STAT  - Statistical field
%       'Z' - Gaussian field
%       'T' - T - field
%       'X' - Chi squared field
%       'F' - F - field
%       'P' - Posterior probability
% R     - RESEL Count {defining search volume}
% n     - number of component SPMs in conjunction
% S     - Voxel count
%
% P     - corrected   P value  - P(n > kmax}
% p     - uncorrected P value  - P(n > k}
% Em    - expected total number of maxima {m}
% En    - expected total number of resels per cluster {n}
% EN    - expected total number of voxels {N}
%
%__________________________________________________________________________
%
% spm_P determines corrected and uncorrected p values, using the minimum
% of different valid methods. 
%
% See the individual methods for details
%
%     spm_P_RF
%     spm_P_Bonf
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Thomas Nichols
% $Id: spm_P.m 2690 2009-02-04 21:44:28Z guillaume $


% set global var NOBONF to 1 to turn off Bonferroni
%--------------------------------------------------------------------------
global NOBONF; if ~isempty(NOBONF) && NOBONF, S = []; end

if (nargin < 8), S = []; end

[P,p,Em,En,EN] = spm_P_RF(c,k,Z,df,STAT,R,n);

% Use lower Bonferroni P value (if possible)
%==========================================================================
if ~isempty(S) && (c == 1 && k == 0) && ~(length(R) == 1 && R == 1)
    P = min(P,spm_P_Bonf(Z,df,STAT,S,n));
end

function [P,p,Em,En,EN] = spm_P_RF(c,k,Z,df,STAT,R,n)
% Returns the [un]corrected P value using unifed EC theory
% FORMAT [P p Em En EN] = spm_P_RF(c,k,Z,df,STAT,R,n)
%
% c     - cluster number 
% k     - extent {RESELS}
% Z     - height {minimum over n values}
% df    - [df{interest} df{error}]
% STAT  - Statistical field
%       'Z' - Gaussian field
%       'T' - T - field
%       'X' - Chi squared field
%       'F' - F - field
% R     - RESEL Count {defining search volume}
% n     - number of component SPMs in conjunction
%
% P     - corrected   P value  - P(n > kmax}
% p     - uncorrected P value  - P(n > k}
% Em    - expected total number of maxima {m}
% En    - expected total number of resels per cluster {n}
% EN    - expected total number of voxels {N}
%
%___________________________________________________________________________
%
% spm_P_RF returns the probability of c or more clusters with more than
% k voxels in volume process of R RESELS thresholded at u.  All p values
% can be considered special cases:
%
% spm_P_RF(1,0,Z,df,STAT,1,n) = uncorrected p value
% spm_P_RF(1,0,Z,df,STAT,R,n) = corrected p value {based on height Z)
% spm_P_RF(1,k,u,df,STAT,R,n) = corrected p value {based on extent k at u)
% spm_P_RF(c,k,u,df,STAT,R,n) = corrected p value {based on number c at k and u)
% spm_P_RF(c,0,u,df,STAT,R,n) = omnibus   p value {based on number c at u)
%
% If n > 1 a conjunction probility over the n values of the statistic
% is returned
%
% Ref: Hasofer AM (1978) Upcrossings of random fields
% Suppl Adv Appl Prob 10:14-21
% Ref: Friston et al (1993) Comparing functional images: Assessing
% the spatial extent of activation foci
% Ref: Worsley KJ et al 1996, Hum Brain Mapp. 4:58-73
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_P_RF.m 2691 2009-02-04 21:46:04Z guillaume $



% get expectations
%===========================================================================

% get EC densities
%---------------------------------------------------------------------------
D       = find(R, 1, 'last' );
R       = R(1:D);
G       = sqrt(pi)./gamma(([1:D])/2);
EC      = spm_ECdensity(STAT,Z,df);
EC      = EC([1:D]) + eps;

% corrected p value
%---------------------------------------------------------------------------
P       = triu(toeplitz(EC'.*G))^n;
P       = P(1,:);
EM      = (R./G).*P;            % <maxima> over D dimensions
Em      = sum(EM);          % <maxima>
EN      = P(1)*R(D);            % <voxels>
En      = EN/EM(D);         % En = EN/EM(D);

% get P{n > k}
%===========================================================================

% assume a Gaussian form for P{n > k} ~ exp(-beta*k^(2/D))
% Appropriate for SPM{Z} and high d.f. SPM{T}
%---------------------------------------------------------------------------
D       = D - 1;
if     ~k || ~D

    p    = 1;

elseif STAT == 'Z'

    beta = (gamma(D/2 + 1)/En)^(2/D);
    p    = exp(-beta*(k^(2/D)));

elseif STAT == 'T'

    beta = (gamma(D/2 + 1)/En)^(2/D);
    p    = exp(-beta*(k^(2/D)));

elseif STAT == 'X'

    beta = (gamma(D/2 + 1)/En)^(2/D);
    p    = exp(-beta*(k^(2/D)));

elseif STAT == 'F'

    beta = (gamma(D/2 + 1)/En)^(2/D);
    p    = exp(-beta*(k^(2/D)));

end

% Poisson clumping heuristic {for multiple clusters}
%===========================================================================
P       = 1 - spm_Pcdf(c - 1,(Em + eps)*p);


% set P and p = [] for non-implemented cases
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if k > 0 && n > 1
    P    = []; p = [];
end
if k > 0 && (STAT == 'X' || STAT == 'F')
    P    = []; p = [];
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




% spm_ECdensity
%===========================================================================

function [EC] = spm_ECdensity(STAT,t,df)
% Returns the EC density
%___________________________________________________________________________
%
% Reference : Worsley KJ et al 1996, Hum Brain Mapp. 4:58-73
%
%---------------------------------------------------------------------------

% EC densities (EC}
%---------------------------------------------------------------------------
t     = t(:)';
if      STAT == 'Z'

    % Gaussian Field
    %-------------------------------------------------------------------
    a       = 4*log(2);
    b       = exp(-t.^2/2);

    EC(1,:) = 1 - spm_Ncdf(t);
    EC(2,:) = a^(1/2)/(2*pi)*b;
    EC(3,:) = a/((2*pi)^(3/2))*b.*t;
    EC(4,:) = a^(3/2)/((2*pi)^2)*b.*(t.^2 - 1);

elseif  STAT == 'T'

    % T - Field
    %-------------------------------------------------------------------
    v       = df(2);
    a       = 4*log(2);
    b       = exp(gammaln((v+1)/2) - gammaln(v/2));
    c       = (1+t.^2/v).^((1-v)/2);

    EC(1,:) = 1 - spm_Tcdf(t,v);
    EC(2,:) = a^(1/2)/(2*pi)*c;
    EC(3,:) = a/((2*pi)^(3/2))*c.*t/((v/2)^(1/2))*b;
    EC(4,:) = a^(3/2)/((2*pi)^2)*c.*((v-1)*(t.^2)/v - 1);

elseif  STAT == 'X'

    % X - Field
    %-------------------------------------------------------------------
    v       = df(2);
    a       = (4*log(2))/(2*pi);
    b       = t.^(1/2*(v - 1)).*exp(-t/2-gammaln(v/2))/2^((v-2)/2);

    EC(1,:) = 1 - spm_Xcdf(t,v);
    EC(2,:) = a^(1/2)*b;
    EC(3,:) = a*b.*(t-(v-1));
    EC(4,:) = a^(3/2)*b.*(t.^2-(2*v-1)*t+(v-1)*(v-2));

elseif  STAT == 'F'

    % F Field
    %-------------------------------------------------------------------
    k       = df(1);
    v       = df(2);
    a       = (4*log(2))/(2*pi);
    b       = gammaln(v/2) + gammaln(k/2);

    EC(1,:) = 1 - spm_Fcdf(t,df);
    EC(2,:) = a^(1/2)*exp(gammaln((v+k-1)/2)-b)*2^(1/2)...
        *(k*t/v).^(1/2*(k-1)).*(1+k*t/v).^(-1/2*(v+k-2));
    EC(3,:) = a*exp(gammaln((v+k-2)/2)-b)*(k*t/v).^(1/2*(k-2))...
            .*(1+k*t/v).^(-1/2*(v+k-2)).*((v-1)*k*t/v-(k-1));
    EC(4,:) = a^(3/2)*exp(gammaln((v+k-3)/2)-b)...
        *2^(-1/2)*(k*t/v).^(1/2*(k-3)).*(1+k*t/v).^(-1/2*(v+k-2))...
            .*((v-1)*(v-2)*(k*t/v).^2-(2*v*k-v-k-1)*(k*t/v)+(k-1)*(k-2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% spm_list (copied from spm 5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = spm_list(varargin)
% Display and analysis of SPM{.}
% FORMAT TabDat = spm_list('List',SPM,hReg,[Num,Dis,Str])
% Summary list of local maxima for entire volume of interest
% FORMAT TabDat = spm_list('ListCluster',SPM,hReg,[Num,Dis,Str])
% List of local maxima for a single suprathreshold cluster
%
% SPM    - structure containing SPM, distribution & filtering details
%        - required fields are:
% .swd   - SPM working directory - directory containing current SPM.mat
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests        
% .STAT  - distribution {Z, T, X or F}     
% .df    - degrees of freedom [df{interest}, df{residual}]
% .u     - height threshold
% .k     - extent threshold {voxels}
% .XYZ   - location of voxels {voxel coords}
% .XYZmm - location of voxels {mm}
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}     
% .M     - voxels - > mm matrix
% .VOX   - voxel dimensions {mm}
% .Vspm  - mapped statistic image(s)
% .Ps    - uncorrected P values in searched volume (for FDR)
%
% (see spm_getSPM for further details of xSPM structures)
%
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Num    - number of maxima per cluster
% Dis    - distance among clusters (mm)
% Str    - header string
%
% TabDat - Structure containing table data
%        - fields are
% .tit   - table Title (string)
% .hdr   - table header (2x11 cell array)
% .fmt   - fprintf format strings for table data (1x11 cell array)
% .str   - table filtering note (string)
% .ftr   - table footnote information (4x2 cell array)
% .dat   - table data (Nx11 cell array)
%
%                           ----------------
%
% FORMAT spm_list('TxtList',TabDat,c)
% Prints a tab-delimited text version of the table
% TabDat - Structure containing table data (format as above)
% c      - Column of table data to start text table at
%          (E.g. c=3 doesn't print set-level results contained in columns 1 & 2)
%                           ----------------
%
% FORMAT spm_list('SetCoords',xyz,hAx,hC)
% Highlighting of table co-ordinates (used by results section registry)
% xyz    - 3-vector of new co-ordinate
% hAx    - table axis (the registry object for tables)
% hReg   - Handle of caller (not used)
%_______________________________________________________________________
%
% spm_list characterizes SPMs (thresholded at u and k) in terms of
% excursion sets (a collection of face, edge and vertex connected
% subsets or clusters).  The currected significance of the results are
% based on set, cluster and voxel-level inferences using distributional
% approximations from the Theory of Gaussian Fields.  These
% distributions assume that the SPM is a reasonable lattice
% approximation of a continuous random field with known component field
% smoothness.
%
% The p values are based on the probability of obtaining c, or more,
% clusters of k, or more, resels above u, in the volume S analysed =
% P(u,k,c).  For specified thresholds u, k, the set-level inference is
% based on the observed number of clusters C, = P(u,k,C).  For each
% cluster of size K the cluster-level inference is based on P(u,K,1)
% and for each voxel (or selected maxima) of height U, in that cluster,
% the voxel-level inference is based on P(U,0,1).  All three levels of
% inference are supported with a tabular presentation of the p values
% and the underlying statistic:
%
% Set-level     - c    = number of suprathreshold clusters
%               - P    = prob(c or more clusters in the search volume)
%
% Cluster-level - k    = number of voxels in this cluster
%               - Pc   = prob(k or more voxels in the search volume)
%               - Pu   = prob(k or more voxels in a cluster)
%
% Voxel-level   - T/F  = Statistic upon which the SPM is based
%               - Ze   = The eqivalent Z score - prob(Z > Ze) = prob(t > T)
%               - Pc   = prob(Ze or higher in the search volume)
%               - Qu   = Expd(Prop of false positives among voxels >= Ze)
%               - Pu   = prob(Ze or higher at that voxel)
%
% x,y,z (mm)    - Coordinates of the voxel
%
% The table is grouped by regions and sorted on the Ze-variate of the
% primary maxima.  Ze-variates (based on the uncorrected p value) are the
% Z score equivalent of the statistic. Volumes are expressed in voxels.
%
% Clicking on values in the table returns the value to the Matlab
% workspace. In addition, clicking on the co-ordinates jumps the
% results section cursor to that location. The table has a context menu
% (obtained by right-clicking in the background of the table),
% providing options to print the current table as a text table, or to
% extract the table data to the Matlab workspace.
%
%_______________________________________________________________________
% @(#)spm_list.m	2.43 Karl Friston, Andrew Holmes 02/10/31


% satellite figure global variable
%-----------------------------------------------------------------------
global SatWindow

%=======================================================================
switch lower(varargin{1}), case 'list'                            %-List
%=======================================================================
% FORMAT TabDat = spm_list('list',SPM,hReg)

%-Tolerance for p-value underflow, when computing equivalent Z's
%-----------------------------------------------------------------------
tol = eps*10;

%-Parse arguments and set maxima number and separation
%-----------------------------------------------------------------------
if nargin < 2,	error('insufficient arguments'),     end
if nargin < 3,	hReg = []; else, hReg = varargin{3}; end


%-Get current location (to highlight selected voxel in table)
%-----------------------------------------------------------------------
%xyzmm     = spm_results_ui('GetCoords');
xyzmm = [0 0 0]';

%-Extract data from xSPM
%-----------------------------------------------------------------------
S     = varargin{2}.S;
R     = varargin{2}.R;
FWHM  = varargin{2}.FWHM;
VOX   = varargin{2}.VOX;
n     = varargin{2}.n;
STAT  = varargin{2}.STAT;
df    = varargin{2}.df;
u     = varargin{2}.u;
M     = varargin{2}.M;
v2r   = 1/prod(FWHM(~isinf(FWHM)));			%-voxels to resels
k     = varargin{2}.k*v2r;
QPs   = varargin{2}.Ps;					% Needed for FDR
QPs   = sort(QPs(:));

%-get number and separation for maxima to be reported
%-----------------------------------------------------------------------
if length(varargin) > 3

	Num    = varargin{4};		% number of maxima per cluster
	Dis    = varargin{5};		% distance among clusters (mm)
else
	Num    = 3;
	Dis    = 8;
end
if length(varargin) > 5

	Title  = varargin{6};
else
	Title  = 'p-values adjusted for search volume';
end


%-Setup graphics panel
%-----------------------------------------------------------------------
spm('Pointer','Watch')
if SatWindow
	Fgraph = SatWindow;
	figure(Fgraph);
else
	Fgraph = figure('unit','normalized','position',[0.4,0.1,0.55,0.5],'Color','w',...
        'Name','volume', 'NumberTitle','off','resize','off','MenuBar','none');
    Fgraph = gcf;
end
%spm_results_ui('Clear',Fgraph)
FS    = spm('FontSizes');			%-Scaled font sizes
PF    = spm_platform('fonts');			%-Font names (for this platform)


%-Table header & footer
%=======================================================================

%-Table axes & Title
%----------------------------------------------------------------------
if SatWindow, ht = 0.85; bot = .14; else, ht = 0.8; bot = 0.15; end;

if STAT == 'P'
	Title = 'Posterior Probabilities';
end
	
hAx   = axes('Position',[0.025 bot 0.9 ht],...
	'DefaultTextFontSize',FS(8),...
	'DefaultTextInterpreter','Tex',...
	'DefaultTextVerticalAlignment','Baseline',...
	'Units','points',...
	'Visible','off');

AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)])
dy    = FS(9);
y     = floor(AxPos(4)) - dy;

text(0,y,['Statistics:  \it\fontsize{',num2str(FS(9)),'}',Title],...
	'FontSize',FS(11),'FontWeight','Bold');	y = y - dy/2;
line([0 1],[y y],'LineWidth',3,'Color','r'),	y = y - 9*dy/8;


%-Construct table header
%-----------------------------------------------------------------------
set(gca,'DefaultTextFontName',PF.helvetica,'DefaultTextFontSize',FS(8))

Hc = [];
Hp = [];
h  = text(0.01,y,	'set-level','FontSize',FS(9));		Hc = [Hc,h];
h  = line([0,0.11],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');	Hc = [Hc,h];
h  = text(0.08,y-9*dy/8,	'\itc ');			Hc = [Hc,h];
h  = text(0.02,y-9*dy/8,	'\itp ');			Hc = [Hc,h];
								Hp = [Hp,h];
text(0.22,y,		'cluster-level','FontSize',FS(9));
line([0.15,0.41],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');
h  = text(0.16,y-9*dy/8,	'\itp \rm_{corrected}');	Hp = [Hp,h];
h  = text(0.33,y-9*dy/8,	'\itp \rm_{uncorrected}');	Hp = [Hp,h];
h  = text(0.26,y-9*dy/8,	'\itk \rm_E');

text(0.60,y,		'voxel-level','FontSize',FS(9));
line([0.46,0.86],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');
h  = text(0.46,y-9*dy/8,	'\itp \rm_{FWE-corr}');		Hp = [Hp,h];
h  = text(0.55,y-9*dy/8,        '\itp \rm_{FDR-corr}');		Hp = [Hp,h];
h  = text(0.79,y-9*dy/8,	'\itp \rm_{uncorrected}');	Hp = [Hp,h];
h  = text(0.64,y-9*dy/8,	 sprintf('\\it%c',STAT));
h  = text(0.72,y-9*dy/8,	'(\itZ\rm_\equiv)');

text(0.93,y - dy/2,['x,y,z \fontsize{',num2str(FS(8)),'}\{mm\}']);


%-Headers for text table...
%-----------------------------------------------------------------------
TabDat.tit = Title;
TabDat.hdr = {	'set',		'c';...
		'set',		'p';...
		'cluster',	'p(cor)';...
		'cluster',	'equivk';...
		'cluster',	'p(unc)';...
		'voxel',	'p(FWE-cor)';...
		'voxel',	'p(FDR-cor)';...
		'voxel',	 STAT;...
		'voxel',	'equivZ';...
		'voxel',	'p(unc)';...
		'',		'x,y,z {mm}'}';...
		
TabDat.fmt = {	'%-0.3f','%g',...				%-Set
		'%0.3f', '%0.0f', '%0.3f',...			%-Cluster
		'%0.3f', '%0.3f', '%6.2f', '%5.2f', '%0.3f',...	%-Voxel
		'%3.0f %3.0f %3.0f'};				%-XYZ

%-Column Locations
%-----------------------------------------------------------------------
tCol       = [  0.00      0.07 ...				%-Set
	        0.16      0.26      0.34 ...			%-Cluster
	        0.46      0.55      0.62      0.71      0.80 ...%-Voxel
                0.92];						%-XYZ

% move to next vertial postion marker
%-----------------------------------------------------------------------
y     = y - 7*dy/4;
line([0 1],[y y],'LineWidth',1,'Color','r')
y     = y - 5*dy/4;
y0    = y;


%-Table filtering note
%-----------------------------------------------------------------------
if isinf(Num)
	TabDat.str = sprintf('table shows all local maxima > %.1fmm apart',Dis);
else
	TabDat.str = sprintf(['table shows %d local maxima ',...
		'more than %.1fmm apart'],Num,Dis);
end
text(0.5,4,TabDat.str,'HorizontalAlignment','Center','FontName',PF.helvetica,...
    'FontSize',FS(8),'FontAngle','Italic')


%-Volume, resels and smoothness (if classical inference)
%-----------------------------------------------------------------------
line([0 1],[0 0],'LineWidth',1,'Color','r')
if STAT ~= 'P'
%-----------------------------------------------------------------------
FWHMmm          = FWHM.*VOX; 				% FWHM {mm}
Pz              = spm_P(1,0,u,df,STAT,1,n,S);
Pu              = spm_P(1,0,u,df,STAT,R,n,S);
Qu              = spm_P_FDR(u,df,STAT,n,QPs);
[P Pn Em En EN] = spm_P(1,k,u,df,STAT,R,n,S);


%-Footnote with SPM parameters
%-----------------------------------------------------------------------
set(gca,'DefaultTextFontName',PF.helvetica,...
	'DefaultTextInterpreter','None','DefaultTextFontSize',FS(8))
TabDat.ftr    = cell(5,2);
TabDat.ftr{1} = ...
	sprintf('Height threshold: %c = %0.2f, p = %0.3f (%0.3f)',...
		 STAT,u,Pz,Pu);
TabDat.ftr{2} = ...
	sprintf('Extent threshold: k = %0.0f voxels, p = %0.3f (%0.3f)',...
	         k/v2r,Pn,P);
TabDat.ftr{3} = ...
	sprintf('Expected voxels per cluster, <k> = %0.3f',En/v2r);
TabDat.ftr{4} = ...
	sprintf('Expected number of clusters, <c> = %0.2f',Em*Pn);
TabDat.ftr{5} = ...
	sprintf('Expected false discovery rate, <= %0.2f',Qu);
TabDat.ftr{6} = ...
	sprintf('Degrees of freedom = [%0.1f, %0.1f]',df);
TabDat.ftr{7} = ...
	sprintf(['Smoothness FWHM = %0.1f %0.1f %0.1f {mm} ',...
		 ' = %0.1f %0.1f %0.1f {voxels}'],FWHMmm,FWHM);
TabDat.ftr{8} = ...
	sprintf('Search vol: %0.0f cmm; %0.0f voxels; %0.1f resels',S*prod(VOX),S,R(end));
TabDat.ftr{9} = ...
	sprintf(['Voxel size: [%0.1f, %0.1f, %0.1f] mm ',...
		' (1 resel = %0.2f voxels)'],VOX,prod(FWHM));

text(0.0,-1*dy,TabDat.ftr{1},...
	'UserData',[u,Pz,Pu,Qu],'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.0,-2*dy,TabDat.ftr{2},...
	'UserData',[k/v2r,Pn,P],'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.0,-3*dy,TabDat.ftr{3},...
	'UserData',En/v2r,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.0,-4*dy,TabDat.ftr{4},...
	'UserData',Em*Pn,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.0,-5*dy,TabDat.ftr{5},...
	'UserData',Qu,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.5,-1*dy,TabDat.ftr{6},...
	'UserData',df,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.5,-2*dy,TabDat.ftr{7},...
	'UserData',FWHMmm,'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.5,-3*dy,TabDat.ftr{8},...
	'UserData',[S*prod(VOX),S,R(end)],...
	'ButtonDownFcn','get(gcbo,''UserData'')')
text(0.5,-4*dy,TabDat.ftr{9},...
	'UserData',[VOX,prod(FWHM)],...
	'ButtonDownFcn','get(gcbo,''UserData'')')

end % Classical


%-Characterize excursion set in terms of maxima
% (sorted on Z values and grouped by regions)
%=======================================================================
if ~length(varargin{2}.Z)
	text(0.5,y-6*dy,'no suprathreshold clusters',...
		'HorizontalAlignment','Center',...
		'FontAngle','Italic','FontWeight','Bold',...
		'FontSize',FS(16),'Color',[1,1,1]*.5);
	TabDat.dat = cell(0,11);
	varargout  = {TabDat};
	spm('Pointer','Arrow')
	return
end

% Includes Darren Gitelman's code for working around
% spm_max for conjunctions with negative thresholds
%-----------------------------------------------------------------------
minz        = abs(min(min(varargin{2}.Z)));
zscores     = 1 + minz + varargin{2}.Z;
[N Z XYZ A] = spm_max(zscores,varargin{2}.XYZ);
Z           = Z - minz - 1;

%-Convert cluster sizes from voxels to resels
%-----------------------------------------------------------------------
if isfield(varargin{2},'VRvp')
	V2R = spm_get_data(varargin{2}.VRvp,XYZ);
else
	V2R = v2r;
end
N           = N.*V2R;

%-Convert maxima locations from voxels to mm
%-----------------------------------------------------------------------
XYZmm = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];



%-Table proper (& note all data in cell array)
%=======================================================================

%-Pagination variables
%-----------------------------------------------------------------------
hPage = [];
set(gca,'DefaultTextFontName',PF.courier,'DefaultTextFontSize',FS(7))


%-Set-level p values {c} - do not display if reporting a single cluster
%-----------------------------------------------------------------------
c     = max(A);					%-Number of clusters
if STAT ~= 'P'
	Pc    = spm_P(c,k,u,df,STAT,R,n,S);	%-Set-level p-value
else
	Pc    = [];
	set(Hp,'Visible','off')
end

if c > 1;
	h     = text(tCol(1),y,sprintf(TabDat.fmt{1},Pc),'FontWeight','Bold',...
		'UserData',Pc,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];
	h     = text(tCol(2),y,sprintf(TabDat.fmt{2},c),'FontWeight','Bold',...
		'UserData',c,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];
else
	set(Hc,'Visible','off')
end

TabDat.dat = {Pc,c};				%-Table data
TabLin     = 1;					%-Table data line


%-Local maxima p-values & statistics
%-----------------------------------------------------------------------
HlistXYZ = [];
while prod(size(find(isfinite(Z))))

	% Paginate if necessary
	%---------------------------------------------------------------
	if y < min(Num + 1,3)*dy

		% added Fgraph term to paginate on Satellite window
		%-------------------------------------------------------
		h     = text(0.5,-5*dy,...
			sprintf('Page %d',spm_figure('#page',Fgraph)),...
			'FontName',PF.helvetica,'FontAngle','Italic',...
			'FontSize',FS(8));

		spm_figure('NewPage',[hPage,h])
		hPage = [];
		y     = y0;
	end

    	%-Find largest remaining local maximum
    	%---------------------------------------------------------------
	[U,i]   = max(Z);			% largest maxima
	j       = find(A == A(i));		% maxima in cluster


    	%-Compute cluster {k} and voxel-level {u} p values for this cluster
    	%---------------------------------------------------------------
	Nv      = N(i)/v2r;			% extent        {voxels}


	if STAT ~= 'P'
	Pz      = spm_P(1,0,   U,df,STAT,1,n,S);% uncorrected p value
	Pu      = spm_P(1,0,   U,df,STAT,R,n,S);% FWE-corrected {based on Z)
	Qu      = spm_P_FDR(   U,df,STAT,n,QPs);% FDR-corrected {based on Z)
	[Pk Pn] = spm_P(1,N(i),u,df,STAT,R,n,S);% [un]corrected {based on k)

	if Pz < tol				% Equivalent Z-variate
	    Ze  = Inf;	 			% (underflow => can't compute)
	else
	    Ze  = spm_invNcdf(1 - Pz);
	end
	else
		Pz	= [];
		Pu      = [];
		Qu      = [];
		Pk	= [];
		Pn	= [];
		Ze      = spm_invNcdf(U);
	end


	%-Print cluster and maximum voxel-level p values {Z}
    	%---------------------------------------------------------------
	h     = text(tCol(3),y,sprintf(TabDat.fmt{3},Pk),'FontWeight','Bold',...
		'UserData',Pk,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];

	h     = text(tCol(4),y,sprintf(TabDat.fmt{4},Nv),'FontWeight','Bold',...
		'UserData',Nv,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];
	h     = text(tCol(5),y,sprintf(TabDat.fmt{5},Pn),'FontWeight','Bold',...
		'UserData',Pn,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];

	h     = text(tCol(6),y,sprintf(TabDat.fmt{6},Pu),'FontWeight','Bold',...
		'UserData',Pu,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];
	h     = text(tCol(7),y,sprintf(TabDat.fmt{7},Qu),'FontWeight','Bold',...
		'UserData',Qu,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];
	h     = text(tCol(8),y,sprintf(TabDat.fmt{8},U),'FontWeight','Bold',...
		'UserData',U,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];
	h     = text(tCol(9),y,sprintf(TabDat.fmt{9},Ze),'FontWeight','Bold',...
		'UserData',Ze,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];
	h     = ...
	text(tCol(10),y,sprintf(TabDat.fmt{10},Pz),'FontWeight','Bold',...
		'UserData',Pz,'ButtonDownFcn','get(gcbo,''UserData'')');
	hPage = [hPage, h];

	% Specifically changed so it properly finds hMIPax
	%---------------------------------------------------------------------
	h     = text(tCol(11),y,sprintf(TabDat.fmt{11},XYZmm(:,i)),...
		'FontWeight','Bold',...
		'Tag','ListXYZ',...
		'ButtonDownFcn',[...
		'hMIPax = findobj(''tag'',''hMIPax'');',...
		'spm_mip_ui(''SetCoords'',get(gcbo,''UserData''),hMIPax);'],...
		'Interruptible','off','BusyAction','Cancel',...
		'UserData',XYZmm(:,i));

	HlistXYZ = [HlistXYZ, h];
	if spm_XYZreg('Edist',xyzmm,XYZmm(:,i))<tol & ~isempty(hReg)
		set(h,'Color','r')
	end
	hPage  = [hPage, h];
 
	y      = y - dy;
	
	[TabDat.dat{TabLin,3:11}] = deal(Pk,Nv,Pn,Pu,Qu,U,Ze,Pz,XYZmm(:,i));
	TabLin = TabLin + 1;

	%-Print Num secondary maxima (> Dis mm apart)
    	%---------------------------------------------------------------
	[l q] = sort(-Z(j));				% sort on Z value
	D     = i;
	for i = 1:length(q)
	    d = j(q(i));
	    if min(sqrt(sum((XYZmm(:,D)-XYZmm(:,d)*ones(1,size(D,2))).^2)))>Dis;

		if length(D) < Num
			
			% Paginate if necessary
			%-----------------------------------------------
			if y < dy	
				h = text(0.5,-5*dy,sprintf('Page %d',...
					spm_figure('#page',Fgraph)),...
					'FontName',PF.helvetica,...
					'FontAngle','Italic',...
					'FontSize',FS(8));

				spm_figure('NewPage',[hPage,h])
				hPage = [];
				y     = y0;
			end

			% voxel-level p values {Z}
			%-----------------------------------------------
			if STAT ~= 'P'
				Pz    = spm_P(1,0,Z(d),df,STAT,1,n,S);
				Pu    = spm_P(1,0,Z(d),df,STAT,R,n,S);
				Qu    = spm_P_FDR(Z(d),df,STAT,n,QPs);
				if Pz < tol
					Ze = Inf;
				else,   Ze = spm_invNcdf(1 - Pz); end
			else
				Pz    = [];
				Pu    = [];
				Qu    = [];
				Ze    = spm_invNcdf(Z(d));
			end

			h     = text(tCol(6),y,sprintf(TabDat.fmt{6},Pu),...
				'UserData',Pu,...
				'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];

			h     = text(tCol(7),y,sprintf(TabDat.fmt{7},Qu),...
				'UserData',Qu,...
				'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];
			h     = text(tCol(8),y,sprintf(TabDat.fmt{8},Z(d)),...
				'UserData',Z(d),...
				'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];
			h     = text(tCol(9),y,sprintf(TabDat.fmt{9},Ze),...
				'UserData',Ze,...
				'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];
			h     = text(tCol(10),y,sprintf(TabDat.fmt{10},Pz),...
				'UserData',Pz,...
				'ButtonDownFcn','get(gcbo,''UserData'')');
			hPage = [hPage, h];

			% specifically modified line to use hMIPax
			%-----------------------------------------------
			h     = text(tCol(11),y,...
				sprintf(TabDat.fmt{11},XYZmm(:,d)),...
				'Tag','ListXYZ',...
				'ButtonDownFcn',[...
				    'hMIPax = findobj(''tag'',''hMIPax'');',...
				    'spm_mip_ui(''SetCoords'',',...
				    'get(gcbo,''UserData''),hMIPax);'],...
				'Interruptible','off','BusyAction','Cancel',...
				'UserData',XYZmm(:,d));

			HlistXYZ = [HlistXYZ, h];
			if spm_XYZreg('Edist',xyzmm,XYZmm(:,d))<tol & ...
				~isempty(hReg)
				set(h,'Color','r')
			end
			hPage = [hPage, h];
			D     = [D d];
			y     = y - dy;
			[TabDat.dat{TabLin,6:11}] = ...
				deal(Pu,Qu,Z(d),Ze,Pz,XYZmm(:,d));
			TabLin = TabLin+1;
		end
	    end
	end
	Z(j) = NaN;		% Set local maxima to NaN
end				% end region


%-Number and register last page (if paginated)
%-Changed to use Fgraph for numbering
%-----------------------------------------------------------------------
if spm_figure('#page',Fgraph)>1
	h = text(0.5,-5*dy,sprintf('Page %d/%d',spm_figure('#page',Fgraph)*[1,1]),...
		'FontName',PF.helvetica,'FontSize',FS(8),'FontAngle','Italic');
	spm_figure('NewPage',[hPage,h])
end

%-End: Store TabDat in UserData of axes & reset pointer
%=======================================================================
h      = uicontextmenu('Tag','TabDat',...
		'UserData',TabDat);
set(gca,'UIContextMenu',h,...
	'Visible','on',...
	'XColor','w','YColor','w')
uimenu(h,'Label','Table')
uimenu(h,'Separator','on','Label','Print text table',...
	'Tag','TD_TxtTab',...
	'CallBack',...
	'spm_list(''txtlist'',get(get(gcbo,''Parent''),''UserData''),3)',...
	'Interruptible','off','BusyAction','Cancel');
uimenu(h,'Separator','off','Label','Extract table data structure',...
	'Tag','TD_Xdat',...
	'CallBack','get(get(gcbo,''Parent''),''UserData'')',...
	'Interruptible','off','BusyAction','Cancel');
uimenu(h,'Separator','on','Label','help',...
	'Tag','TD_Xdat',...
	'CallBack','spm_help(''spm_list'')',...
	'Interruptible','off','BusyAction','Cancel');

%-Setup registry
%-----------------------------------------------------------------------
set(hAx,'UserData',struct('hReg',hReg,'HlistXYZ',HlistXYZ))
spm_XYZreg('Add2Reg',hReg,hAx,'spm_list');

%-Return TabDat structure & reset pointer
%-----------------------------------------------------------------------
varargout = {TabDat};
spm('Pointer','Arrow')





%=======================================================================
case 'listcluster'                       %-List for current cluster only
%=======================================================================
% FORMAT TabDat = spm_list('listcluster',SPM,hReg)

spm('Pointer','Watch')

%-Parse arguments
%-----------------------------------------------------------------------
if nargin < 2,	error('insufficient arguments'),     end
if nargin < 3,	hReg = []; else, hReg = varargin{3}; end
SPM    = varargin{2};

%-get number and separation for maxima to be reported
%-----------------------------------------------------------------------
if length(varargin) > 3

	Num    = varargin{4};		% number of maxima per cluster
	Dis    = varargin{5};		% distance among clusters (mm)
else
	Num    = 32;
	Dis    = 4;
end


%-if there are suprathreshold voxels, filter out all but current cluster
%-----------------------------------------------------------------------
if length(SPM.Z)

	%-Jump to voxel nearest current location
	%--------------------------------------------------------------
	[xyzmm,i] = spm_XYZreg('NearestXYZ',...
			spm_results_ui('GetCoords'),SPM.XYZmm);
	spm_results_ui('SetCoords',SPM.XYZmm(:,i));
	
	%-Find selected cluster
	%--------------------------------------------------------------
	A         = spm_clusters(SPM.XYZ);
	j         = find(A == A(i));
	SPM.Z     = SPM.Z(j);
	SPM.XYZ   = SPM.XYZ(:,j);
	SPM.XYZmm = SPM.XYZmm(:,j);
	if isfield(SPM,'Rd'), SPM.Rd = SPM.Rd(:,j); end
end

%-Call 'list' functionality to produce table
%-----------------------------------------------------------------------
varargout = {spm_list('list',SPM,hReg,Num,Dis)};





%=======================================================================
case 'txtlist'                                  %-Print ASCII text table
%=======================================================================
% FORMAT spm_list('TxtList',TabDat,c)

if nargin<2, error('Insufficient arguments'), end
if nargin<3, c=1; else, c=varargin{3}; end
TabDat = varargin{2};

%-Table Title
%-----------------------------------------------------------------------
fprintf('\n\nSTATISTICS: %s\n',TabDat.tit)
fprintf('%c','='*ones(1,80)), fprintf('\n')

%-Table header
%-----------------------------------------------------------------------
fprintf('%s\t',TabDat.hdr{1,c:end-1}), fprintf('%s\n',TabDat.hdr{1,end})
fprintf('%s\t',TabDat.hdr{2,c:end-1}), fprintf('%s\n',TabDat.hdr{2,end})
fprintf('%c','-'*ones(1,80)), fprintf('\n')

%-Table data
%-----------------------------------------------------------------------
for i = 1:size(TabDat.dat,1)
	for j=c:size(TabDat.dat,2)
		fprintf(TabDat.fmt{j},TabDat.dat{i,j})
		fprintf('\t')
	end
	fprintf('\n')
end
for i=1:max(1,11-size(TabDat.dat,1)), fprintf('\n'), end
fprintf('%s\n',TabDat.str)
fprintf('%c','-'*ones(1,80)), fprintf('\n')

%-Table footer
%-----------------------------------------------------------------------
fprintf('%s\n',TabDat.ftr{:})
fprintf('%c','='*ones(1,80)), fprintf('\n\n')



%=======================================================================
case 'setcoords'                                    %-Co-ordinate change
%=======================================================================
% FORMAT spm_list('SetCoords',xyz,hAx,hReg)
if nargin<3, error('Insufficient arguments'), end
hAx      = varargin{3};
xyz      = varargin{2};
UD       = get(hAx,'UserData');
HlistXYZ = UD.HlistXYZ(ishandle(UD.HlistXYZ));

%-Set all co-ord strings to black
%-----------------------------------------------------------------------
set(HlistXYZ,'Color','k')

%-If co-ord matches a string, highlight it in red
%-----------------------------------------------------------------------
XYZ      = get(HlistXYZ,'UserData');
if iscell(XYZ), XYZ = cat(2,XYZ{:}); end
[null,i,d] = spm_XYZreg('NearestXYZ',xyz,XYZ);
if d<eps
	set(HlistXYZ(i),'Color','r')
end


%=======================================================================
otherwise                                        %-Unknown action string
%=======================================================================
error('Unknown action string')


%=======================================================================
end
%=======================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% spm_list8 (copied from spm 8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = spm_list8(varargin)
% Display and analysis of SPM{.}
% FORMAT TabDat = spm_list('List',SPM,hReg,[Num,Dis,Str])
% Summary list of local maxima for entire volume of interest
% FORMAT TabDat = spm_list('ListCluster',SPM,hReg,[Num,Dis,Str])
% List of local maxima for a single suprathreshold cluster
%
% SPM    - structure containing SPM, distribution & filtering details
%        - required fields are:
% .swd   - SPM working directory - directory containing current SPM.mat
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests
% .STAT  - distribution {Z, T, X or F}
% .df    - degrees of freedom [df{interest}, df{residual}]
% .u     - height threshold
% .k     - extent threshold {voxels}
% .XYZ   - location of voxels {voxel coords}
% .XYZmm - location of voxels {mm}
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}
% .M     - voxels - > mm matrix
% .VOX   - voxel dimensions {mm}
% .Vspm  - mapped statistic image(s)
% .Ps    - uncorrected P values in searched volume (for voxel FDR)
% .Pp    - uncorrected P values of peaks (for peak FDR)
% .Pc    - uncorrected P values of cluster extents (for cluster FDR)
% .uc    - 0.05 critical thresholds for FWEp, FDRp, FWEc, FDRc
% .thresDesc - description of height threshold (string)
%
% (see spm_getSPM for further details of xSPM structures)
%
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Num    - number of maxima per cluster
% Dis    - distance among clusters (mm)
% Str    - header string
%
% TabDat - Structure containing table data
%        - fields are
% .tit   - table Title (string)
% .hdr   - table header (2x12 cell array)
% .fmt   - fprintf format strings for table data (1x12 cell array)
% .str   - table filtering note (string)
% .ftr   - table footnote information (5x2 cell array)
% .dat   - table data (Nx12 cell array)
%
%                           ----------------
%
% FORMAT spm_list('TxtList',TabDat,c)
% Prints a tab-delimited text version of the table
% TabDat - Structure containing table data (format as above)
% c      - Column of table data to start text table at
%          (E.g. c=3 doesn't print set-level results contained in columns 1 & 2)
%                           ----------------
%
% FORMAT spm_list('SetCoords',xyz,hAx,hC)
% Highlighting of table co-ordinates (used by results section registry)
% xyz    - 3-vector of new co-ordinate
% hAx    - table axis (the registry object for tables)
% hReg   - Handle of caller (not used)
%__________________________________________________________________________
%
% spm_list characterizes SPMs (thresholded at u and k) in terms of
% excursion sets (a collection of face, edge and vertex connected
% subsets or clusters).  The corrected significance of the results are
% based on set, cluster and voxel-level inferences using distributional
% approximations from the Theory of Gaussian Fields.  These
% distributions assume that the SPM is a reasonable lattice
% approximation of a continuous random field with known component field
% smoothness.
%
% The p values are based on the probability of obtaining c, or more,
% clusters of k, or more, resels above u, in the volume S analysed =
% P(u,k,c).  For specified thresholds u, k, the set-level inference is
% based on the observed number of clusters C, = P(u,k,C).  For each
% cluster of size K the cluster-level inference is based on P(u,K,1)
% and for each voxel (or selected maxima) of height U, in that cluster,
% the voxel-level inference is based on P(U,0,1).  All three levels of
% inference are supported with a tabular presentation of the p values
% and the underlying statistic:
%
% Set-level     - c    = number of suprathreshold clusters
%               - P    = prob(c or more clusters in the search volume)
%
% Cluster-level - k    = number of voxels in this cluster
%               - Pc   = prob(k or more voxels in the search volume)
%               - Pu   = prob(k or more voxels in a cluster)
%               - Qc   = lowest FDR bound for which this cluster would be
%                        declared positive
%
% Peak-level    - T/F  = Statistic upon which the SPM is based
%               - Ze   = The equivalent Z score - prob(Z > Ze) = prob(t > T)
%               - Pc   = prob(Ze or higher in the search volume)
%               - Qp   = lowest FDR bound for which this peak would be
%                        declared positive
%               - Pu   = prob(Ze or higher at that voxel)
%
% Voxel-level   - Qu   = Expd(Prop of false positives among voxels >= Ze)
%
% x,y,z (mm)    - Coordinates of the voxel
%
% The table is grouped by regions and sorted on the Ze-variate of the
% primary maxima.  Ze-variates (based on the uncorrected p value) are the
% Z score equivalent of the statistic. Volumes are expressed in voxels.
%
% Clicking on values in the table returns the value to the Matlab
% workspace. In addition, clicking on the co-ordinates jumps the
% results section cursor to that location. The table has a context menu
% (obtained by right-clicking in the background of the table),
% providing options to print the current table as a text table, or to
% extract the table data to the Matlab workspace.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston & Andrew Holmes
% $Id: spm_list.m 2821 2009-03-03 19:54:19Z guillaume $


% satellite figure global variable
%--------------------------------------------------------------------------
global SatWindow

% Choose between voxel-wise and topological FDR
%--------------------------------------------------------------------------
defaults = spm('GetGlobal','defaults');
try
    topoFDR = defaults.stats.topoFDR;
catch
    topoFDR = true;
end

%==========================================================================
switch lower(varargin{1}), case 'list'                            %-List
%==========================================================================
% FORMAT TabDat = spm_list('list',SPM,hReg)

    %-Tolerance for p-value underflow, when computing equivalent Z's
    %----------------------------------------------------------------------
    tol = eps*10;

    %-Parse arguments and set maxima number and separation
    %----------------------------------------------------------------------
    if nargin < 2,  error('insufficient arguments'),     end
    if nargin < 3,  hReg = []; else  hReg = varargin{3}; end


    %-Get current location (to highlight selected voxel in table)
    %----------------------------------------------------------------------
    %xyzmm     = spm_results_ui('GetCoords');
    xyzmm = spm_XYZreg('GetCoords',hReg);% added by xu cui
    

    %-Extract data from xSPM
    %----------------------------------------------------------------------
    S     = varargin{2}.S;
    VOX   = varargin{2}.VOX;
    DIM   = varargin{2}.DIM;
    n     = varargin{2}.n;
    STAT  = varargin{2}.STAT;
    df    = varargin{2}.df;
    u     = varargin{2}.u;
    M     = varargin{2}.M;
    k     = varargin{2}.k;
    try, QPs = varargin{2}.Ps; end
    try, QPp = varargin{2}.Pp; end
    try, QPc = varargin{2}.Pc; end
    try
        thresDesc = sprintf('{%s}', varargin{2}.thresDesc);
    catch
        thresDesc = '';
    end
    
    if STAT~='P'
        R     = varargin{2}.R;
        FWHM  = varargin{2}.FWHM;
    end
    try
        units = varargin{2}.units;
    catch
        units = {'mm' 'mm' 'mm'};
    end
    units{1}  = [units{1} ' '];
    units{2}  = [units{2} ' '];

    DIM       = DIM > 1;              % dimensions
    VOX       = VOX(DIM);             % scaling

    if STAT~='P'
        FWHM  = FWHM(DIM);            % Full width at max/2
        FWmm  = FWHM.*VOX;            % FWHM {units}
        v2r   = 1/prod(FWHM);         % voxels to resels
        k     = k*v2r;                % extent threshold in resels
        R(find(~DIM) + 1) = [];       % eliminate null resel counts
        try, QPs = sort(QPs(:)); end  % Needed for voxel FDR
        try, QPp = sort(QPp(:)); end  % Needed for peak FDR
        try, QPc = sort(QPc(:)); end  % Needed for cluster FDR
    end

    %-get number and separation for maxima to be reported
    %----------------------------------------------------------------------
    if length(varargin) > 3
        Num    = varargin{4};         % number of maxima per cluster
        Dis    = varargin{5};         % distance among clusters (mm)
    else
        Num    = 3;
        Dis    = 8;
    end
    if length(varargin) > 5
        Title  = varargin{6};
    else
        Title  = 'p-values adjusted for search volume';
    end

    %-Setup graphics panel
    %----------------------------------------------------------------------
    spm('Pointer','Watch')
    if SatWindow
        Fgraph = SatWindow;
        figure(Fgraph);
    else%edit by xu cui
        Fgraph = figure('unit','normalized','position',[0.4,0.1,0.55,0.5],'Color','w',...
        'Name','volume', 'NumberTitle','off','resize','on','MenuBar','none');
        Fgraph = gcf;
        %Fgraph = spm_figure('GetWin','Graphics');
    end
    spm_results_ui('Clear',Fgraph)
    FS    = spm('FontSizes');           %-Scaled font sizes
    PF    = spm_platform('fonts');      %-Font names (for this platform)


    %-Table header & footer
    %======================================================================

    %-Table axes & Title
    %----------------------------------------------------------------------
    if SatWindow, ht = 0.85; bot = 0.14; else ht = 0.8; bot = 0.15; end % changed by xu cui

    if STAT == 'P'
        Title = 'Posterior Probabilities';
    end

    hAx   = axes('Position',[0.025 bot 0.9 ht],...
                'DefaultTextFontSize',FS(8),...
                'DefaultTextInterpreter','Tex',...
                'DefaultTextVerticalAlignment','Baseline',...
                'Units','points',...
                'Visible','off');

    AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)])
    dy    = FS(9);
    y     = floor(AxPos(4)) - dy;

    text(0,y,['Statistics:  \it\fontsize{',num2str(FS(9)),'}',Title],...
              'FontSize',FS(11),'FontWeight','Bold');   y = y - dy/2;
    line([0 1],[y y],'LineWidth',3,'Color','r'),        y = y - 9*dy/8;

    %-Construct table header
    %----------------------------------------------------------------------
    set(gca,'DefaultTextFontName',PF.helvetica,'DefaultTextFontSize',FS(8))

    Hc = [];
    Hp = [];
    h  = text(0.01,y,   'set-level','FontSize',FS(9));      Hc = [Hc,h];
    h  = line([0,0.11],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r'); Hc = [Hc,h];
    h  = text(0.08,y-9*dy/8,    '\itc ');                   Hc = [Hc,h];
    h  = text(0.02,y-9*dy/8,    '\itp ');                   Hc = [Hc,h];
    Hp = [Hp,h];
    text(0.22,y,        'cluster-level','FontSize',FS(9));
    line([0.14,0.44],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');
    h  = text(0.15,y-9*dy/8,    '\itp\rm_{FWE-corr}');     Hp = [Hp,h];
    h  = text(0.24,y-9*dy/8,    '\itq\rm_{FDR-corr}');     Hp = [Hp,h];
    h  = text(0.39,y-9*dy/8,    '\itp\rm_{uncorr}');       Hp = [Hp,h];
    h  = text(0.34,y-9*dy/8,    '\itk\rm_E');

    text(0.64,y,        'peak-level','FontSize',FS(9));
    line([0.48,0.88],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');
    h  = text(0.49,y-9*dy/8,    '\itp\rm_{FWE-corr}');     Hp = [Hp,h];
    h  = text(0.58,y-9*dy/8,        '\itq\rm_{FDR-corr}'); Hp = [Hp,h];
    h  = text(0.82,y-9*dy/8,    '\itp\rm_{uncorr}');       Hp = [Hp,h];
    h  = text(0.67,y-9*dy/8,     sprintf('\\it%c',STAT));
    h  = text(0.75,y-9*dy/8,    '(\itZ\rm_\equiv)');

    text(0.92,y - dy/2,[units{:}],'Fontsize',FS(8));


    %-Headers for text table...
    %-----------------------------------------------------------------------
    TabDat.tit = Title;
    TabDat.hdr = {  'set',      'c';...
        'set',      'p';...
        'cluster',  'p(FWE-cor)';...
        'cluster',  'p(FDR-cor)';...
        'cluster',  'equivk';...
        'cluster',  'p(unc)';...
        'peak',     'p(FWE-cor)';...
        'peak',     'p(FDR-cor)';...
        'peak',      STAT;...
        'peak',     'equivZ';...
        'peak',     'p(unc)';...
        '',         'x,y,z {mm}'}';...

    TabDat.fmt = {  '%-0.3f','%g',...                          %-Set
        '%0.3f', '%0.3f','%0.0f', '%0.3f',...                  %-Cluster
        '%0.3f', '%0.3f', '%6.2f', '%5.2f', '%0.3f',...        %-Peak
        '%3.0f %3.0f %3.0f'};                                  %-XYZ

    %-Column Locations
    %----------------------------------------------------------------------
    tCol = [ 0.01      0.08 ...                                %-Set
             0.15      0.24      0.33      0.39 ...            %-Cluster
             0.49      0.58      0.65      0.74      0.83 ...  %-Peak
             0.92];                                            %-XYZ

    %-Move to next vertical position marker
    %----------------------------------------------------------------------
    y     = y - 7*dy/4;
    line([0 1],[y y],'LineWidth',1,'Color','r')
    y     = y - 5*dy/4;
    y0    = y;


    %-Table filtering note
    %----------------------------------------------------------------------
    if isinf(Num)
        TabDat.str = sprintf('table shows all local maxima > %.1fmm apart',Dis);
    else
        TabDat.str = sprintf(['table shows %d local maxima ',...
            'more than %.1fmm apart'],Num,Dis);
    end
    text(0.5,4,TabDat.str,'HorizontalAlignment','Center','FontName',PF.helvetica,...
        'FontSize',FS(8),'FontAngle','Italic')


    %-Volume, resels and smoothness (if classical inference)
    %----------------------------------------------------------------------
    line([0 1],[0 0],'LineWidth',1,'Color','r')
    if STAT ~= 'P'
        %------------------------------------------------------------------
        Pz              = spm_P(1,0,u,df,STAT,1,n,S);
        Pu              = spm_P(1,0,u,df,STAT,R,n,S);
        %Qu              = spm_P_FDR(u,df,STAT,n,QPs);
        [P Pn Em En EN] = spm_P(1,k,u,df,STAT,R,n,S);
        
        %-Footnote with SPM parameters
        %------------------------------------------------------------------
        set(gca,'DefaultTextFontName',PF.helvetica,...
            'DefaultTextInterpreter','None','DefaultTextFontSize',FS(8))
        TabDat.ftr    = cell(5,2);
        TabDat.ftr{1} = ...
            sprintf('Height threshold: %c = %0.2f, p = %0.3f (%0.3f)',...
            STAT,u,Pz,Pu);
        TabDat.ftr{2} = ...
            sprintf('Extent threshold: k = %0.0f voxels, p = %0.3f (%0.3f)',...
            k/v2r,Pn,P);
        TabDat.ftr{3} = ...
            sprintf('Expected voxels per cluster, <k> = %0.3f',En/v2r);
        TabDat.ftr{4} = ...
            sprintf('Expected number of clusters, <c> = %0.2f',Em*Pn);
        if any(isnan(varargin{2}.uc))
            TabDat.ftr{5} = ...
            sprintf('FWEp: %0.3f, FDRp: %0.3f',varargin{2}.uc(1:2));
        else
            TabDat.ftr{5} = ...
            sprintf('FWEp: %0.3f, FDRp: %0.3f, FWEc: %0.0f, FDRc: %0.0f',...
            varargin{2}.uc);
        end
        TabDat.ftr{6} = ...
            sprintf('Degrees of freedom = [%0.1f, %0.1f]',df);
        TabDat.ftr{7} = ...
            ['FWHM = ' sprintf('%0.1f ', FWmm) units{:} '; ' ...
            sprintf('%0.1f ', FWHM) '{voxels}'];
        TabDat.ftr{8} = ...
            sprintf('Volume: %0.0f = %0.0f voxels = %0.1f resels', ...
            S*prod(VOX),S,R(end));
        TabDat.ftr{9} = ...
            ['Voxel size: ' sprintf('%0.1f ',VOX) units{:} '; ' ...
            sprintf('(resel = %0.2f voxels)',prod(FWHM))];

        text(0.0,-1*dy,TabDat.ftr{1},...
            'UserData',[u,Pz,Pu],'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.0,-2*dy,TabDat.ftr{2},...
            'UserData',[k/v2r,Pn,P],'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.0,-3*dy,TabDat.ftr{3},...
            'UserData',En/v2r,'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.0,-4*dy,TabDat.ftr{4},...
            'UserData',Em*Pn,'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.0,-5*dy,TabDat.ftr{5},...
            'UserData',varargin{2}.uc,'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.5,-1*dy,TabDat.ftr{6},...
            'UserData',df,'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.5,-2*dy,TabDat.ftr{7},...
            'UserData',FWmm,'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.5,-3*dy,TabDat.ftr{8},...
            'UserData',[S*prod(VOX),S,R(end)],...
            'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.5,-4*dy,TabDat.ftr{9},...
            'UserData',[VOX,prod(FWHM)],...
            'ButtonDownFcn','get(gcbo,''UserData'')')
    else
        TabDat.ftr = {};
    end


    %-Characterize excursion set in terms of maxima
    % (sorted on Z values and grouped by regions)
    %======================================================================
    if isempty(varargin{2}.Z)
        text(0.5,y-6*dy,'no suprathreshold clusters',...
            'HorizontalAlignment','Center',...
            'FontAngle','Italic','FontWeight','Bold',...
            'FontSize',FS(16),'Color',[1,1,1]*.5);
        TabDat.dat = cell(0,12);
        varargout  = {TabDat};
        spm('Pointer','Arrow')
        return
    end

    % Includes Darren Gitelman's code for working around
    % spm_max for conjunctions with negative thresholds
    %----------------------------------------------------------------------
    minz        = abs(min(min(varargin{2}.Z)));
    zscores     = 1 + minz + varargin{2}.Z;
    [N Z XYZ A] = spm_max(zscores,varargin{2}.XYZ);
    Z           = Z - minz - 1;

    %-Convert cluster sizes from voxels to resels
    %----------------------------------------------------------------------
    if STAT~='P'
        if isfield(varargin{2},'VRvp')
            V2R = spm_get_data(varargin{2}.VRvp,XYZ);
        else
            V2R = v2r;
        end
        N       = N.*V2R;
    end

    %-Convert maxima locations from voxels to mm
    %----------------------------------------------------------------------
    XYZmm = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];



    %-Table proper (& note all data in cell array)
    %======================================================================

    %-Pagination variables
    %----------------------------------------------------------------------
    hPage = [];
    set(gca,'DefaultTextFontName',PF.courier,'DefaultTextFontSize',FS(7))


    %-Set-level p values {c} - do not display if reporting a single cluster
    %----------------------------------------------------------------------
    c     = max(A);                                    %-Number of clusters
    if STAT ~= 'P'
        Pc    = spm_P(c,k,u,df,STAT,R,n,S);            %-Set-level p-value
    else
        Pc    = [];
        set(Hp,'Visible','off')
    end

    if c > 1;
        h     = text(tCol(1),y,sprintf(TabDat.fmt{1},Pc),'FontWeight','Bold',...
            'UserData',Pc,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(2),y,sprintf(TabDat.fmt{2},c),'FontWeight','Bold',...
            'UserData',c,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
    else
        set(Hc,'Visible','off')
    end

    TabDat.dat = {Pc,c};            %-Table data
    TabLin     = 1;                 %-Table data line


    %-Local maxima p-values & statistics
    %----------------------------------------------------------------------
    HlistXYZ = [];
    while numel(find(isfinite(Z)))

        % Paginate if necessary
        %------------------------------------------------------------------
        if y < min(Num + 1,3)*dy

            % added Fgraph term to paginate on Satellite window
            %--------------------------------------------------------------
            h     = text(0.5,-5*dy,...
                sprintf('Page %d',spm_figure('#page',Fgraph)),...
                'FontName',PF.helvetica,'FontAngle','Italic',...
                'FontSize',FS(8));

            spm_figure('NewPage',[hPage,h])
            hPage = [];
            y     = y0;
        end

        %-Find largest remaining local maximum
        %------------------------------------------------------------------
        [U,i]   = max(Z);           % largest maxima
        j       = find(A == A(i));  % maxima in cluster


        %-Compute cluster {k} and peak-level {u} p values for this cluster
        %------------------------------------------------------------------
        if STAT ~= 'P'
            Nv      = N(i)/v2r;                       % extent {voxels}
            
            Pz      = spm_P(1,0,   U,df,STAT,1,n,S);  % uncorrected p value
            Pu      = spm_P(1,0,   U,df,STAT,R,n,S);  % FWE-corrected {based on Z}
            [Pk Pn] = spm_P(1,N(i),u,df,STAT,R,n,S);  % [un]corrected {based on k}
            if topoFDR
                Qc  = spm_P_clusterFDR(N(i),df,STAT,R,n,u,QPc); % cluster FDR-corrected {based on k}
                Qp  = spm_P_peakFDR(U,df,STAT,R,n,u,QPp); % peak FDR-corrected {based on Z}
                Qu  = [];
            else
                Qu  = spm_P_FDR(   U,df,STAT,n,QPs);  % voxel FDR-corrected {based on Z}
                Qc  = [];
                Qp  = [];
            end

            if Pz < tol                               % Equivalent Z-variate
                Ze  = Inf;                            % (underflow => can't compute)
            else
                Ze  = spm_invNcdf(1 - Pz);
            end
        else
            Nv      = N(i);
            
            Pz      = [];
            Pu      = [];
            Qu      = [];
            Pk      = [];
            Pn      = [];
            Qc      = [];
            Qp      = [];
            Ze      = spm_invNcdf(U);
        end


        %-Print cluster and maximum peak-level p values {Z}
        %------------------------------------------------------------------
        h     = text(tCol(3),y,sprintf(TabDat.fmt{3},Pk),'FontWeight','Bold',...
            'UserData',Pk,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(4),y,sprintf(TabDat.fmt{4},Qc),'FontWeight','Bold',...
            'UserData',Qc,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(5),y,sprintf(TabDat.fmt{5},Nv),'FontWeight','Bold',...
            'UserData',Nv,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(6),y,sprintf(TabDat.fmt{6},Pn),'FontWeight','Bold',...
            'UserData',Pn,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];

        h     = text(tCol(7),y,sprintf(TabDat.fmt{7},Pu),'FontWeight','Bold',...
            'UserData',Pu,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        if topoFDR
        h     = text(tCol(8),y,sprintf(TabDat.fmt{8},Qp),'FontWeight','Bold',...
            'UserData',Qp,'ButtonDownFcn','get(gcbo,''UserData'')');
        else
        h     = text(tCol(8),y,sprintf(TabDat.fmt{8},Qu),'FontWeight','Bold',...
            'UserData',Qu,'ButtonDownFcn','get(gcbo,''UserData'')');
        end
        hPage = [hPage, h];
        h     = text(tCol(9),y,sprintf(TabDat.fmt{9},U),'FontWeight','Bold',...
            'UserData',U,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(10),y,sprintf(TabDat.fmt{10},Ze),'FontWeight','Bold',...
            'UserData',Ze,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = ...
            text(tCol(11),y,sprintf(TabDat.fmt{11},Pz),'FontWeight','Bold',...
            'UserData',Pz,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];

        % Specifically changed so it properly finds hMIPax
        %------------------------------------------------------------------
        h     = text(tCol(12),y,sprintf(TabDat.fmt{12},XYZmm(:,i)),...
            'FontWeight','Bold',...
            'Tag','ListXYZ',...
            'ButtonDownFcn',[...
            'hMIPax = findobj(''tag'',''hMIPax'');',...
            'spm_mip_ui(''SetCoords'',get(gcbo,''UserData''),hMIPax);'],...
            'Interruptible','off','BusyAction','Cancel',...
            'UserData',XYZmm(:,i));

        HlistXYZ = [HlistXYZ, h];
        if spm_XYZreg('Edist',xyzmm,XYZmm(:,i))<tol && ~isempty(hReg)
            set(h,'Color','r')
        end
        hPage  = [hPage, h];

        y      = y - dy;

        if topoFDR
        [TabDat.dat{TabLin,3:12}] = deal(Pk,Qc,Nv,Pn,Pu,Qp,U,Ze,Pz,XYZmm(:,i));
        else
        [TabDat.dat{TabLin,3:12}] = deal(Pk,Qc,Nv,Pn,Pu,Qu,U,Ze,Pz,XYZmm(:,i));
        end
        TabLin = TabLin + 1;

        %-Print Num secondary maxima (> Dis mm apart)
        %------------------------------------------------------------------
        [l q] = sort(-Z(j));                % sort on Z value
        D     = i;
        for i = 1:length(q)
            d = j(q(i));
            if min(sqrt(sum((XYZmm(:,D)-XYZmm(:,d)*ones(1,size(D,2))).^2)))>Dis;

                if length(D) < Num

                    % Paginate if necessary
                    %------------------------------------------------------
                    if y < dy
                        h = text(0.5,-5*dy,sprintf('Page %d',...
                            spm_figure('#page',Fgraph)),...
                            'FontName',PF.helvetica,...
                            'FontAngle','Italic',...
                            'FontSize',FS(8));

                        spm_figure('NewPage',[hPage,h])
                        hPage = [];
                        y     = y0;
                    end

                    % voxel-level p values {Z}
                    %------------------------------------------------------
                    if STAT ~= 'P'
                        Pz    = spm_P(1,0,Z(d),df,STAT,1,n,S);
                        Pu    = spm_P(1,0,Z(d),df,STAT,R,n,S);
                        if topoFDR
                            Qp = spm_P_peakFDR(Z(d),df,STAT,R,n,u,QPp);
                            Qu = [];
                        else
                            Qu = spm_P_FDR(Z(d),df,STAT,n,QPs);
                            Qp = [];
                        end
                        if Pz < tol
                            Ze = Inf;
                        else
                            Ze = spm_invNcdf(1 - Pz); 
                        end
                    else
                        Pz    = [];
                        Pu    = [];
                        Qu    = [];
                        Qp    = [];
                        Ze    = spm_invNcdf(Z(d));
                    end

                    h     = text(tCol(7),y,sprintf(TabDat.fmt{7},Pu),...
                        'UserData',Pu,...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    hPage = [hPage, h];

                    if topoFDR
                    h     = text(tCol(8),y,sprintf(TabDat.fmt{8},Qp),...
                        'UserData',Qp,...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    else
                    h     = text(tCol(8),y,sprintf(TabDat.fmt{8},Qu),...
                        'UserData',Qu,...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    end
                    hPage = [hPage, h];
                    h     = text(tCol(9),y,sprintf(TabDat.fmt{9},Z(d)),...
                        'UserData',Z(d),...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    hPage = [hPage, h];
                    h     = text(tCol(10),y,sprintf(TabDat.fmt{10},Ze),...
                        'UserData',Ze,...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    hPage = [hPage, h];
                    h     = text(tCol(11),y,sprintf(TabDat.fmt{11},Pz),...
                        'UserData',Pz,...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    hPage = [hPage, h];

                    % specifically modified line to use hMIPax
                    %------------------------------------------------------
                    h     = text(tCol(12),y,...
                        sprintf(TabDat.fmt{12},XYZmm(:,d)),...
                        'Tag','ListXYZ',...
                        'ButtonDownFcn',[...
                        'hMIPax = findobj(''tag'',''hMIPax'');',...
                        'spm_mip_ui(''SetCoords'',',...
                        'get(gcbo,''UserData''),hMIPax);'],...
                        'Interruptible','off','BusyAction','Cancel',...
                        'UserData',XYZmm(:,d));

                    HlistXYZ = [HlistXYZ, h];
                    if spm_XYZreg('Edist',xyzmm,XYZmm(:,d))<tol && ...
                            ~isempty(hReg)
                        set(h,'Color','r')
                    end
                    hPage = [hPage, h];
                    D     = [D d];
                    y     = y - dy;
                    if topoFDR
                    [TabDat.dat{TabLin,7:12}] = ...
                        deal(Pu,Qp,Z(d),Ze,Pz,XYZmm(:,d));
                    else
                    [TabDat.dat{TabLin,7:12}] = ...
                        deal(Pu,Qu,Z(d),Ze,Pz,XYZmm(:,d));
                    end
                    TabLin = TabLin+1;
                end
            end
        end
        Z(j) = NaN;     % Set local maxima to NaN
    end             % end region


    %-Number and register last page (if paginated)
    %-Changed to use Fgraph for numbering
    %----------------------------------------------------------------------
    if spm_figure('#page',Fgraph)>1
        h = text(0.5,-5*dy,sprintf('Page %d/%d',spm_figure('#page',Fgraph)*[1,1]),...
            'FontName',PF.helvetica,'FontSize',FS(8),'FontAngle','Italic');
        spm_figure('NewPage',[hPage,h])
    end

    %-End: Store TabDat in UserData of axes & reset pointer
    %======================================================================
    h      = uicontextmenu('Tag','TabDat',...
        'UserData',TabDat);
    set(gca,'UIContextMenu',h,...
        'Visible','on',...
        'XColor','w','YColor','w')
    uimenu(h,'Label','Table')
    uimenu(h,'Separator','on','Label','Print text table',...
        'Tag','TD_TxtTab',...
        'CallBack',...
        'spm_list(''txtlist'',get(get(gcbo,''Parent''),''UserData''),3)',...
        'Interruptible','off','BusyAction','Cancel');
    uimenu(h,'Separator','off','Label','Extract table data structure',...
        'Tag','TD_Xdat',...
        'CallBack','get(get(gcbo,''Parent''),''UserData'')',...
        'Interruptible','off','BusyAction','Cancel');
    uimenu(h,'Separator','on','Label','help',...
        'Tag','TD_Xdat',...
        'CallBack','spm_help(''spm_list'')',...
        'Interruptible','off','BusyAction','Cancel');

    %-Setup registry
    %----------------------------------------------------------------------
    set(hAx,'UserData',struct('hReg',hReg,'HlistXYZ',HlistXYZ))
    spm_XYZreg('Add2Reg',hReg,hAx,'spm_list');

    %-Return TabDat structure & reset pointer
    %----------------------------------------------------------------------
    varargout = {TabDat};
    spm('Pointer','Arrow')


    %======================================================================
    case 'listcluster'                      %-List for current cluster only
    %======================================================================
    % FORMAT TabDat = spm_list('listcluster',SPM,hReg)

        spm('Pointer','Watch')

        %-Parse arguments
        %------------------------------------------------------------------
        if nargin < 2,  error('insufficient arguments'),     end
        if nargin < 3,  hReg = []; else hReg = varargin{3}; end
        SPM    = varargin{2};

        %-get number and separation for maxima to be reported
        %------------------------------------------------------------------
        if length(varargin) > 3

            Num    = varargin{4};       % number of maxima per cluster
            Dis    = varargin{5};       % distance among clusters (mm)
        else
            Num    = 32;
            Dis    = 4;
        end


        %-if there are suprathreshold voxels, filter out all but current cluster
        %------------------------------------------------------------------
        if ~isempty(SPM.Z)

            %-Jump to voxel nearest current location
            %--------------------------------------------------------------
            [xyzmm,i] = spm_XYZreg('NearestXYZ',...
                spm_results_ui('GetCoords'),SPM.XYZmm);
            spm_results_ui('SetCoords',SPM.XYZmm(:,i));

            %-Find selected cluster
            %--------------------------------------------------------------
            A         = spm_clusters(SPM.XYZ);
            j         = find(A == A(i));
            SPM.Z     = SPM.Z(j);
            SPM.XYZ   = SPM.XYZ(:,j);
            SPM.XYZmm = SPM.XYZmm(:,j);
            if isfield(SPM,'Rd'), SPM.Rd = SPM.Rd(:,j); end
        end

        %-Call 'list' functionality to produce table
        %------------------------------------------------------------------
        varargout = {spm_list('list',SPM,hReg,Num,Dis)};


    %======================================================================
    case 'txtlist'                                 %-Print ASCII text table
    %======================================================================
    % FORMAT spm_list('TxtList',TabDat,c)

        if nargin<2, error('Insufficient arguments'), end
        if nargin<3, c=1; else c=varargin{3}; end
        TabDat = varargin{2};

        %-Table Title
        %------------------------------------------------------------------
        fprintf('\n\nSTATISTICS: %s\n',TabDat.tit)
        fprintf('%c',repmat('=',1,80)), fprintf('\n')

        %-Table header
        %------------------------------------------------------------------
        fprintf('%s\t',TabDat.hdr{1,c:end-1}), fprintf('%s\n',TabDat.hdr{1,end})
        fprintf('%s\t',TabDat.hdr{2,c:end-1}), fprintf('%s\n',TabDat.hdr{2,end})
        fprintf('%c',repmat('-',1,80)), fprintf('\n')

        %-Table data
        %------------------------------------------------------------------
        for i = 1:size(TabDat.dat,1)
            for j=c:size(TabDat.dat,2)
                fprintf(TabDat.fmt{j},TabDat.dat{i,j})
                fprintf('\t')
            end
            fprintf('\n')
        end
        for i=1:max(1,12-size(TabDat.dat,1)), fprintf('\n'), end
        fprintf('%s\n',TabDat.str)
        fprintf('%c',repmat('-',1,80)), fprintf('\n')

        %-Table footer
        %------------------------------------------------------------------
        fprintf('%s\n',TabDat.ftr{:})
        fprintf('%c',repmat('=',1,80)), fprintf('\n\n')



        %==================================================================
    case 'setcoords'                                   %-Co-ordinate change
        %==================================================================
        % FORMAT spm_list('SetCoords',xyz,hAx,hReg)
        if nargin<3, error('Insufficient arguments'), end
        hAx      = varargin{3};
        xyz      = varargin{2};
        UD       = get(hAx,'UserData');
        HlistXYZ = UD.HlistXYZ(ishandle(UD.HlistXYZ));

        %-Set all co-ord strings to black
        %------------------------------------------------------------------
        set(HlistXYZ,'Color','k')

        %-If co-ord matches a string, highlight it in red
        %------------------------------------------------------------------
        XYZ      = get(HlistXYZ,'UserData');
        if iscell(XYZ), XYZ = cat(2,XYZ{:}); end
        [null,i,d] = spm_XYZreg('NearestXYZ',xyz,XYZ);
        if d<eps
            set(HlistXYZ(i),'Color','r')
        end

        %==================================================================
    otherwise                                       %-Unknown action string
        %==================================================================
        error('Unknown action string')
end
%==========================================================================


function []=setcolormap(what)

switch lower(what), case 'gray'
	colormap(gray(64))
case 'hot'
	colormap(hot(64))
case 'pink'
	colormap(pink(64))
case 'gray-hot'
	tmp = hot(64 + 16);  tmp = tmp([1:64] + 16,:);
	colormap([gray(64); tmp])
case 'gray-cold'
	tmp = jet(64 + 48);  tmp = tmp([1:64] + 16,:);
	colormap([gray(64); tmp])    
case 'gray-hot-cold'
	tmp = jet(64 + 16);  tmp = tmp([1:64] + 16,:);
	colormap([gray(64); tmp])
case 'gray-pink'
	tmp = pink(64 + 16); tmp = tmp([1:64] + 16,:);
	colormap([gray(64); tmp])
case 'invert'
	colormap(flipud(colormap))
case 'brighten'
	colormap(brighten(colormap, 0.2))
case 'darken'
	colormap(brighten(colormap, -0.2))
otherwise
	error('Illegal ColAction specification')
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% spm_orthviews
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = spm_orthviews(action,varargin)
% Display Orthogonal Views of a Normalized Image
% FORMAT H = spm_orthviews('Image',filename[,position])
% filename - name of image to display
% area     - position of image
%            -  area(1) - position x
%            -  area(2) - position y
%            -  area(3) - size x
%            -  area(4) - size y
% H        - handle for ortho sections
% FORMAT spm_orthviews('BB',bb)
% bb       - bounding box
%            [loX loY loZ
%             hiX hiY hiZ]
%
% FORMAT spm_orthviews('Redraw')
% Redraws the images
%
% FORMAT spm_orthviews('Reposition',centre)
% centre   - X, Y & Z coordinates of centre voxel
%
% FORMAT spm_orthviews('Space'[,handle])
% handle   - the view to define the space by
% with no arguments - puts things into mm space
%
% FORMAT spm_orthviews('MaxBB')
% sets the bounding box big enough display the whole of all images
%
% FORMAT spm_orthviews('Resolution',res)
% res      - resolution (mm)
%
% FORMAT spm_orthviews('Delete', handle)
% handle   - image number to delete
%
% FORMAT spm_orthviews('Reset')
% clears the orthogonal views
%
% FORMAT spm_orthviews('Pos')
% returns the co-ordinate of the crosshairs in millimetres in the
% standard space.
%
% FORMAT spm_orthviews('Pos', i)
% returns the voxel co-ordinate of the crosshairs in the image in the
% ith orthogonal section.
%
% FORMAT spm_orthviews('Xhairs','off') OR spm_orthviews('Xhairs')
% disables the cross-hairs on the display.
%
% FORMAT spm_orthviews('Xhairs','on')
% enables the cross-hairs.
%
% FORMAT spm_orthviews('Interp',hld)
% sets the hold value to hld (see spm_slice_vol).
%
% FORMAT spm_orthviews('AddBlobs',handle,XYZ,Z,mat)
% Adds blobs from a pointlist to the image specified by the handle(s).
% handle   - image number to add blobs to
% XYZ      - blob voxel locations (currently in millimeters)
% Z        - blob voxel intensities
% mat      - matrix from millimeters to voxels of blob.
% This method only adds one set of blobs, and displays them using a
% split colour table.
%
% FORMAT spm_orthviews('AddColouredBlobs',handle,XYZ,Z,mat,colour)
% Adds blobs from a pointlist to the image specified by the handle(s).
% handle   - image number to add blobs to
% XYZ      - blob voxel locations (currently in millimeters)
% Z        - blob voxel intensities
% mat      - matrix from millimeters to voxels of blob.
% colour   - the 3 vector containing the colour that the blobs should be
% Several sets of blobs can be added in this way, and it uses full colour.
% Although it may not be particularly attractive on the screen, the colour
% blobs print well.
%
% FORMAT spm_orthviews('AddColourBar',handle,blobno)
% Adds colourbar for a specified blob set
% handle   - image number
% blobno   - blob number
%
% FORMAT spm_orthviews('Register',hReg)
% See spm_XYZreg for more information.
%
% FORMAT spm_orthviews('RemoveBlobs',handle)
% Removes all blobs from the image specified by the handle(s).
%
% spm_orthviews('Register',hReg)
% hReg      - Handle of HandleGraphics object to build registry in.
% See spm_XYZreg for more information.
%
% spm_orthviews('AddContext',handle)
% handle   - image number to add context menu to
%
% spm_orthviews('RemoveContext',handle)
% handle   - image number to remove context menu from
%
% CONTEXT MENU
% spm_orthviews offers many of its features in a context menu, which is
% accessible via the right mouse button in each displayed image.
%
% PLUGINS
% The display capabilities of spm_orthviews can be extended with
% plugins. These are located in the spm_orthviews subdirectory of the SPM
% distribution. Currently there are 3 plugins available:
% quiver    Add Quiver plots to a displayed image
% quiver3d  Add 3D Quiver plots to a displayed image
% roi       ROI creation and modification
% The functionality of plugins can be accessed via calls to
% spm_orthviews('plugin_name', plugin_arguments). For detailed descriptions
% of each plugin see help spm_orthviews/spm_ov_'plugin_name'.
%
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner, Matthew Brett, Tom Nichols and Volkmar Glauche
% $Id: spm_orthviews.m 601 2006-08-22 08:34:24Z volkmar $



% The basic fields of st are:
%         n        - the number of images currently being displayed
%         vols     - a cell array containing the data on each of the
%                    displayed images.
%         Space    - a mapping between the displayed images and the
%                    mm space of each image.
%         bb       - the bounding box of the displayed images.
%         centre   - the current centre of the orthogonal views
%         callback - a callback to be evaluated on a button-click.
%         xhairs   - crosshairs off/on
%         hld      - the interpolation method
%         fig      - the figure that everything is displayed in
%         mode     - the position/orientation of the sagittal view.
%                    - currently always 1
% 
%         st.registry.hReg \_ See spm_XYZreg for documentation
%         st.registry.hMe  /
% 
% For each of the displayed images, there is a non-empty entry in the
% vols cell array.  Handles returned by "spm_orthviews('Image',.....)"
% indicate the position in the cell array of the newly created ortho-view.
% Operations on each ortho-view require the handle to be passed.
% 
% When a new image is displayed, the cell entry contains the information
% returned by spm_vol (type help spm_vol for more info).  In addition,
% there are a few other fields, some of which I will document here:
% 
%         premul - a matrix to premultiply the .mat field by.  Useful
%                  for re-orienting images.
%         window - either 'auto' or an intensity range to display the
%                  image with.
%         mapping- Mapping of image intensities to grey values. Currently
%                  one of 'linear', 'histeq', loghisteq',
%                  'quadhisteq'. Default is 'linear'.
%                  Histogram equalisation depends on the image toolbox
%                  and is only available if there is a license available
%                  for it.
%         ax     - a cell array containing an element for the three
%                  views.  The fields of each element are handles for
%                  the axis, image and crosshairs.
% 
%         blobs - optional.  Is there for using to superimpose blobs.
%                 vol     - 3D array of image data
%                 mat     - a mapping from vox-to-mm (see spm_vol, or
%                           help on image formats).
%                 max     - maximum intensity for scaling to.  If it
%                           does not exist, then images are auto-scaled.
% 
%                 There are two colouring modes: full colour, and split
%                 colour.  When using full colour, there should be a
%                 'colour' field for each cell element.  When using
%                 split colourscale, there is a handle for the colorbar
%                 axis.
% 
%                 colour  - if it exists it contains the
%                           red,green,blue that the blobs should be
%                           displayed in.
%                 cbar    - handle for colorbar (for split colourscale).
%
% PLUGINS
% The plugin concept has been developed to extend the display capabilities
% of spm_orthviews without the need to rewrite parts of it. Interaction
% between spm_orthviews and plugins takes place
% a) at startup: The subfunction 'reset_st' looks for files with a name
%                spm_ov_PLUGINNAME.m in the directory 'SWD/spm_orthviews'.
%                For each such file, PLUGINNAME will be added to the list
%                st.plugins{:}.
%                The subfunction 'add_context' calls each plugin with
%                feval(['spm_ov_', st.plugins{k}], ...
%			  'context_menu', i, parent_menu)
%                Each plugin may add its own submenu to the context
%                menu.
% b) at redraw:  After images and blobs of st.vols{i} are drawn, the
%                struct st.vols{i} is checked for field names that occur in
%                the plugin list st.plugins{:}. For each matching entry, the
%                corresponding plugin is called with the command 'redraw':
%                feval(['spm_ov_', st.plugins{k}], ...
%			  'redraw', i, TM0, TD, CM0, CD, SM0, SD);
%                The values of TM0, TD, CM0, CD, SM0, SD are defined in the
%                same way as in the redraw subfunction of spm_orthviews.
%                It is up to the plugin to do all necessary redraw
%                operations for its display contents. Each displayed item
%                must have set its property 'HitTest' to 'off' to let events
%                go through to the underlying axis, which is responsible for
%                callback handling. The order in which plugins are called is
%                undefined.

global st;

if isempty(st), reset_st; end;

spm('Pointer','watch');

if nargin == 0, action = ''; end;
action = lower(action);

switch lower(action),
case 'image',
	H = specify_image(varargin{1});
	if ~isempty(H)
		st.vols{H}.area = [0 0 1 1];
		if length(varargin)>=2, st.vols{H}.area = varargin{2}; end;
		if isempty(st.bb), st.bb = maxbb; end;
		bbox;
		cm_pos;
	end;
	varargout{1} = H;
	st.centre    = mean(maxbb);
	redraw_all

case 'bb',
	if length(varargin)> 0 & all(size(varargin{1})==[2 3]), st.bb = varargin{1}; end;
	bbox;
	redraw_all;

case 'redraw',
	redraw_all;
	eval(st.callback);
	if isfield(st,'registry'),
		spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
	end;

case 'reposition',
	if length(varargin)<1, tmp = findcent;
	else, tmp = varargin{1}; end;
	if length(tmp)==3
                h = valid_handles(st.snap);
                if ~isempty(h)
                        tmp=st.vols{h(1)}.mat*...
                            round(inv(st.vols{h(1)}.mat)*[tmp; ...
                                            1]);
                end;
                st.centre = tmp(1:3); 
        end;
	redraw_all;
	eval(st.callback);
	if isfield(st,'registry'),
		spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
	end;
	cm_pos;

case 'setcoords',
	st.centre = varargin{1};
	st.centre = st.centre(:);
	redraw_all;
	eval(st.callback);
	cm_pos;

case 'space',
	if length(varargin)<1,
		st.Space = eye(4);
		st.bb = maxbb;
		bbox;
		redraw_all;
	else,
		space(varargin{1});
		bbox;
		redraw_all;
	end;

case 'maxbb',
	st.bb = maxbb;
	bbox;
	redraw_all;

case 'resolution',
	resolution(varargin{1});
	bbox;
	redraw_all;

case 'window',
	if length(varargin)<2,
		win = 'auto';
	elseif length(varargin{2})==2,
		win = varargin{2};
	end;
	for i=valid_handles(varargin{1}),
		st.vols{i}.window = win;
	end;
	redraw(varargin{1});

case 'delete',
	my_delete(varargin{1});

case 'move',
	move(varargin{1},varargin{2});
	% redraw_all;

case 'reset',
	my_reset;

case 'pos',
	if isempty(varargin),
		H = st.centre(:);
	else,
		H = pos(varargin{1});
	end;
	varargout{1} = H;

case 'interp',
	st.hld = varargin{1};
	redraw_all;

case 'xhairs',
	xhairs(varargin{1});

case 'register',
	register(varargin{1});

case 'addblobs',
	addblobs(varargin{1}, varargin{2},varargin{3},varargin{4});
	% redraw(varargin{1});

case 'addcolouredblobs',
	addcolouredblobs(varargin{1}, varargin{2},varargin{3},varargin{4},varargin{5});
	% redraw(varargin{1});

case 'addimage',
	addimage(varargin{1}, varargin{2});
	% redraw(varargin{1});

case 'addcolouredimage',
	addcolouredimage(varargin{1}, varargin{2},varargin{3});
	% redraw(varargin{1});

case 'addtruecolourimage',
	% spm_orthviews('Addtruecolourimage',handle,filename,colourmap,prop,mx,mn)
	% Adds blobs from an image in true colour
	% handle   - image number to add blobs to [default 1]
	% filename of image containing blob data [default - request via GUI]
	% colourmap - colormap to display blobs in [GUI input]
	% prop - intensity proportion of activation cf grayscale [0.4]
	% mx   - maximum intensity to scale to [maximum value in activation image]
	% mn   - minimum intensity to scale to [minimum value in activation image]
	%
	if nargin < 2
		varargin(1) = {1};
	end
	if nargin < 3
		varargin(2) = {spm_select(1, 'image', 'Image with activation signal')};
	end
	if nargin < 4
		actc = [];
		while isempty(actc)
			actc = getcmap(spm_input('Colourmap for activation image', '+1','s'));
		end
		varargin(3) = {actc};
	end
	if nargin < 5
		varargin(4) = {0.4};
	end
	if nargin < 6
		actv = spm_vol(varargin{2});
		varargin(5) = {max([eps maxval(actv)])};
	end
	if nargin < 7
		varargin(6) = {min([0 minval(actv)])};
	end

	addtruecolourimage(varargin{1}, varargin{2},varargin{3}, varargin{4}, ...
	                   varargin{5}, varargin{6});
	% redraw(varargin{1});

case 'addcolourbar',
    addcolourbar(varargin{1}, varargin{2});
        
case 'rmblobs',
	rmblobs(varargin{1});
	% redraw(varargin{1});

case 'addcontext',
	if nargin == 1,
		handles = 1:24;
	else,
		handles = varargin{1};
	end;
	addcontexts(handles);

case 'rmcontext',
	if nargin == 1,
		handles = 1:24;
	else,
		handles = varargin{1};
	end;
	rmcontexts(handles);

case 'context_menu',
	c_menu(varargin{:});

case 'valid_handles',
	if nargin == 1
		handles = 1:24;
	else,
		handles = varargin{1};
	end;
	varargout{1} = valid_handles(handles);

otherwise,
  addonaction = strcmp(st.plugins,action);
  if any(addonaction)
    feval(['spm_ov_' st.plugins{addonaction}],varargin{:});
  else
    warning('Unknown action string')
  end;
end;

spm('Pointer');
return;


%_______________________________________________________________________
%_______________________________________________________________________
function addblobs(handle, xyz, t, mat)
global st
global TMAX_
global TMIN_
for i=valid_handles(handle),
	if ~isempty(xyz),
		rcp      = round(xyz);
		dim      = max(rcp,[],2)';
		off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
		vol      = zeros(dim)+NaN;
		vol(off) = t;
		vol      = reshape(vol,dim);
		st.vols{i}.blobs=cell(1,1);
		if st.mode == 0,
			axpos = get(st.vols{i}.ax{2}.ax,'Position');
		else,
			axpos = get(st.vols{i}.ax{1}.ax,'Position');
		end;
		mx = max([eps max(t)]);
		mn = min([0 min(t)]);        
        if ~strcmp(TMAX_, 'auto')
            mx = str2num(TMAX_);
            %mn = -str2num(TMAX_);
        end
        if ~strcmp(TMIN_, 'auto')
            mn = str2num(TMIN_);
        end
		st.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'max',mx, 'min',mn);
		addcolourbar(handle,1);
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addimage(handle, fname)
global st
for i=valid_handles(handle),
	if isstruct(fname),
		vol = fname(1);
	else,
		vol = spm_vol(fname);
	end;
	mat = vol.mat;
	st.vols{i}.blobs=cell(1,1);
	mx = max([eps maxval(vol)]);
	mn = min([0 minval(vol)]);
	st.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'max',mx,'min',mn);
	addcolourbar(handle,1);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addcolouredblobs(handle, xyz, t, mat,colour)
global st
for i=valid_handles(handle),
	if ~isempty(xyz),
		rcp      = round(xyz);
		dim      = max(rcp,[],2)';
		off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
		vol      = zeros(dim)+NaN;
		vol(off) = t;
		vol      = reshape(vol,dim);
		if ~isfield(st.vols{i},'blobs'),
			st.vols{i}.blobs=cell(1,1);
			bset = 1;
		else,
			bset = length(st.vols{i}.blobs)+1;
		end;
                mx = max([eps maxval(vol)]);
                mn = min([0 minval(vol)]);
		st.vols{i}.blobs{bset} = struct('vol',vol,'mat',mat,'max',mx,'min',mn,'colour',colour);
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addcolouredimage(handle, fname,colour)
global st
for i=valid_handles(handle),
	if isstruct(fname),
		vol = fname(1);
	else,
		vol = spm_vol(fname);
	end;
	mat = vol.mat;
	if ~isfield(st.vols{i},'blobs'),
		st.vols{i}.blobs=cell(1,1);
		bset = 1;
	else,
		bset = length(st.vols{i}.blobs)+1;
	end;
	mx = max([eps maxval(vol)]);
	mn = min([0 minval(vol)]);
	st.vols{i}.blobs{bset} = struct('vol',vol,'mat',mat,'max',mx,'min',mn,'colour',colour);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addtruecolourimage(handle,fname,colourmap,prop,mx,mn)
% adds true colour image to current displayed image  
global st
for i=valid_handles(handle),
	if isstruct(fname),
		vol = fname(1);
	else,
		vol = spm_vol(fname);
	end;
	mat = vol.mat;
	if ~isfield(st.vols{i},'blobs'),
		st.vols{i}.blobs=cell(1,1);
		bset = 1;
	else,
		bset = length(st.vols{i}.blobs)+1;
	end;
	c = struct('cmap', colourmap,'prop',prop);
	st.vols{i}.blobs{bset} = struct('vol',vol,'mat',mat,'max',mx, ...
                                        'min',mn,'colour',c);
	addcolourbar(handle,bset);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addcolourbar(vh,bh)
global st
if st.mode == 0,
    axpos = get(st.vols{vh}.ax{2}.ax,'Position');
else,
    axpos = get(st.vols{vh}.ax{1}.ax,'Position');
end;
st.vols{vh}.blobs{bh}.cbar = axes('Parent',st.fig,...
          'Position',[(axpos(1)+axpos(3)+0.05+(bh-1)*.1) (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
          'Box','on', 'YDir','normal', 'XTickLabel',[], 'XTick',[]);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function rmblobs(handle)
global st
for i=valid_handles(handle),
	if isfield(st.vols{i},'blobs'),
		for j=1:length(st.vols{i}.blobs),
			if isfield(st.vols{i}.blobs{j},'cbar') & ishandle(st.vols{i}.blobs{j}.cbar),
				delete(st.vols{i}.blobs{j}.cbar);
			end;
		end;
		st.vols{i} = rmfield(st.vols{i},'blobs');
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function register(hreg)
global st
tmp = uicontrol('Position',[0 0 1 1],'Visible','off','Parent',st.fig);
h   = valid_handles(1:24);
if ~isempty(h),
	tmp = st.vols{h(1)}.ax{1}.ax;
	st.registry = struct('hReg',hreg,'hMe', tmp);
	spm_XYZreg('Add2Reg',st.registry.hReg,st.registry.hMe, 'spm_orthviews');
else,
	warning('Nothing to register with');
end;
st.centre = spm_XYZreg('GetCoords',st.registry.hReg);
st.centre = st.centre(:);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function xhairs(arg1),
global st
st.xhairs = 0;
opt = 'on';
if ~strcmp(arg1,'on'),
	opt = 'off';
else,
	st.xhairs = 1;
end;
for i=valid_handles(1:24),
	for j=1:3,
		set(st.vols{i}.ax{j}.lx,'Visible',opt);
		set(st.vols{i}.ax{j}.ly,'Visible',opt);  
	end; 
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function H = pos(arg1)
global st
H = [];
for arg1=valid_handles(arg1),
	is = inv(st.vols{arg1}.premul*st.vols{arg1}.mat);
	H = is(1:3,1:3)*st.centre(:) + is(1:3,4);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function my_reset
global st
if ~isempty(st) & isfield(st,'registry') & ishandle(st.registry.hMe),
	delete(st.registry.hMe); st = rmfield(st,'registry');
end;
my_delete(1:24);
reset_st;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function my_delete(arg1)
global st
for i=valid_handles(arg1),
	kids = get(st.fig,'Children');
	for j=1:3,
		if any(kids == st.vols{i}.ax{j}.ax),
			set(get(st.vols{i}.ax{j}.ax,'Children'),'DeleteFcn','');
			delete(st.vols{i}.ax{j}.ax);
		end;
	end;
	st.vols{i} = [];
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function resolution(arg1)
global st
res      = arg1/mean(svd(st.Space(1:3,1:3)));
Mat      = diag([res res res 1]);
st.Space = st.Space*Mat;
st.bb    = st.bb/res;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function move(handle,pos)
global st
for handle = valid_handles(handle),
	st.vols{handle}.area = pos;
end;
bbox;
% redraw(valid_handles(handle));
return;
%_______________________________________________________________________
%_______________________________________________________________________
function bb = maxbb
global st
mn = [Inf Inf Inf];
mx = -mn;
for i=valid_handles(1:24),
	bb = [[1 1 1];st.vols{i}.dim(1:3)];
	c = [	bb(1,1) bb(1,2) bb(1,3) 1
		bb(1,1) bb(1,2) bb(2,3) 1
		bb(1,1) bb(2,2) bb(1,3) 1
		bb(1,1) bb(2,2) bb(2,3) 1
		bb(2,1) bb(1,2) bb(1,3) 1
		bb(2,1) bb(1,2) bb(2,3) 1
		bb(2,1) bb(2,2) bb(1,3) 1
		bb(2,1) bb(2,2) bb(2,3) 1]';
	tc = st.Space\(st.vols{i}.premul*st.vols{i}.mat)*c;
	tc = tc(1:3,:)';
	mx = max([tc ; mx]);
	mn = min([tc ; mn]);
end;
bb = [mn ; mx];
return;
%_______________________________________________________________________
%_______________________________________________________________________
function space(arg1)
global st
if ~isempty(st.vols{arg1})
	num = arg1;
	Mat = st.vols{num}.premul(1:3,1:3)*st.vols{num}.mat(1:3,1:3);
	vox = sqrt(sum(Mat.^2));
	if det(Mat(1:3,1:3))<0, vox(1) = -vox(1); end;
	Mat = diag([vox 1]);
	Space = (st.vols{num}.mat)/Mat;
	bb = [1 1 1;st.vols{num}.dim(1:3)];
	bb = [bb [1;1]];
	bb=bb*Mat';
	bb=bb(:,1:3);
	bb=sort(bb);
	st.Space  = Space;
	st.bb = bb;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function H = specify_image(arg1, arg2)
global st
H=[];
ok = 1;
if isstruct(arg1),
	V = arg1(1);
else,
	try,
		V = spm_vol(arg1);
	catch,
		fprintf('Can not use image "%s"\n', arg1);
		return;
	end;
end;

ii = 1;
while ~isempty(st.vols{ii}), ii = ii + 1; end;

DeleteFcn = ['spm_orthviews(''Delete'',' num2str(ii) ');'];
V.ax = cell(3,1);
for i=1:3,
	ax = axes('Visible','off','DrawMode','fast','Parent',st.fig,'DeleteFcn',DeleteFcn,...
		'YDir','normal','ButtonDownFcn',...
		['if strcmp(get(gcf,''SelectionType''),''normal''),spm_orthviews(''Reposition'');',...
		'elseif strcmp(get(gcf,''SelectionType''),''extend''),spm_orthviews(''Reposition'');',...
		'spm_orthviews(''context_menu'',''ts'',1);end;']);
	d  = image(0,'Tag','Transverse','Parent',ax,...
		'DeleteFcn',DeleteFcn);
	set(ax,'Ydir','normal','ButtonDownFcn',...
		['if strcmp(get(gcf,''SelectionType''),''normal''),spm_orthviews(''Reposition'');',...
		'elseif strcmp(get(gcf,''SelectionType''),''extend''),spm_orthviews(''reposition'');',...
		'spm_orthviews(''context_menu'',''ts'',1);end;']);

	lx = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
	ly = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
	if ~st.xhairs,
		set(lx,'Visible','off');
		set(ly,'Visible','off');
	end;
	V.ax{i} = struct('ax',ax,'d',d,'lx',lx,'ly',ly);
end;
V.premul    = eye(4);
V.window    = 'auto';
V.mapping   = 'linear';
st.vols{ii} = V;

H = ii;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addcontexts(handles)
global st
for ii = valid_handles(handles),
	cm_handle = addcontext(ii);
	for i=1:3,
		set(st.vols{ii}.ax{i}.ax,'UIcontextmenu',cm_handle);
		st.vols{ii}.ax{i}.cm = cm_handle;
	end;
end;
spm_orthviews('reposition',spm_orthviews('pos'));
return;
%_______________________________________________________________________
%_______________________________________________________________________
function rmcontexts(handles)
global st
for ii = valid_handles(handles),
	for i=1:3,
		set(st.vols{ii}.ax{i}.ax,'UIcontextmenu',[]);
		st.vols{ii}.ax{i} = rmfield(st.vols{ii}.ax{i},'cm');
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function bbox
global st
Dims = diff(st.bb)'+1;

TD = Dims([1 2])';
CD = Dims([1 3])';
if st.mode == 0, SD = Dims([3 2])'; else, SD = Dims([2 3])'; end;

un    = get(st.fig,'Units');set(st.fig,'Units','Pixels');
sz    = get(st.fig,'Position');set(st.fig,'Units',un);
sz    = sz(3:4);
sz(2) = sz(2)-40;

for i=valid_handles(1:24),
	area = st.vols{i}.area(:);
	area = [area(1)*sz(1) area(2)*sz(2) area(3)*sz(1) area(4)*sz(2)];
	if st.mode == 0,
		sx   = area(3)/(Dims(1)+Dims(3))/1.02;
	else,
		sx   = area(3)/(Dims(1)+Dims(2))/1.02;
	end;
	sy   = area(4)/(Dims(2)+Dims(3))/1.02;
	s    = min([sx sy]);

	offy = (area(4)-(Dims(2)+Dims(3))*1.02*s)/2 + area(2);
	sky = s*(Dims(2)+Dims(3))*0.02;
	if st.mode == 0,
		offx = (area(3)-(Dims(1)+Dims(3))*1.02*s)/2 + area(1);
		skx = s*(Dims(1)+Dims(3))*0.02;
	else,
		offx = (area(3)-(Dims(1)+Dims(2))*1.02*s)/2 + area(1);
		skx = s*(Dims(1)+Dims(2))*0.02;
	end;

	DeleteFcn = ['spm_orthviews(''Delete'',' num2str(i) ');'];

	% Transverse
	set(st.vols{i}.ax{1}.ax,'Units','pixels', ...
		'Position',[offx offy s*Dims(1) s*Dims(2)],...
		'Units','normalized','Xlim',[0 TD(1)]+0.5,'Ylim',[0 TD(2)]+0.5,...
		'Visible','on','XTick',[],'YTick',[]);

	% Coronal
	set(st.vols{i}.ax{2}.ax,'Units','Pixels',...
		'Position',[offx offy+s*Dims(2)+sky s*Dims(1) s*Dims(3)],...
		'Units','normalized','Xlim',[0 CD(1)]+0.5,'Ylim',[0 CD(2)]+0.5,...
		'Visible','on','XTick',[],'YTick',[]);

	% Sagittal
	if st.mode == 0,
		set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
			'Position',[offx+s*Dims(1)+skx offy s*Dims(3) s*Dims(2)],...
			'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
			'Visible','on','XTick',[],'YTick',[]);
	else,
		set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
			'Position',[offx+s*Dims(1)+skx offy+s*Dims(2)+sky s*Dims(2) s*Dims(3)],...
			'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
			'Visible','on','XTick',[],'YTick',[]);
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function redraw_all
global st
redraw(1:24);
return;
%_______________________________________________________________________
function mx = maxval(vol)
if isstruct(vol),
	mx = -Inf;
	for i=1:vol.dim(3),
		tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
		imx = max(tmp(find(isfinite(tmp))));
		if ~isempty(imx),mx = max(mx,imx);end
	end;
else,
	mx = max(vol(find(isfinite(vol))));
end;
%_______________________________________________________________________
function mn = minval(vol)
if isstruct(vol),
        mn = Inf;
        for i=1:vol.dim(3),
                tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
		imn = min(tmp(find(isfinite(tmp))));
		if ~isempty(imn),mn = min(mn,imn);end
        end;
else,
        mn = min(vol(find(isfinite(vol))));
end;

%_______________________________________________________________________
%_______________________________________________________________________
function redraw(arg1)
global st
bb   = st.bb;
Dims = round(diff(bb)'+1);
is   = inv(st.Space);
cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);
% cent(1)=-40
% cent(2)=-40
% cent(3)=60
for i = valid_handles(arg1),
	M = st.vols{i}.premul*st.vols{i}.mat;
	TM0 = [	1 0 0 -bb(1,1)+1
		0 1 0 -bb(1,2)+1
		0 0 1 -cent(3)
		0 0 0 1];
	TM = inv(TM0*(st.Space\M));
	TD = Dims([1 2]);

	CM0 = [	1 0 0 -bb(1,1)+1
		0 0 1 -bb(1,3)+1
		0 1 0 -cent(2)
		0 0 0 1];
	CM = inv(CM0*(st.Space\M));
	CD = Dims([1 3]);

	if st.mode ==0,
		SM0 = [	0 0 1 -bb(1,3)+1
			0 1 0 -bb(1,2)+1
			1 0 0 -cent(1)
			0 0 0 1];
		SM = inv(SM0*(st.Space\M)); SD = Dims([3 2]);
	else,
		SM0 = [	0  1 0 -bb(1,2)+1
			0  0 1 -bb(1,3)+1
			1  0 0 -cent(1)
			0  0 0 1];
		SM0 = [	0 -1 0 +bb(2,2)+1
			0  0 1 -bb(1,3)+1
			1  0 0 -cent(1)
			0  0 0 1];
		SM = inv(SM0*(st.Space\M));
		SD = Dims([2 3]);
	end;

	ok=1;
	eval('imgt  = (spm_slice_vol(st.vols{i},TM,TD,st.hld))'';','ok=0;');
	eval('imgc  = (spm_slice_vol(st.vols{i},CM,CD,st.hld))'';','ok=0;');
	eval('imgs  = (spm_slice_vol(st.vols{i},SM,SD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.vols{i}.fname);
	else,
                % get min/max threshold
                if strcmp(st.vols{i}.window,'auto')
                        mn = -Inf;
                        mx = Inf;
                else
                        mn = min(st.vols{i}.window);
                        mx = max(st.vols{i}.window);
                end;
                % threshold images
                imgt = max(imgt,mn); imgt = min(imgt,mx);
                imgc = max(imgc,mn); imgc = min(imgc,mx);
                imgs = max(imgs,mn); imgs = min(imgs,mx);
                % compute intensity mapping, if histeq is available
                if license('test','image_toolbox') == 0
                    st.vols{i}.mapping = 'linear';
                end;
                switch st.vols{i}.mapping,
                 case 'linear',
                 case 'histeq',
                  % scale images to a range between 0 and 1
                  imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                  imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                  imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                  img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                  imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                  imgc = reshape(img(numel(imgt1)+[1:numel(imgc1)]),size(imgc1));
                  imgs = reshape(img(numel(imgt1)+numel(imgc1)+[1:numel(imgs1)]),size(imgs1));
                  mn = 0;
                  mx = 1;
                 case 'quadhisteq',
                  % scale images to a range between 0 and 1
                  imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                  imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                  imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                  img  = histeq([imgt1(:).^2; imgc1(:).^2; imgs1(:).^2],1024);
                  imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                  imgc = reshape(img(numel(imgt1)+[1:numel(imgc1)]),size(imgc1));
                  imgs = reshape(img(numel(imgt1)+numel(imgc1)+[1:numel(imgs1)]),size(imgs1));
                  mn = 0;
                  mx = 1;
                 case 'loghisteq',
                  warning off % messy - but it may avoid extra queries
                  imgt = log(imgt-min(imgt(:)));
                  imgc = log(imgc-min(imgc(:)));
                  imgs = log(imgs-min(imgs(:)));
                  warning on
                  imgt(~isfinite(imgt)) = 0;
                  imgc(~isfinite(imgc)) = 0;
                  imgs(~isfinite(imgs)) = 0;
                  % scale log images to a range between 0 and 1
                  imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                  imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                  imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                  img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                  imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                  imgc = reshape(img(numel(imgt1)+[1:numel(imgc1)]),size(imgc1));
                  imgs = reshape(img(numel(imgt1)+numel(imgc1)+[1:numel(imgs1)]),size(imgs1));
                  mn = 0;
                  mx = 1;
                end;
                % recompute min/max for display
                if strcmp(st.vols{i}.window,'auto')
                    mx = -inf; mn = inf;
                end;
                if ~isempty(imgt),
			tmp = imgt(isfinite(imgt));
                        mx = max([mx max(max(tmp))]);
                        mn = min([mn min(min(tmp))]);
                end;
                if ~isempty(imgc),
			tmp = imgc(isfinite(imgc));
                        mx = max([mx max(max(tmp))]);
                        mn = min([mn min(min(tmp))]);
                end;
                if ~isempty(imgs),
			tmp = imgs(isfinite(imgs));
                        mx = max([mx max(max(tmp))]);
                        mn = min([mn min(min(tmp))]);
                end;
                if mx==mn, mx=mn+eps; end;

		if isfield(st.vols{i},'blobs'),
			if ~isfield(st.vols{i}.blobs{1},'colour'),
				% Add blobs for display using the split colourmap
				scal = 64/(mx-mn);
				dcoff = -mn*scal;
				imgt = imgt*scal+dcoff;
				imgc = imgc*scal+dcoff;
				imgs = imgs*scal+dcoff;

				if isfield(st.vols{i}.blobs{1},'max'),
					mx = st.vols{i}.blobs{1}.max;
				else,
					mx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
					st.vols{i}.blobs{1}.max = mx;
				end;
				if isfield(st.vols{i}.blobs{1},'min'),
					mn = st.vols{i}.blobs{1}.min;
				else,
					mn = min([0 minval(st.vols{i}.blobs{1}.vol)]);
					st.vols{i}.blobs{1}.min = mn;
				end;

				vol  = st.vols{i}.blobs{1}.vol;
				M    = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				tmpt = spm_slice_vol(vol,inv(TM0*(st.Space\M)),TD,[0 NaN])';
				tmpc = spm_slice_vol(vol,inv(CM0*(st.Space\M)),CD,[0 NaN])';
				tmps = spm_slice_vol(vol,inv(SM0*(st.Space\M)),SD,[0 NaN])';

				%tmpt_z = find(tmpt==0);tmpt(tmpt_z) = NaN;
				%tmpc_z = find(tmpc==0);tmpc(tmpc_z) = NaN;
				%tmps_z = find(tmps==0);tmps(tmps_z) = NaN;

				sc   = 64/(mx-mn);
				off  = 65.51-mn*sc;
				msk  = find(isfinite(tmpt)); imgt(msk) = off+tmpt(msk)*sc;
				msk  = find(isfinite(tmpc)); imgc(msk) = off+tmpc(msk)*sc;
				msk  = find(isfinite(tmps)); imgs(msk) = off+tmps(msk)*sc;

				cmap = get(st.fig,'Colormap');

                figure(st.fig)
                if mn*mx < 0
                    setcolormap('gray-hot-cold')                    
                elseif mx > 0
                    setcolormap('gray-hot'); % can change this to cold color
                else
                    setcolormap('gray-cold')
                end                
                                redraw_colourbar(i,1,[mn mx],[1:64]'+64); 
			elseif isstruct(st.vols{i}.blobs{1}.colour),
				% Add blobs for display using a defined
                                % colourmap

				% colourmaps
				gryc = [0:63]'*ones(1,3)/63;
				actc = ...
				    st.vols{1}.blobs{1}.colour.cmap;
				actp = ...
				    st.vols{1}.blobs{1}.colour.prop;
				
				% scale grayscale image, not finite -> black
				imgt = scaletocmap(imgt,mn,mx,gryc,65);
				imgc = scaletocmap(imgc,mn,mx,gryc,65);
				imgs = scaletocmap(imgs,mn,mx,gryc,65);
				gryc = [gryc; 0 0 0];
				
				% get max for blob image
				vol = st.vols{i}.blobs{1}.vol;
				mat = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				if isfield(st.vols{i}.blobs{1},'max'),
					cmx = st.vols{i}.blobs{1}.max;
				else,
					cmx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
				end;
				if isfield(st.vols{i}.blobs{1},'min'),
					cmn = st.vols{i}.blobs{1}.min;
				else,
					cmn = -cmx;
				end;

				% get blob data
				vol  = st.vols{i}.blobs{1}.vol;
				M    = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				tmpt = spm_slice_vol(vol,inv(TM0*(st.Space\M)),TD,[0 NaN])';
				tmpc = spm_slice_vol(vol,inv(CM0*(st.Space\M)),CD,[0 NaN])';
				tmps = spm_slice_vol(vol,inv(SM0*(st.Space\M)),SD,[0 NaN])';
				
				% actimg scaled round 0, black NaNs
				topc = size(actc,1)+1;
				tmpt = scaletocmap(tmpt,cmn,cmx,actc,topc);
				tmpc = scaletocmap(tmpc,cmn,cmx,actc,topc);
				tmps = scaletocmap(tmps,cmn,cmx,actc,topc);
				actc = [actc; 0 0 0];
				
				% combine gray and blob data to
				% truecolour
				imgt = reshape(actc(tmpt(:),:)*actp+ ...
					       gryc(imgt(:),:)*(1-actp), ...
					       [size(imgt) 3]);
				imgc = reshape(actc(tmpc(:),:)*actp+ ...
					       gryc(imgc(:),:)*(1-actp), ...
					       [size(imgc) 3]);
				imgs = reshape(actc(tmps(:),:)*actp+ ...
					       gryc(imgs(:),:)*(1-actp), ...
					       [size(imgs) 3]);
				
                                redraw_colourbar(i,1,[cmn cmx],[1:64]'+64); 
				
			else,
				% Add full colour blobs - several sets at once
				scal  = 1/(mx-mn);
				dcoff = -mn*scal;

				wt = zeros(size(imgt));
				wc = zeros(size(imgc));
				ws = zeros(size(imgs));

				imgt  = repmat(imgt*scal+dcoff,[1,1,3]);
				imgc  = repmat(imgc*scal+dcoff,[1,1,3]);
				imgs  = repmat(imgs*scal+dcoff,[1,1,3]);

				cimgt = zeros(size(imgt));
				cimgc = zeros(size(imgc));
				cimgs = zeros(size(imgs));

				for j=1:length(st.vols{i}.blobs), % get colours of all images first
					if isfield(st.vols{i}.blobs{j},'colour'),
						colour(j,:) = reshape(st.vols{i}.blobs{j}.colour, [1 3]);
					else,
						colour(j,:) = [1 0 0];
					end;
				end;
				%colour = colour/max(sum(colour));

				for j=1:length(st.vols{i}.blobs),
					if isfield(st.vols{i}.blobs{j},'max'),
						mx = st.vols{i}.blobs{j}.max;
					else,
						mx = max([eps max(st.vols{i}.blobs{j}.vol(:))]);
						st.vols{i}.blobs{j}.max = mx;
					end;
					if isfield(st.vols{i}.blobs{j},'min'),
						mn = st.vols{i}.blobs{j}.min;
					else,
						mn = min([0 min(st.vols{i}.blobs{j}.vol(:))]);
						st.vols{i}.blobs{j}.min = mn;
					end;

					vol  = st.vols{i}.blobs{j}.vol;
					M    = st.Space\st.vols{i}.premul*st.vols{i}.blobs{j}.mat;
                                        tmpt = spm_slice_vol(vol,inv(TM0*M),TD,[0 NaN])';
                                        tmpc = spm_slice_vol(vol,inv(CM0*M),CD,[0 NaN])';
                                        tmps = spm_slice_vol(vol,inv(SM0*M),SD,[0 NaN])';
                                        % check min/max of sampled image
                                        % against mn/mx as given in st
                                        tmpt(tmpt(:)<mn) = mn;
                                        tmpc(tmpc(:)<mn) = mn;
                                        tmps(tmps(:)<mn) = mn;
                                        tmpt(tmpt(:)>mx) = mx;
                                        tmpc(tmpc(:)>mx) = mx;
                                        tmps(tmps(:)>mx) = mx;
                                        tmpt = (tmpt-mn)/(mx-mn);
					tmpc = (tmpc-mn)/(mx-mn);
					tmps = (tmps-mn)/(mx-mn);
					tmpt(~isfinite(tmpt)) = 0;
					tmpc(~isfinite(tmpc)) = 0;
					tmps(~isfinite(tmps)) = 0;

					cimgt = cimgt + cat(3,tmpt*colour(j,1),tmpt*colour(j,2),tmpt*colour(j,3));
					cimgc = cimgc + cat(3,tmpc*colour(j,1),tmpc*colour(j,2),tmpc*colour(j,3));
					cimgs = cimgs + cat(3,tmps*colour(j,1),tmps*colour(j,2),tmps*colour(j,3));

					wt = wt + tmpt;
					wc = wc + tmpc;
					ws = ws + tmps;
                                        cdata=permute(shiftdim([1/64:1/64:1]'* ...
                                                               colour(j,:),-1),[2 1 3]);
                                        redraw_colourbar(i,j,[mn mx],cdata);
				end;

				imgt = repmat(1-wt,[1 1 3]).*imgt+cimgt;
				imgc = repmat(1-wc,[1 1 3]).*imgc+cimgc;
				imgs = repmat(1-ws,[1 1 3]).*imgs+cimgs;

				imgt(imgt<0)=0; imgt(imgt>1)=1;
				imgc(imgc<0)=0; imgc(imgc>1)=1;
				imgs(imgs<0)=0; imgs(imgs>1)=1;
			end;
		else,
			scal = 64/(mx-mn);
			dcoff = -mn*scal;
			imgt = imgt*scal+dcoff;
			imgc = imgc*scal+dcoff;
			imgs = imgs*scal+dcoff;
		end;

		set(st.vols{i}.ax{1}.d,'HitTest','off', 'Cdata',imgt);
		set(st.vols{i}.ax{1}.lx,'HitTest','off',...
			'Xdata',[0 TD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
		set(st.vols{i}.ax{1}.ly,'HitTest','off',...
			'Ydata',[0 TD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));

		set(st.vols{i}.ax{2}.d,'HitTest','off', 'Cdata',imgc);
		set(st.vols{i}.ax{2}.lx,'HitTest','off',...
			'Xdata',[0 CD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
		set(st.vols{i}.ax{2}.ly,'HitTest','off',...
			'Ydata',[0 CD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));

		set(st.vols{i}.ax{3}.d,'HitTest','off','Cdata',imgs);
		if st.mode ==0,
			set(st.vols{i}.ax{3}.lx,'HitTest','off',...
				'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
			set(st.vols{i}.ax{3}.ly,'HitTest','off',...
				'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(cent(3)-bb(1,3)+1));
		else,
			set(st.vols{i}.ax{3}.lx,'HitTest','off',...
				'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
			set(st.vols{i}.ax{3}.ly,'HitTest','off',...
				'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(bb(2,2)+1-cent(2)));
		end;

		if ~isempty(st.plugins) % process any addons
			for k = 1:prod(size(st.plugins))
				if isfield(st.vols{i},st.plugins{k})
					feval(['spm_ov_', st.plugins{k}], ...
						'redraw', i, TM0, TD, CM0, CD, SM0, SD);
				end;
			end;
		end;
	end;
end;
drawnow;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function redraw_colourbar(vh,bh,interval,cdata)
global st
if isfield(st.vols{vh}.blobs{bh},'cbar')
    if st.mode == 0,
        axpos = get(st.vols{vh}.ax{2}.ax,'Position');
    else,
        axpos = get(st.vols{vh}.ax{1}.ax,'Position');
    end;
    % only scale cdata if we have out-of-range truecolour values
    if ndims(cdata)==3 && max(cdata(:))>1
        cdata=cdata./max(cdata(:));
    end;
    image([0 1],interval,cdata,'Parent',st.vols{vh}.blobs{bh}.cbar);
    set(st.vols{vh}.blobs{bh}.cbar, ...
        'Position',[(axpos(1)+axpos(3)+0.05+(bh-1)*.1)...
                    (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
        'YDir','normal','XTickLabel',[],'XTick',[]);
end;
%_______________________________________________________________________
%_______________________________________________________________________
function centre = findcent
global st
obj    = get(st.fig,'CurrentObject');
centre = [];
cent   = [];
cp     = [];
for i=valid_handles(1:24),
	for j=1:3,
		if ~isempty(obj),
			if (st.vols{i}.ax{j}.ax == obj),
				cp = get(obj,'CurrentPoint');
			end;
		end;
		if ~isempty(cp),
			cp   = cp(1,1:2);
			is   = inv(st.Space);
			cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);
			switch j,
				case 1,
				cent([1 2])=[cp(1)+st.bb(1,1)-1 cp(2)+st.bb(1,2)-1];
				case 2,
				cent([1 3])=[cp(1)+st.bb(1,1)-1 cp(2)+st.bb(1,3)-1];
				case 3,
				if st.mode ==0,
					cent([3 2])=[cp(1)+st.bb(1,3)-1 cp(2)+st.bb(1,2)-1];
				else,
					cent([2 3])=[st.bb(2,2)+1-cp(1) cp(2)+st.bb(1,3)-1];
				end;
			end;
			break;
		end;
	end;
	if ~isempty(cent), break; end;
end;
if ~isempty(cent), centre = st.Space(1:3,1:3)*cent(:) + st.Space(1:3,4); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function handles = valid_handles(handles)
global st;
handles = handles(:)';
handles = handles(find(handles<=24 & handles>=1 & ~rem(handles,1)));
for h=handles,
	if isempty(st.vols{h}), handles(find(handles==h))=[]; end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function reset_st
global st
fig     = spm_figure('FindWin','Graphics');
bb      = []; %[ [-78 78]' [-112 76]' [-50 85]' ];
st      = struct('n', 0, 'vols',[], 'bb',bb,'Space',eye(4),'centre',[0 0 0],'callback',';','xhairs',1,'hld',1,'fig',fig,'mode',1,'plugins',[],'snap',[]);
st.vols = cell(24,1);

pluginpath = fullfile(spm('Dir'),'spm_orthviews');
if isdir(pluginpath)
	pluginfiles = dir(fullfile(pluginpath,'spm_ov_*.m'));
	if ~isempty(pluginfiles)
		addpath(pluginpath);
		% fprintf('spm_orthviews: Using Plugins in %s\n', pluginpath);
		for k = 1:length(pluginfiles)
			[p pluginname e] = fileparts(pluginfiles(k).name);
			st.plugins{k} = strrep(pluginname, 'spm_ov_','');
			% fprintf('%s\n',st.plugins{k});
		end;
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function img = scaletocmap(inpimg,mn,mx,cmap,miscol)
if nargin < 5, miscol=1;end
cml = size(cmap,1);
scf = (cml-1)/(mx-mn);
img = round((inpimg-mn)*scf)+1;
img(find(img<1))   = 1; 
img(find(img>cml)) = cml;
img(~isfinite(img))  = miscol;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function cmap = getcmap(acmapname)
% get colormap of name acmapname
if ~isempty(acmapname),
	cmap = evalin('base',acmapname,'[]');
	if isempty(cmap), % not a matrix, is .mat file?
		[p f e] = fileparts(acmapname);
		acmat   = fullfile(p, [f '.mat']);
		if exist(acmat, 'file'),
			s    = struct2cell(load(acmat));
			cmap = s{1};
		end;
	end;
end;
if size(cmap, 2)~=3,
	warning('Colormap was not an N by 3 matrix')
	cmap = [];
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function item_parent = addcontext(volhandle)
global st;
%create context menu
fg = spm_figure('Findwin','Graphics');set(0,'CurrentFigure',fg);
%contextmenu
item_parent = uicontextmenu;

%contextsubmenu 0
item00  = uimenu(item_parent, 'Label','unknown image', 'Separator','on');
spm_orthviews('context_menu','image_info',item00,volhandle);
item0a    = uimenu(item_parent, 'UserData','pos_mm',     'Callback','spm_orthviews(''context_menu'',''repos_mm'');','Separator','on');
item0b    = uimenu(item_parent, 'UserData','pos_vx',     'Callback','spm_orthviews(''context_menu'',''repos_vx'');');
item0c    = uimenu(item_parent, 'UserData','v_value');

%contextsubmenu 1
item1     = uimenu(item_parent,'Label','Zoom');
item1_1   = uimenu(item1,      'Label','Full Volume',   'Callback','spm_orthviews(''context_menu'',''zoom'',6);', 'Checked','on');
item1_2   = uimenu(item1,      'Label','160x160x160mm', 'Callback','spm_orthviews(''context_menu'',''zoom'',5);');
item1_3   = uimenu(item1,      'Label','80x80x80mm',    'Callback','spm_orthviews(''context_menu'',''zoom'',4);');
item1_4   = uimenu(item1,      'Label','40x40x40mm',    'Callback','spm_orthviews(''context_menu'',''zoom'',3);');
item1_5   = uimenu(item1,      'Label','20x20x20mm',    'Callback','spm_orthviews(''context_menu'',''zoom'',2);');
item1_6   = uimenu(item1,      'Label','10x10x10mm',    'Callback','spm_orthviews(''context_menu'',''zoom'',1);');

%contextsubmenu 2
checked={'off','off'};
checked{st.xhairs+1} = 'on';
item2     = uimenu(item_parent,'Label','Crosshairs');
item2_1   = uimenu(item2,      'Label','on',  'Callback','spm_orthviews(''context_menu'',''Xhair'',''on'');','Checked',checked{2});
item2_2   = uimenu(item2,      'Label','off', 'Callback','spm_orthviews(''context_menu'',''Xhair'',''off'');','Checked',checked{1});

%contextsubmenu 3
if st.Space == eye(4)
	checked = {'off', 'on'};
else
	checked = {'on', 'off'};
end;
item3     = uimenu(item_parent,'Label','Orientation');
item3_1   = uimenu(item3,      'Label','World space', 'Callback','spm_orthviews(''context_menu'',''orientation'',3);','Checked',checked{2});
item3_2   = uimenu(item3,      'Label','Voxel space (1st image)', 'Callback','spm_orthviews(''context_menu'',''orientation'',2);','Checked',checked{1});
item3_3   = uimenu(item3,      'Label','Voxel space (this image)', 'Callback','spm_orthviews(''context_menu'',''orientation'',1);','Checked','off');

%contextsubmenu 3
if isempty(st.snap)
	checked = {'off', 'on'};
else
	checked = {'on', 'off'};
end;
item3     = uimenu(item_parent,'Label','Snap to Grid');
item3_1   = uimenu(item3,      'Label','Don''t snap', 'Callback','spm_orthviews(''context_menu'',''snap'',3);','Checked',checked{2});
item3_2   = uimenu(item3,      'Label','Snap to 1st image', 'Callback','spm_orthviews(''context_menu'',''snap'',2);','Checked',checked{1});
item3_3   = uimenu(item3,      'Label','Snap to this image', 'Callback','spm_orthviews(''context_menu'',''snap'',1);','Checked','off');

%contextsubmenu 4
if st.hld == 0,
	checked = {'off', 'off', 'on'};
elseif st.hld > 0,
	checked = {'off', 'on', 'off'};
else,
	checked = {'on', 'off', 'off'};
end;
item4     = uimenu(item_parent,'Label','Interpolation');
item4_1   = uimenu(item4,      'Label','NN',    'Callback','spm_orthviews(''context_menu'',''interpolation'',3);', 'Checked',checked{3});
item4_2   = uimenu(item4,      'Label','Bilin', 'Callback','spm_orthviews(''context_menu'',''interpolation'',2);','Checked',checked{2});
item4_3   = uimenu(item4,      'Label','Sinc',  'Callback','spm_orthviews(''context_menu'',''interpolation'',1);','Checked',checked{1});

%contextsubmenu 5
% item5     = uimenu(item_parent,'Label','Position', 'Callback','spm_orthviews(''context_menu'',''position'');');

%contextsubmenu 6
item6       = uimenu(item_parent,'Label','Image','Separator','on');
item6_1     = uimenu(item6,      'Label','Window');
item6_1_1   = uimenu(item6_1,    'Label','local');
item6_1_1_1 = uimenu(item6_1_1,  'Label','auto',       'Callback','spm_orthviews(''context_menu'',''window'',2);');
item6_1_1_2 = uimenu(item6_1_1,  'Label','manual',     'Callback','spm_orthviews(''context_menu'',''window'',1);');
item6_1_2   = uimenu(item6_1,    'Label','global');
item6_1_2_1 = uimenu(item6_1_2,  'Label','auto',       'Callback','spm_orthviews(''context_menu'',''window_gl'',2);');
item6_1_2_2 = uimenu(item6_1_2,  'Label','manual',     'Callback','spm_orthviews(''context_menu'',''window_gl'',1);');
if license('test','image_toolbox') == 1
    offon = {'off', 'on'};
    checked = offon(strcmp(st.vols{volhandle}.mapping, ...
                           {'linear', 'histeq', 'loghisteq', 'quadhisteq'})+1);
    item6_2     = uimenu(item6,      'Label','Intensity mapping');
    item6_2_1   = uimenu(item6_2,    'Label','local');
    item6_2_1_1 = uimenu(item6_2_1,  'Label','Linear', 'Checked',checked{1}, ...
                         'Callback','spm_orthviews(''context_menu'',''mapping'',''linear'');');
    item6_2_1_2 = uimenu(item6_2_1,  'Label','Equalised histogram', 'Checked',checked{2}, ...
                         'Callback','spm_orthviews(''context_menu'',''mapping'',''histeq'');');
    item6_2_1_3 = uimenu(item6_2_1,  'Label','Equalised log-histogram', 'Checked',checked{3}, ...
                         'Callback','spm_orthviews(''context_menu'',''mapping'',''loghisteq'');');
    item6_2_1_4 = uimenu(item6_2_1,  'Label','Equalised squared-histogram', 'Checked',checked{4}, ...
                         'Callback','spm_orthviews(''context_menu'',''mapping'',''quadhisteq'');');
    item6_2_2   = uimenu(item6_2,    'Label','global');
    item6_2_2_1 = uimenu(item6_2_2,  'Label','Linear', 'Checked',checked{1}, ...
                         'Callback','spm_orthviews(''context_menu'',''mapping_gl'',''linear'');');
    item6_2_2_2 = uimenu(item6_2_2,  'Label','Equalised histogram', 'Checked',checked{2}, ...
                         'Callback','spm_orthviews(''context_menu'',''mapping_gl'',''histeq'');');
    item6_2_2_3 = uimenu(item6_2_2,  'Label','Equalised log-histogram', 'Checked',checked{3}, ...
                         'Callback','spm_orthviews(''context_menu'',''mapping_gl'',''loghisteq'');');
    item6_2_2_4 = uimenu(item6_2_2,  'Label','Equalised squared-histogram', 'Checked',checked{4}, ...
                         'Callback','spm_orthviews(''context_menu'',''mapping_gl'',''quadhisteq'');');
end;
%contextsubmenu 7
item7     = uimenu(item_parent,'Label','Blobs');
item7_1   = uimenu(item7,      'Label','Add blobs');
item7_1_1 = uimenu(item7_1,    'Label','local',  'Callback','spm_orthviews(''context_menu'',''add_blobs'',2);');
item7_1_2 = uimenu(item7_1,    'Label','global', 'Callback','spm_orthviews(''context_menu'',''add_blobs'',1);');
item7_2   = uimenu(item7,      'Label','Add image');
item7_2_1 = uimenu(item7_2,    'Label','local',  'Callback','spm_orthviews(''context_menu'',''add_image'',2);');
item7_2_2 = uimenu(item7_2,    'Label','global', 'Callback','spm_orthviews(''context_menu'',''add_image'',1);');
item7_3   = uimenu(item7,      'Label','Add colored blobs','Separator','on');
item7_3_1 = uimenu(item7_3,    'Label','local',  'Callback','spm_orthviews(''context_menu'',''add_c_blobs'',2);');
item7_3_2 = uimenu(item7_3,    'Label','global', 'Callback','spm_orthviews(''context_menu'',''add_c_blobs'',1);');
item7_4   = uimenu(item7,      'Label','Add colored image');
item7_4_1 = uimenu(item7_4,    'Label','local',  'Callback','spm_orthviews(''context_menu'',''add_c_image'',2);');
item7_4_2 = uimenu(item7_4,    'Label','global', 'Callback','spm_orthviews(''context_menu'',''add_c_image'',1);');
item7_5   = uimenu(item7,      'Label','Remove blobs',        'Visible','off','Separator','on');
item7_6   = uimenu(item7,      'Label','Remove colored blobs','Visible','off');
item7_6_1 = uimenu(item7_6,    'Label','local', 'Visible','on');
item7_6_2 = uimenu(item7_6,    'Label','global','Visible','on');

if ~isempty(st.plugins) % process any plugins
	for k = 1:prod(size(st.plugins)),
		feval(['spm_ov_', st.plugins{k}], ...
			'context_menu', volhandle, item_parent);
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function c_menu(varargin)
global st

switch lower(varargin{1}),
case 'image_info',
	if nargin <3,
		current_handle = get_current_handle;
	else
		current_handle = varargin{3};
	end;
	if isfield(st.vols{current_handle},'fname'),
		[p,n,e,v] = spm_fileparts(st.vols{current_handle}.fname);
                if isfield(st.vols{current_handle},'n')
                    v = sprintf(',%d',st.vols{current_handle}.n);
                end;
		set(varargin{2}, 'Label',[n e v]);
	end;
	delete(get(varargin{2},'children'));
	if exist('p','var')
		item1 = uimenu(varargin{2}, 'Label', p);
	end;
	if isfield(st.vols{current_handle},'descrip'),
		item2 = uimenu(varargin{2}, 'Label',...
		st.vols{current_handle}.descrip);
	end;
	dt = st.vols{current_handle}.dt(1);
	item3 = uimenu(varargin{2}, 'Label', sprintf('Data type: %s', spm_type(dt)));
	str   = 'Intensity: varied';
	if size(st.vols{current_handle}.pinfo,2) == 1,
		if st.vols{current_handle}.pinfo(2),
			str = sprintf('Intensity: Y = %g X + %g',...
				st.vols{current_handle}.pinfo(1:2)');
		else,
			str = sprintf('Intensity: Y = %g X', st.vols{current_handle}.pinfo(1)');
		end;
	end;
	item4  = uimenu(varargin{2}, 'Label',str);
	item5  = uimenu(varargin{2}, 'Label', 'Image dims', 'Separator','on');
	item51 = uimenu(varargin{2}, 'Label',...
		sprintf('%dx%dx%d', st.vols{current_handle}.dim(1:3)));
	prms   = spm_imatrix(st.vols{current_handle}.mat);
	item6  = uimenu(varargin{2}, 'Label','Voxel size', 'Separator','on');
	item61 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', prms(7:9)));
	item7  = uimenu(varargin{2}, 'Label','Origin', 'Separator','on');
	item71 = uimenu(varargin{2}, 'Label',...
		sprintf('%.2f %.2f %.2f', prms(1:3)));
	R      = spm_matrix([0 0 0 prms(4:6)]);
	item8  = uimenu(varargin{2}, 'Label','Rotations', 'Separator','on');
	item81 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(1,1:3)));
	item82 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(2,1:3)));
	item83 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(3,1:3)));
	item9  = uimenu(varargin{2},...
		'Label','Specify other image...',...
		'Callback','spm_orthviews(''context_menu'',''swap_img'');',...
		'Separator','on');

case 'repos_mm',
	oldpos_mm = spm_orthviews('pos');
	newpos_mm = spm_input('New Position (mm)','+1','r',sprintf('%.2f %.2f %.2f',oldpos_mm),3);
	spm_orthviews('reposition',newpos_mm);

case 'repos_vx'
	current_handle = get_current_handle;
	oldpos_vx = spm_orthviews('pos', current_handle);
	newpos_vx = spm_input('New Position (voxels)','+1','r',sprintf('%.2f %.2f %.2f',oldpos_vx),3);
	newpos_mm = st.vols{current_handle}.mat*[newpos_vx;1];
	spm_orthviews('reposition',newpos_mm(1:3));

case 'zoom'
	zoom_all(varargin{2});
	bbox;
	redraw_all;

case 'xhair',
	spm_orthviews('Xhairs',varargin{2});
	cm_handles = get_cm_handles;
	for i = 1:length(cm_handles),
		z_handle = get(findobj(cm_handles(i),'label','Crosshairs'),'Children');
		set(z_handle,'Checked','off'); %reset check
		if strcmp(varargin{2},'off'), op = 1; else op = 2; end
		set(z_handle(op),'Checked','on');
	end;

case 'orientation',
	cm_handles = get_cm_handles;
	for i = 1:length(cm_handles),
		z_handle = get(findobj(cm_handles(i),'label','Orientation'),'Children');
		set(z_handle,'Checked','off');
	end;
	if varargin{2} == 3,
		spm_orthviews('Space');
	elseif varargin{2} == 2,
		spm_orthviews('Space',1);
	else,
		spm_orthviews('Space',get_current_handle);
		z_handle = get(findobj(st.vols{get_current_handle}.ax{1}.cm,'label','Orientation'),'Children');
		set(z_handle(1),'Checked','on');
		return;
	end;
	for i = 1:length(cm_handles),
		z_handle = get(findobj(cm_handles(i),'label','Orientation'),'Children');
		set(z_handle(varargin{2}),'Checked','on');
	end;

case 'snap',
	cm_handles = get_cm_handles;
	for i = 1:length(cm_handles),
		z_handle = get(findobj(cm_handles(i),'label','Snap to Grid'),'Children');
		set(z_handle,'Checked','off');
	end;
	if varargin{2} == 3,
		st.snap = [];
	elseif varargin{2} == 2,
		st.snap = 1;
	else,
		st.snap = get_current_handle;
		z_handle = get(findobj(st.vols{get_current_handle}.ax{1}.cm,'label','Snap to Grid'),'Children');
		set(z_handle(1),'Checked','on');
		return;
	end;
	for i = 1:length(cm_handles),
		z_handle = get(findobj(cm_handles(i),'label','Snap to Grid'),'Children');
		set(z_handle(varargin{2}),'Checked','on');
	end;

case 'interpolation',
	tmp        = [-4 1 0];
	st.hld     = tmp(varargin{2});
	cm_handles = get_cm_handles;
	for i = 1:length(cm_handles),
		z_handle = get(findobj(cm_handles(i),'label','Interpolation'),'Children');
		set(z_handle,'Checked','off');
		set(z_handle(varargin{2}),'Checked','on');
	end;
	redraw_all;

case 'window',
	current_handle = get_current_handle;
	if varargin{2} == 2,
		spm_orthviews('window',current_handle);
	else
		if isnumeric(st.vols{current_handle}.window)
			defstr = sprintf('%.2f %.2f', st.vols{current_handle}.window);
		else
			defstr = '';
		end;
		spm_orthviews('window',current_handle,spm_input('Range','+1','e',defstr,2));
	end;

case 'window_gl',
	if varargin{2} == 2,
		for i = 1:length(get_cm_handles),
			st.vols{i}.window = 'auto';
		end;
	else,
		current_handle = get_current_handle;
		if isnumeric(st.vols{current_handle}.window)
			defstr = sprintf('%d %d', st.vols{current_handle}.window);
		else
			defstr = '';
		end;
		data = spm_input('Range','+1','e',defstr,2);

		for i = 1:length(get_cm_handles),
			st.vols{i}.window = data;
		end;
	end;
	redraw_all;
        
case 'mapping',
        checked = strcmp(varargin{2}, ...
                         {'linear', 'histeq', 'loghisteq', ...
                          'quadhisteq'});
        checked = checked(end:-1:1); % Handles are stored in inverse order
	current_handle = get_current_handle;        
        cm_handles = get_cm_handles;
        st.vols{current_handle}.mapping = varargin{2};
        z_handle = get(findobj(cm_handles(current_handle), ...
                               'label','Intensity mapping'),'Children');
        for k = 1:numel(z_handle)
                c_handle = get(z_handle(k), 'Children');
                set(c_handle, 'checked', 'off');
                set(c_handle(checked), 'checked', 'on');
        end;
        redraw_all;
        
case 'mapping_gl',
        checked = strcmp(varargin{2}, ...
                         {'linear', 'histeq', 'loghisteq', 'quadhisteq'});
        checked = checked(end:-1:1); % Handles are stored in inverse order
        cm_handles = get_cm_handles;
        for k = valid_handles(1:24),
                st.vols{k}.mapping = varargin{2};
                z_handle = get(findobj(cm_handles(k), ...
                                       'label','Intensity mapping'),'Children');
                for l = 1:numel(z_handle)
                        c_handle = get(z_handle(l), 'Children');
                        set(c_handle, 'checked', 'off');
                        set(c_handle(checked), 'checked', 'on');
                end;
        end;
        redraw_all;
        
case 'swap_img',
	current_handle = get_current_handle;
	new_info = spm_vol(spm_select(1,'image','select new image'));
        fn = fieldnames(new_info);
        for k=1:numel(fn)
                st.vols{current_handle}.(fn{k}) = new_info.(fn{k});
        end;
	spm_orthviews('context_menu','image_info',get(gcbo, 'parent'));
	redraw_all;

case 'add_blobs',
	% Add blobs to the image - in split colortable
	cm_handles = valid_handles(1:24);
	if varargin{2} == 2, cm_handles = get_current_handle; end;
	spm_figure('Clear','Interactive');
	[SPM,VOL] = spm_getSPM;
	for i = 1:length(cm_handles),
		addblobs(cm_handles(i),VOL.XYZ,VOL.Z,VOL.M);
		c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove blobs');
		set(c_handle,'Visible','on');
		delete(get(c_handle,'Children'));
		item7_3_1 = uimenu(c_handle,'Label','local','Callback','spm_orthviews(''context_menu'',''remove_blobs'',2);');
		if varargin{2} == 1,
			item7_3_2 = uimenu(c_handle,'Label','global','Callback','spm_orthviews(''context_menu'',''remove_blobs'',1);');
		end;
	end;
	redraw_all;

case 'remove_blobs',
	cm_handles = valid_handles(1:24);
	if varargin{2} == 2, cm_handles = get_current_handle; end;
	for i = 1:length(cm_handles),
		rmblobs(cm_handles(i));
		c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove blobs');
		delete(get(c_handle,'Children'));
		set(c_handle,'Visible','off');
	end;
	redraw_all;

case 'add_image',
	% Add blobs to the image - in split colortable
	cm_handles = valid_handles(1:24);
	if varargin{2} == 2, cm_handles = get_current_handle; end;
	spm_figure('Clear','Interactive');
	fname = spm_select(1,'image','select image');
	for i = 1:length(cm_handles),
		addimage(cm_handles(i),fname);
		c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove blobs');
		set(c_handle,'Visible','on');
		delete(get(c_handle,'Children'));
		item7_3_1 = uimenu(c_handle,'Label','local','Callback','spm_orthviews(''context_menu'',''remove_blobs'',2);');
		if varargin{2} == 1,
			item7_3_2 = uimenu(c_handle,'Label','global','Callback','spm_orthviews(''context_menu'',''remove_blobs'',1);');
		end;
	end;
	redraw_all;

case 'add_c_blobs',
	% Add blobs to the image - in full colour
	cm_handles = valid_handles(1:24);
	if varargin{2} == 2, cm_handles = get_current_handle; end;
	spm_figure('Clear','Interactive');
	[SPM,VOL] = spm_getSPM;
	c         = spm_input('Colour','+1','m',...
		'Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
	colours   = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
	c_names   = {'red';'yellow';'green';'cyan';'blue';'magenta'};
        hlabel = sprintf('%s (%s)',VOL.title,c_names{c});
	for i = 1:length(cm_handles),
		addcolouredblobs(cm_handles(i),VOL.XYZ,VOL.Z,VOL.M,colours(c,:));
                addcolourbar(cm_handles(i),numel(st.vols{cm_handles(i)}.blobs));
		c_handle    = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove colored blobs');
		ch_c_handle = get(c_handle,'Children');
		set(c_handle,'Visible','on');
		%set(ch_c_handle,'Visible',on');
		item7_4_1   = uimenu(ch_c_handle(2),'Label',hlabel,'ForegroundColor',colours(c,:),...
			'Callback','c = get(gcbo,''UserData'');spm_orthviews(''context_menu'',''remove_c_blobs'',2,c);',...
			'UserData',c);
		if varargin{2} == 1,
			item7_4_2 = uimenu(ch_c_handle(1),'Label',hlabel,'ForegroundColor',colours(c,:),...
				'Callback','c = get(gcbo,''UserData'');spm_orthviews(''context_menu'',''remove_c_blobs'',1,c);',...
				'UserData',c);
		end;
	end;
	redraw_all;

case 'remove_c_blobs',
    cm_handles = valid_handles(1:24);
    if varargin{2} == 2, cm_handles = get_current_handle; end;
    colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
    c_names = {'red';'yellow';'green';'cyan';'blue';'magenta'};
    for i = 1:length(cm_handles),
        if isfield(st.vols{cm_handles(i)},'blobs'),
            for j = 1:length(st.vols{cm_handles(i)}.blobs),
                if st.vols{cm_handles(i)}.blobs{j}.colour == colours(varargin{3},:);
                    if isfield(st.vols{cm_handles(i)}.blobs{j},'cbar')
                        delete(st.vols{cm_handles(i)}.blobs{j}.cbar);
                    end
                    st.vols{cm_handles(i)}.blobs(j) = [];
                    break;
                end;
            end;
            rm_c_menu = findobj(st.vols{cm_handles(i)}.ax{1}.cm,'Label','Remove colored blobs');
            delete(findobj(rm_c_menu,'Label',c_names{varargin{3}}));
            if isempty(st.vols{cm_handles(i)}.blobs),
                st.vols{cm_handles(i)} = rmfield(st.vols{cm_handles(i)},'blobs');
                set(rm_c_menu, 'Visible', 'off');
            end;
        end;
    end;
    redraw_all;

case 'add_c_image',
	% Add truecolored image
	cm_handles = valid_handles(1:24);
	if varargin{2} == 2, cm_handles = get_current_handle;end;
	spm_figure('Clear','Interactive');
	fname   = spm_select(1,'image','select image');
	c       = spm_input('Colour','+1','m','Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
	colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
	c_names = {'red';'yellow';'green';'cyan';'blue';'magenta'};
        hlabel = sprintf('%s (%s)',fname,c_names{c});
	for i = 1:length(cm_handles),
		addcolouredimage(cm_handles(i),fname,colours(c,:));
                addcolourbar(cm_handles(i),numel(st.vols{cm_handles(i)}.blobs));
		c_handle    = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove colored blobs');
		ch_c_handle = get(c_handle,'Children');
		set(c_handle,'Visible','on');
		%set(ch_c_handle,'Visible',on');
		item7_4_1 = uimenu(ch_c_handle(2),'Label',hlabel,'ForegroundColor',colours(c,:),...
			'Callback','c = get(gcbo,''UserData'');spm_orthviews(''context_menu'',''remove_c_blobs'',2,c);','UserData',c);
		if varargin{2} == 1
			item7_4_2 = uimenu(ch_c_handle(1),'Label',hlabel,'ForegroundColor',colours(c,:),...
				'Callback','c = get(gcbo,''UserData'');spm_orthviews(''context_menu'',''remove_c_blobs'',1,c);',...
				'UserData',c);
		end
	end
	redraw_all;
end;
%_______________________________________________________________________
%_______________________________________________________________________
function current_handle = get_current_handle
global st
cm_handle      = get(gca,'UIContextMenu');
cm_handles     = get_cm_handles;
current_handle = find(cm_handles==cm_handle);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function cm_pos
global st
for i = 1:length(valid_handles(1:24)),
	if isfield(st.vols{i}.ax{1},'cm')
		set(findobj(st.vols{i}.ax{1}.cm,'UserData','pos_mm'),...
			'Label',sprintf('mm:  %.1f %.1f %.1f',spm_orthviews('pos')));
		pos = spm_orthviews('pos',i);
		set(findobj(st.vols{i}.ax{1}.cm,'UserData','pos_vx'),...
			'Label',sprintf('vx:  %.1f %.1f %.1f',pos));
		set(findobj(st.vols{i}.ax{1}.cm,'UserData','v_value'),...
			'Label',sprintf('Y = %g',spm_sample_vol(st.vols{i},pos(1),pos(2),pos(3),st.hld)));
	end
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function cm_handles = get_cm_handles
global st
cm_handles = [];
for i=valid_handles(1:24),
	cm_handles = [cm_handles st.vols{i}.ax{1}.cm];
end
return;
%_______________________________________________________________________
%_______________________________________________________________________
function zoom_all(op)
global st
cm_handles = get_cm_handles;
res = [.125 .125 .25 .5 1 1];
if op==6,
	st.bb = maxbb;
else,
	vx = sqrt(sum(st.Space(1:3,1:3).^2));
	vx = vx.^(-1);
	pos = spm_orthviews('pos');
	pos = st.Space\[pos ; 1];
	pos = pos(1:3)';
	if     op == 5, st.bb = [pos-80*vx ; pos+80*vx] ;
	elseif op == 4, st.bb = [pos-40*vx ; pos+40*vx] ;
	elseif op == 3, st.bb = [pos-20*vx ; pos+20*vx] ;
	elseif op == 2, st.bb = [pos-10*vx ; pos+10*vx] ;
	elseif op == 1; st.bb = [pos- 5*vx ; pos+ 5*vx] ;
	else disp('no Zoom possible');
	end;
end
resolution(res(op));
redraw_all;
for i = 1:length(cm_handles)
	z_handle = get(findobj(cm_handles(i),'label','Zoom'),'Children');
	set(z_handle,'Checked','off');
	set(z_handle(op),'Checked','on');
end
return;


function display_slices(imgs, colormaps, ranges)
% FORMAT display_slices(imgs, dispf)
%
% request a some parameters for slice_overlay routine
% while accepting many defaults
% SO structure contains all the parameters for display
% See slice_overlay.m for detailed comments
%
% Defaults to GUI if no arguments passed
% imgs  - string or cell array of image names to display
% dispf - flag, if set, displays overlay (default = 1)
%
% Matthew Brett 5/00 V 0.3
% revised by Xu Cui 1/4/2008
% 
% The first image file has to be a structural image
 
if ischar(imgs)
  imgs = cellstr(imgs);
end
  
clear global SO
global SO  

% load images
nimgs = size(imgs);

% process names
nchars = 20;
imgns = spm_str_manip(imgs, ['rck' num2str(nchars)]);

% identify image types
cscale = [];
deftype = 1;
SO.cbar = [];
for i = 1:nimgs
  SO.img(i).vol = spm_vol(imgs{i});
  options = {'Structural','Truecolour', ...
	     'Blobs','Negative blobs'};
  % if there are SPM results in the workspace, add this option
  if evalin('base','exist(''SPM'', ''var'')')
    options = {'Structural with SPM blobs', options{:}};
  end
%   itype = spm_input(sprintf('Img %d: %s - image type?', i, imgns{i}), '+1', ...
% 		    'm', char(options),options, deftype);
  if i == 1
      itype = {'Structural'};
  else
      itype = {'Blobs'};
  end        
  %imgns(i) = {sprintf('Img %d (%s)',i,itype{1})};
  [mx mn] = slice_overlay('volmaxmin', SO.img(i).vol);
  if ~isempty(strmatch('Structural', itype))
    SO.img(i).cmap = gray;
    SO.img(i).range = [mn mx];
    deftype = 2;
    cscale = [cscale i];
    if strcmp(itype,'Structural with SPM blobs')
      errstr = sprintf(['Cannot find SPM/VOL structs in the workspace\n'...
		       'Please run SPM results GUI before' ...
			' display_slices']);
      SPM = evalin('base', 'SPM', ['error(' errstr ')']);
      VOL = evalin('base', 'VOL', ['error(' errstr ')']);
      slice_overlay('addspm',SPM,VOL,0);
    end
  else
    SO.cbar = [SO.cbar i];
    cprompt = ['Colormap: ' imgns{i}];
    switch itype{1}
     case 'Truecolour'
      dcmap = 'actc';
      drange = [mn mx];
      cscale = [cscale i];
     case 'Blobs'
      dcmap = 'hot';
      drange = [0 mx];
      SO.img(i).prop = Inf;
     case 'Negative blobs'
      dcmap = 'winter';
      drange = [0 mn];
      SO.img(i).prop = Inf;
    end
    if isempty(colormaps{i})
        SO.img(i).cmap = return_cmap(cprompt, dcmap);
    else
        SO.img(i).cmap = colormaps{i};%return_cmap(cprompt, dcmap);
    end
    SO.img(i).range = ranges{i};%spm_input('Img val range for colormap','+1', 'e', drange, 2);
  end
end
ncmaps=length(cscale);
if ncmaps == 1
  SO.img(cscale).prop = 1;
else
  remcol=1;
  for i = 1:ncmaps
    ino = cscale(i);
    SO.img(ino).prop = spm_input(sprintf('%s intensity',imgns{ino}),...
				 '+1', 'e', ...
				 remcol/(ncmaps-i+1),1);
    remcol = remcol - SO.img(ino).prop;
  end
end
%  
% SO.transform = deblank(spm_input('Image orientation', '+1', ['Axial|' ...
% 		    ' Coronal|Sagittal'], strvcat('axial','coronal','sagittal'), ...
% 		    1));
SO.transform = deblank(questdlg('Image orientation?', ...
                         'Image orientation', ...
                         'Axial', 'Coronal', 'Sagittal', 'Axial'));        

% use SPM figure window
SO.figure = spm_figure('GetWin', 'Graphics'); 

% slices for display
slice_overlay('checkso');
SO.slices = 		      eval(sprintf('%0.0f:%0.0f:%0.0f',...
			      SO.slices(1),...
			      mean(diff(SO.slices)),...
			      SO.slices(end)));

% and do the display

slice_overlay

return

function cmap = return_cmap(prompt,defmapn)

prompt={prompt};
name='Input color';
numlines=1;
defaultanswer={defmapn};

answer=inputdlg(prompt,name,numlines,defaultanswer);
answer = answer{1};

cmap = [];
while isempty(cmap)
  cmap = slice_overlay('getcmap', answer);
end
return

% slice_overlay() - Function to display + manage slice display
% 
% Usage: 
%    >> slice_overlay(action, varargin);
%
% Inputs:
%
% Outputs:
%
% Author: Matthew Brett matthew@mrc-cbu.cam.ac.uk
%
% See also: show_2d(), show_3d()

% $Log: slice_overlay.m,v $
% Revision 1.1  2003/02/06 19:17:48  duann
% Initial revision
%

function varargout = slice_overlay(action, varargin);
% Function to display + manage slice display 
% Slice display works on a global structure SO
% with fields 
%  - img - array of images to display
%        - img structs contain fields
%             vol - vol struct info (see spm_vol)
%                   can also be vol containing image as 3d matrix
%                   set with slice_overlay('AddBlobs'...) call
%             cmap - colormap for this image
%             nancol - color for NaN. If scalar, this is an index into
%                    the image cmap.  If 1x3 vector, it's a colour
%             prop - proportion of intensity for this cmap/img 
%                    if = Inf, gives split cmap effect where values of 
%                    this cmap override previous image cmap values
%             func - function to apply to image before scaling to cmap
%                    (and therefore before min/max thresholding. E.g. a func of
%                    'i1(i1==0)=NaN' would convert zeros to NaNs
%             range - 2x1 vector of values for image to distribute colormap across
%                    the first row of the colormap applies to the first
%                    value in 'range', and the last value to the second
%                    value in 'range'
%             outofrange - behavior for image values to the left and
%                    right of image limits in 'range'.  Left means
%                    colormap values < 1, i.e for image values <
%                    range(1), if (range(1)<range(2)), and image values >
%                    range(1) where (range(1)>range(2)). If missing,
%                    display min (for Left) and max (for Right) value from colormap. 
%                    Otherwise should be a 2 element cell array, where
%                    the first element is the colour value for image values
%                    left of 'range', and the second is for image values
%                    right of 'range'.  Scalar values for
%                    colour index the colormap, 3x1 vectors are colour
%                    values.  An empty array attracts default settings
%                    appropriate to the mode - i.e. transparent colour (where 
%                    SO.prop ~= Inf), or split colour.  Empty cells
%                    default to 0. 0 specifies that voxels with this
%                    colour do not influence the image (split =
%                    background, true = black)
%            hold  - resampling order for image (see spm_sample_vol) -
%                    default 1
%            background - value when resampling outside image - default
%                    NaN
%            
% - transform - either - 4x4 transformation to apply to image slice position,
%             relative to mm given by slicedef, before display
%               or     - text string, one of axial, coronal, sagittal
%                        These orientations assume the image is currently
%                        (after its mat file has been applied) axially
%                        oriented
% - slicedef - 2x3 array specifying dimensions for slice images in mm
%             where rows are x,and y of slice image, and cols are neg max dim,
%             slice separation and pos max dim
% - slices   - vector of slice positions in mm in z (of transformed image)
% - figure    - figure handle for slice display figure
% - refreshf  - flag - if set or empty, refresh axis info for figure
%             else assume this is OK
% - clf       - flag, non zero -> clear figure before display.  Redundant
%               if refreshf == 0
% - area      struct with fields
%                  position - bottom left, x size y size 1x4 vector of
%                      area in which to display slices
%                  units    - one of
%                    inches,centimeters,normalized,points,{pixels}
%                  halign - one of left,{center},right
%                  valign - one of top,{middle},bottom
% - xslices  - no of slices to display across figure (defaults to an optimum)
% - cbar      - if empty, missing, no colourbar.  If an array of integers, then 
%             indexes img array, and makes colourbar for each cmap for
%             that img.  Cbars specified in order of appearance L->R
% - labels - struct can be absent (-> default numerical labels)
%                  empty (SO.labels = []) (no labels) or contain fields 
%                  colour - colour for label text 
%                  size - font size in units normalized to slice axes 
%                  format - if = cell array of strings =
%                  labels for each slice in Z.  If is string, specifies
%                  sprintf format string for labelling in distance of the
%                  origin (Xmm=0, Ymm=0) of each slice from plane containing
%                  the AC, in mm, in the space of the transformed image
% - callback - callback string for button down on image panels.  E.g.
%              setting SO.callback to 'slice_overlay(''getpos'')' prints to
%              the matlab window the equivalent position in mm of the
%              position of a mouse click on one of the image slices
% - printstr - string for printing slice overlay figure window, e.g.
%              'print -dpsc -painters -noui' (the default)
% - printfile - name of file to print output to; default 'slices.ps'
%
% FORMAT slice_overlay
% Checks, fills SO struct (slice_overlay('checkso')), and 
% displays slice overlay (slice_overlay('display'))
%
% FORMAT slice_overlay('checkso')
% Checks SO structure and sets defaults
%
% FORMAT cmap = slice_overlay('getcmap',cmapname)
% Gets colormap named in cmapname string
%
% FORMAT [mx mn] = slice_overlay('volmaxmin', vol)
% Returns maximum and minimum finite values from vol struct 'vol'
%
% FORMAT slice_overlay('addspm',SPM,VOL,dispf)
% Adds SPM blobs as new img to SO struct, split effect, 'hot' colormap, 
% Structures SPM and VOL are generated by calls to SPM results
% if not passed, they are fetched from the workspace
% If dispf is not passed, or nonzero, displays resulting SO figure also
%
% FORMAT slice_overlay('addblobs', imgno, XYZ, vals, mat)
% adds SPM blobs to img no 'imgno', as specified in 
% XYZ  - 3xN voxel coordinates of N blob values
% vals - N blob intensity values
% mat  - 4x4 matrix specifying voxels -> mm
%
% FORMAT vol = slice_overlay('blobs2vol', XYZ, vals, mat)
% returns (pseudo) vol struct for 3d blob volume specified
% in matrices as above
%
% FORMAT slice_overlay('addmatrix', imgno, mat3d, mat)
% adds 3d matrix image vol to img imgno.  Optionally
% mat  - 4x4 matrix specifying voxels -> mm
%
% FORMAT vol = slice_overlay('matrix2vol', mat3d, mat)
% returns (pseudo) vol struct for 3d matrix 
% input matrices as above
%
% FORMAT mmpos = slice_overlay('getpos')
% returns equivalent position in mm of last click on current axes (gca)
% if the axes contain an image slice (empty otherwise)
%
% FORMAT vals = slice_overlay('pointvals', XYZmm, holdlist)
% returns IxN matrix with values of each image 1..I, at each
% point 1..N specified in 3xN mm coordinate matrix XYZmm
% If specified, 'holdlist' contains I values giving hold
% values for resampling for each image (see spm_sample_vol)
%
% FORMAT slice_overlay('display')
% Displays slice overlay from SO struct
% 
% FORMAT slice_overlay('print', filename, printstr) 
% Prints slice overlay figure, usually to file.  If filename is not
% passed/empty gets filename from SO.printfile.  If printstr is not
% passed/empty gets printstr from SO.printstr
% 
% V 0.8 2/8/00  
% More or less  beta - take care.  Please report problems to  
% Matthew Brett matthew@mrc-cbu.cam.ac.uk

global SO

if nargin < 1
  checkso;
  action = 'display';
else
  action = lower(action);
end

switch action
 case 'checkso'
  checkso;
 case 'getcmap'
  varargout = {getcmap(varargin{1})};
 case 'volmaxmin'
  [mx mn] = volmaxmin(varargin{1});
  varargout = {mx, mn};
 case 'addspm'
  if nargin < 2
    varargin(1) = {evalin('base', 'SPM', ['error(''Cannot find SPM' ...
		    ' struct'')'])};
  end
  if nargin < 3
    varargin(2) = {evalin('base', 'VOL', ['error(''Cannot find VOL' ...
		    ' struct'')'])};
  end
  if nargin < 4
    varargin{3} = 1;
  end
  newimg = length(SO.img)+1;
  SO.img(newimg).vol = blobs2vol(varargin{1}.XYZ,varargin{1}.Z, varargin{2}.M);
  SO.img(newimg).prop = Inf;
  SO.img(newimg).cmap = hot;
  SO.img(newimg).range = [0 max(varargin{1}.Z)];
  SO.cbar = [SO.cbar newimg];
  if varargin{3}
    checkso;
    slice_overlay('display');
  end
  
 case 'addblobs'
  addblobs(varargin{1},varargin{2},varargin{3},varargin{4});
 case 'blobs2vol'
  varargout = {blobs2vol(varargin{1},varargin{2},varargin{3})};
 case 'addmatrix'
  if nargin<3,varargin{2}='';end
  if nargin<4,varargin{3}='';end
  addmatrix(varargin{1},varargin{2},varargin{3});
 case 'matrix2vol'
  if nargin<3,varargin{2}=[];end
  varargout = {matrix2vol(varargin{1},varargin{2})};
 case 'getpos'
  varargout = {getpos};
 case 'pointvals'
  varargout = {pointvals(varargin{1})};
 case 'print'
  if nargin<2,varargin{1}='';end
  if nargin<3,varargin{2}='';end
  printfig(varargin{1}, varargin{2});
 case 'display'

% get coordinates for plane
X=1;Y=2;Z=3;
dims = SO.slicedef;
xmm = dims(X,1):dims(X,2):dims(X,3);
ymm = dims(Y,1):dims(Y,2):dims(Y,3);
zmm = SO.slices;
[y x] = meshgrid(ymm,xmm');
vdims = [length(xmm),length(ymm),length(zmm)];

% no of slices, and panels (an extra for colorbars)
nslices = vdims(Z);
minnpanels = nslices;
cbars = 0;
if is_there(SO,'cbar')
  cbars = length(SO.cbar);
  minnpanels = minnpanels+cbars;
end

% get figure data
% if written to, the axes may be specified already
figno = figure(SO.figure);

% (re)initialize axes and stuff

% check if the figure is set up correctly
if ~SO.refreshf
  axisd = flipud(findobj(SO.figure, 'Type','axes','Tag', 'slice overlay panel'));
  npanels = length(axisd);
  if npanels < vdims(Z)+cbars;
    SO.refreshf = 1;
  end
end
if SO.refreshf
  % clear figure, axis store
  if SO.clf, clf; end
  axisd = [];

  % prevent print inversion problems
  set(figno,'InvertHardCopy','off');
  
  % calculate area of display in pixels
  parea = SO.area.position;
  if ~strcmp(SO.area.units, 'pixels')
    ubu = get(SO.figure, 'units');
    set(SO.figure, 'units','pixels');
    tmp = get(SO.figure, 'Position');
    ascf = tmp(3:4);
    if ~strcmp(SO.area.units, 'normalized')
      set(SO.figure, 'units',SO.area.units);
      tmp = get(SO.figure, 'Position');
      ascf = ascf ./ tmp(3:4);
    end
    set(figno, 'Units', ubu);
    parea = parea .* repmat(ascf, 1, 2);
  end
  asz = parea(3:4);
  
  % by default, make most parsimonious fit to figure
  yxratio = length(ymm)*dims(Y,2)/(length(xmm)*dims(X,2));
  if ~is_there(SO, 'xslices')
    % iteration needed to optimize, surprisingly.  Thanks to Ian NS
    axlen(X,:)=asz(1):-1:1;
    axlen(Y,:)=yxratio*axlen(X,:);
    panels = floor(asz'*ones(1,size(axlen,2))./axlen);
    estnpanels = prod(panels);
    tmp = find(estnpanels >= minnpanels);
    if isempty(tmp)
      error('Whoops, cannot fit panels onto figure');
    end
    b = tmp(1); % best fitting scaling
    panels = panels(:,b);
    axlen = axlen(:, b);
  else
    % if xslices is specified, assume X is flush with X figure dimensions
    panels([X:Y],1) = [SO.xslices; 0];
    axlen([X:Y],1) = [asz(X)/panels(X); 0];
  end
  
  % Axis dimensions are in pixels.  This prevents aspect ratio rescaling
  panels(Y) = ceil(minnpanels/panels(X));
  axlen(Y) = axlen(X)*yxratio;
  
  % centre (etc) panels in display area as required
  divs = [Inf 2 1];the_ds = [0;0];
  the_ds(X) = divs(strcmp(SO.area.halign, {'left','center','right'}));
  the_ds(Y) = divs(strcmp(SO.area.valign, {'bottom','middle','top'}));
  startc = parea(1:2)' + (asz'-(axlen.*panels))./the_ds;
  
  % make axes for panels
  r=0;c=1;
  npanels = prod(panels);
  lastempty = npanels-cbars;
  for i = 1:npanels
    % panel userdata
    if i<=nslices
      u.type = 'slice';
      u.no   = zmm(i);
    elseif i > lastempty
      u.type = 'cbar';
      u.no   = i - lastempty;
    else
      u.type = 'empty';
      u.no   = i - nslices;
    end
    axpos = [r*axlen(X)+startc(X) (panels(Y)-c)*axlen(Y)+startc(Y) axlen'];
    axisd(i) = axes(...
	'Parent',figno,...
	'XTick',[],...
	'XTickLabel',[],...
	'YTick',[],...
	'YTickLabel',[],...
	'Box','on',...
	'XLim',[1 vdims(X)],...
	'YLim',[1 vdims(Y)],...
	'Units', 'pixels',...
	'Position',axpos,...
	'Tag','slice overlay panel',...
	'UserData',u);
    r = r+1;
    if r >= panels(X)
      r = 0;
      c = c+1;
    end
  end
end

% sort out labels
if is_there(SO,'labels')
  labels = SO.labels;
  if iscell(labels.format)
    if length(labels.format)~=vdims(Z)
      error(...
	  sprintf('Oh dear, expecting %d labels, but found %d',...
		  vdims(Z), length(labels.contents)));
    end
  else
    % format string for mm from AC labelling
    fstr = labels.format;
    labels.format = cell(vdims(Z),1);
    acpt = SO.transform * [0 0 0 1]';
    for i = 1:vdims(Z)
      labels.format(i) = {sprintf(fstr,zmm(i)-acpt(Z))};
    end
  end
end

% modify colormaps with any new colours
nimgs = length(SO.img);
lrn = zeros(nimgs,3);
cmaps = cell(nimgs);
for i = 1:nimgs
  cmaps(i)={SO.img(i).cmap};
  lrnv = {SO.img(i).outofrange{:}, SO.img(i).nancol};
  for j = 1:length(lrnv)
    if prod(size(lrnv{j}))==1
      lrn(i,j) = lrnv{j};
    else
      cmaps(i) = {[cmaps{i}; lrnv{j}(1:3)]};
      lrn(i,j) = size(cmaps{i},1);
    end
  end
end

% cycle through slices displaying images
nvox = prod(vdims(1:2));
pandims = [vdims([2 1]) 3]; % NB XY transpose for display

zimg = zeros(pandims);
for i = 1:nslices
  ixyzmm = [x(:)';y(:)';ones(1,nvox)*zmm(i);ones(1,nvox)];
  img = zimg;
  for j = 1:nimgs
    thisimg = SO.img(j);
    % to voxel space of image
    vixyz = inv(SO.transform*thisimg.vol.mat)*ixyzmm;
    % raw data 
    if is_there(thisimg.vol, 'imgdata')
      V = thisimg.vol.imgdata;
    else
      V = thisimg.vol;
    end
    i1 = spm_sample_vol(V,vixyz(X,:),vixyz(Y,:),vixyz(Z,:), ...
			 [thisimg.hold thisimg.background]);
    if is_there(thisimg, 'func')
      eval(thisimg.func);
    end
    % transpose to reverse X and Y for figure
    i1 = reshape(i1, vdims(1:2))';
    % rescale to colormap
    [csdata badvals]= scaletocmap(...
	i1,...
	thisimg.range(1),...
	thisimg.range(2),...
	cmaps{j},...
	lrn(j,:));
    % take indices from colormap to make true colour image
    iimg = reshape(cmaps{j}(csdata(:),:),pandims);
    tmp = repmat(logical(~badvals),[1 1 3]);
    if thisimg.prop ~= Inf % truecolor overlay
      img(tmp) = img(tmp) + iimg(tmp)*thisimg.prop;
    else % split colormap effect
      img(tmp) = iimg(tmp);
    end
  end
  % threshold out of range values
  img(img>1) = 1;
  
  image('Parent', axisd(i),...
	'ButtonDownFcn', SO.callback,...
	'CData',img);
  if is_there(SO,'labels')
    text('Parent',axisd(i),...
	 'Color', labels.colour,...
	 'FontUnits', 'normalized',...
	 'VerticalAlignment','bottom',...
	 'HorizontalAlignment','left',...
	 'Position', [1 1],...
	 'FontSize',labels.size,...
	 'ButtonDownFcn', SO.callback,...
	 'String', labels.format{i});
  end
end
for i = (nslices+1):npanels
   set(axisd(i),'Color',[0 0 0]);
end
% add colorbar(s) 
for i = 1:cbars
  axno = axisd(end-cbars+i);
  cbari = SO.img(SO.cbar(i));
  cml = size(cbari.cmap,1);
  p = get(axno, 'Position');; % position of last axis
  cw = p(3)*0.2;
  ch = p(4)*0.75;
  pc = p(3:4)/2;
  [axlims idxs] = sort(cbari.range);
  a=axes(...
      'Parent',figno,...
      'XTick',[],...
      'XTickLabel',[],...
      'Units', 'pixels',...
      'YLim', axlims,...   
      'FontUnits', 'normalized',...
      'FontSize', 0.075,...
      'YColor',[1 1 1],...
      'Tag', 'cbar',...
      'Box', 'off',...
      'Position',[p(1)+pc(1)-cw/2,p(2)+pc(2)-ch/2,cw,ch]...
      );
  ih = image('Parent', a,...
	'YData', axlims(idxs),...     
	'CData', reshape(cbari.cmap,[cml,1,3]));

end

 otherwise
  error(sprintf('Unrecognized action string %s', action));

% end switch action
end

return

function checkso
% checks and fills SO structure
global SO

% figure
if is_there(SO, 'figure')
  try
    if ~strcmp(get(SO.figure,'Type'),'figure')
      error('Figure handle is not a figure')
    end
  catch
    error('Figure handle is not a valid figure')
  end
else
  % no figure handle. Try spm figure, then gcf
  SO.figure = spm_figure('FindWin', 'Graphics'); 
  if isempty(SO.figure)
    SO.figure = gcf;
  end
end
% set defaults for SPM figure 
if strcmp(get(SO.figure, 'Tag'),'Graphics')
  % position figure nicely for SPM
  defstruct = struct('position', [0 0 1 0.92], 'units', 'normalized', ...
		     'valign', 'top');
  SO = set_def(SO, 'area', defstruct);
  SO.area = set_def(SO.area, 'position', defstruct.position); 
  SO.area = set_def(SO.area, 'units', defstruct.units); 
  SO.area = set_def(SO.area, 'valign', defstruct.valign); 
end
SO = set_def(SO, 'clf', 1);

% orientation; string or 4x4 matrix
orientn = [];
SO = set_def(SO, 'transform', 'axial');
if ischar(SO.transform)
  orientn = find(strcmpi(SO.transform, {'axial','coronal','sagittal'}));
  if isempty(orientn)
    error(sprintf('Unexpected orientation %s', SO.transform));
  end
  ts = [0 0 0 0 0 0 1 1 1;...
      0 0 0 pi/2 0 0 1 -1 1;...
      0 0 0 pi/2 0 -pi/2 -1 1 1];
  SO.transform = spm_matrix(ts(orientn,:));
end
% default slice size, slice matrix depends on orientation
if ~is_there(SO,'slicedef' | ~is_there(SO, 'slices'))
  % take image sizes from first image
  V = SO.img(1).vol;
  D = V.dim(1:3);
  T = SO.transform * V.mat;
  vcorners = [1 1 1; D(1) 1 1; 1 D(2) 1; D(1:2) 1; ...
	     1 1 D(3); D(1) 1 D(3); 1 D(2:3) ; D(1:3)]';
  corners = T * [vcorners; ones(1,8)];
  SC = sort(corners');
  vxsz = sqrt(sum(T(1:3,1:3).^2));
  
  SO = set_def(SO, 'slicedef',...
    [SC(1,1) vxsz(1) SC(8,1);SC(1,2) vxsz(2) SC(8,2)]);
  SO = set_def(SO, 'slices',[SC(1,3):vxsz(3):SC(8,3)]);
end

% no colourbars by default
SO = set_def(SO, 'cbars', []);

% always refresh figure window, by default
SO = set_def(SO, 'refreshf', 1);  

% labels
defstruct = struct('colour',[1 1 1],'size',0.075,'format', '%+3.0f');
if ~isfield(SO, 'labels') % no field, -> default
  SO.labels = defstruct;
elseif ~isempty(SO.labels) % empty -> no labels
  % colour for slice labels
  SO.labels = set_def(SO.labels, 'colour', defstruct.colour); 
  % font size normalized to image axis
  SO.labels = set_def(SO.labels, 'size', defstruct.size); 
  % format string for slice labels
  SO.labels = set_def(SO.labels, 'format', defstruct.format); 
end

% callback
SO = set_def(SO, 'callback', ';');

% figure area stuff
defarea = struct('position',[0 0 1 1],'units','normalized');
SO = set_def(SO, 'area', defarea);
if ~is_there(SO.area, 'position')
  SO.area = defarea;
end
if ~is_there(SO.area,'units')
  if (all(SO.area.position>=0 & SO.area.position<=1))
    SO.area.units = 'normalized';
  else
    SO.area.units = 'pixels';
  end
end
SO.area = set_def(SO.area,'halign', 'center');
SO.area = set_def(SO.area,'valign', 'middle');

% printing
SO = set_def(SO, 'printstr', 'print -dpsc -painters -noui');
SO = set_def(SO, 'printfile', 'slices.ps');

% fill various img arguments
% would be nice to use set_def, but we can't

% set colour intensities as we go
remcol = 1;
for i = 1:length(SO.img)
  if ~is_there(SO.img(i),'hold')
    if ~is_there(SO.img(i).vol,'imgdata')
      % normal file vol struct
      SO.img(i).hold = 1;
    else
      % 3d matrix vol struct
      SO.img(i).hold = 0;
    end
  end
  if ~is_there(SO.img(i),'background')
    SO.img(i).background = NaN;
  end
  if ~is_there(SO.img(i),'prop')
    % default is true colour
    SO.img(i).prop = remcol/(length(SO.img)-i+1);
    remcol = remcol - SO.img(i).prop;
  end
  if ~is_there(SO.img(i),'range')
    [mx mn] = volmaxmin(SO.img(i).vol);
    SO.img(i).range = [mn mx];
  end
  if ~is_there(SO.img(i),'cmap')
    if SO.img(i).prop == Inf; % split map
      if SO.range(1)<SO.range(2)
	SO.img(i).cmap = getcmap('hot');
      else
	SO.img(i).cmap = getcmap('winter');
      end
    else                  % true colour
      SO.img(i).cmap = getcmap('actc');
    end
  end  
  if ~is_there(SO.img(i),'outofrange')
    % this can be complex, and depends on split/true colour
    if SO.img(i).prop == Inf % split colour
      if xor(SO.img(i).range(1) < SO.img(i).range(2), ...
	     SO.img(i).range(2) < 0)
	SO.img(i).outofrange = {[0],size(SO.img(i).cmap,1)};
      else
	SO.img(imgno).outofrange={[1], [0]};
      end
    else            % true colour
      SO.img(i).outofrange = {1,size(SO.img(i).cmap,1)};
    end
  end
  for j=1:2
    if isempty(SO.img(i).outofrange{j})
      SO.img(i).outofrange(j) = {0};
    end
  end
  if ~is_there(SO.img(i),'nancol')
    SO.img(i).nancol = 0;
  end
end  
return

function tf = is_there(a, fname)
% returns true if field fname is present in struct a, and not empty
tf = isfield(a, fname);
if tf
  tf = ~isempty(getfield(a, fname));
end
return

% function [img, badvals]=scaletocmap(inpimg,mn,mx,cmap,lrn)
% img = (inpimg-mn)/(mx-mn);  % img normalized to mn=0,mx=1
% cml = size(cmap,1);
% if cml==1 % values between 0 and 1 -> 1
%   img(img>=0 & img<=1)=1;
% else
%   img = img*(cml-1)+1;
% end
% outvals = {img<1, img>cml, isnan(img)};
% img= round(img);
% badvals = zeros(size(img));
% for i = 1:length(lrn)
%   if lrn(i)
%     img(outvals{i}) = lrn(i);
%   else
%     badvals = badvals | outvals{i};
%     img(outvals{i}) = 1;
%   end    
% end
% return

function st = set_def(st, fld, def)
if ~is_there(st, fld)
  st = setfield(st, fld, def);
end
return

% function addblobs(imgno, xyz,vals,mat)
% global SO
% if isempty(imgno)
%   imgno = length(SO.img);
% end
% if ~isempty(xyz)
%   SO.img(imgno).vol = blobs2vol(xyz,vals,mat);
% end

function vol = blobs2vol(xyz,vals,mat)
vol = [];
if ~isempty(xyz),
  rcp      = round(xyz);
  vol.dim  = max(rcp,[],2)';
  off      = rcp(1,:) + vol.dim(1)*(rcp(2,:)-1+vol.dim(2)*(rcp(3,:)-1));
  vol.imgdata = zeros(vol.dim)+NaN;
  vol.imgdata(off) = vals;
  vol.imgdata      = reshape(vol.imgdata,vol.dim);
  vol.mat = mat;
end
return

function addmatrix(imgno,mat3d,mat)
global SO
if isempty(imgno)
  imgno = length(SO.img);
end
if nargin<3
  mat = [];
end
if ~isempty(mat3d)
  SO.img(imgno).vol = matrix2vol(mat3d,mat);
end

function vol = matrix2vol(mat3d,mat)
if nargin < 2
  mat = spm_matrix([]);
end
if isempty(mat)
  mat = spm_matrix([]);
end
vol = [];
if ~isempty(mat3d)
  vol.imgdata = mat3d;
  vol.mat = mat;
  vol.dim = size(mat3d);
end
return

function [mx,mn] = volmaxmin(vol)
if is_there(vol, 'imgdata')
  tmp = vol.imgdata(isfinite(vol.imgdata));
  mx = max(tmp);
  mn = min(tmp);
else
    mx = -Inf;mn=Inf;
    for i=1:vol.dim(3),
      tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),[0 NaN]);
      tmp = tmp(find(isfinite(tmp(:))));
      if ~isempty(tmp)
	mx = max([mx; tmp]);
	mn = min([mn; tmp]);
      end
    end
end
return

% function cmap = getcmap(acmapname)
% % get colormap of name acmapname
% if ~isempty(acmapname)
%   cmap = evalin('base',acmapname,'[]');
%   if isempty(cmap) % not a matrix, is it...
%     % a colour name?
%     tmp = strcmp(acmapname, {'red','green','blue'});
%     if any(tmp)
%       cmap = zeros(64,3);
%       cmap(:,tmp) = ((0:63)/63)';
%     else
%       % a .mat file?
%       [p f e] = fileparts(acmapname);
%       acmat = fullfile(p, [f '.mat']);
%       if exist(acmat, 'file')
% 	s = struct2cell(load(acmat));
% 	cmap = s{1};
%       end
%     end
%   end
% end
% if size(cmap, 2)~=3
%   warning('Colormap was not an N by 3 matrix')
%   cmap = [];
% end
% return

function mmpos = getpos
% returns point location from last click, in mm
global SO
mmpos=[];
pos = get(gca, 'CurrentPoint');
u = get(gca, 'UserData');
if is_there(u, 'type')
  if strcmp(u.type, 'slice') % is slice panel
    mmpos = (pos(1,1:2)'-1).*SO.slicedef(:,2)+SO.slicedef(:,1);
    mmpos = inv(SO.transform) * [mmpos; u.no; 1];
    mmpos = mmpos(1:3,1);
  end
end
return

function vals = pointvals(XYZmm, holdlist)
% returns values from all the images at points given in XYZmm
global SO
if nargin < 2
  holdlist = [SO.img(:).hold];
end
X=1;Y=2;Z=3;
nimgs = length(SO.img);
nvals = size(XYZmm,2);
vals = zeros(nimgs,nvals)+NaN;
if size(XYZmm,1)~=4
  XYZmm = [XYZmm(X:Z,:); ones(1,nvals)];
end
for i = 1:nimgs
  I = SO.img(i);
  XYZ = I.vol.mat\XYZmm;
  if ~is_there(I.vol, 'imgdata')
    vol = I.vol;
  else
    vol = I.vol.imgdata;
  end
  vals(i,:) = spm_sample_vol(vol, XYZ(X,:), XYZ(Y,:),XYZ(Z,:),[holdlist(i) ...
		    I.background]);
end  
return

function printfig(filename,printstr)
% print slice overlay figure
% based on spm_figure print, and including fix from thence for ps printing
global SO;
if nargin < 1
  filename = [];
end
if isempty(filename)
  filename = SO.printfile;
end
if nargin < 2
  printstr = '';
end
if isempty(printstr)
  printstr = SO.printstr;
end

%-Note current figure, & switch to figure to print
cF = get(0,'CurrentFigure');
set(0,'CurrentFigure',SO.figure)

%-Temporarily change all units to normalized prior to printing
% (Fixes bizzarre problem with stuff jumping around!)
%-----------------------------------------------------------------------
H  = findobj(get(SO.figure,'Children'),'flat','Type','axes');
un = cellstr(get(H,'Units'));
set(H,'Units','normalized')

%-Print
%-----------------------------------------------------------------------
err = 0;
try, eval([printstr ' ' filename]), catch, err=1; end
if err
	errstr = lasterr;
	tmp = [find(abs(errstr)==10),length(errstr)+1];
	str = {errstr(1:tmp(1)-1)};
	for i = 1:length(tmp)-1
		if tmp(i)+1 < tmp(i+1) 
			str = [str, {errstr(tmp(i)+1:tmp(i+1)-1)}];
		end
	end
	str = {str{:},	'','- print command is:',['    ',printstr ' ' filename],...
			'','- current directory is:',['    ',pwd],...
			'','            * nothing has been printed *'};
	for i=1:length(str)
	  disp(str{i});end
end

set(H,{'Units'},un)
set(0,'CurrentFigure',cF)

return