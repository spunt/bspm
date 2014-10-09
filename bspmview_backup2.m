function S = bspmview(ol, ul)
% BSPMVIEW Program for viewing fMRI statistical maps
%
%  USAGE: bspmview(ol*, ul*)	*optional input
%
%  Requires that Statistical Parametric Mapping (SPM; Wellcome Trust Centre
%  for Neuroimaging; www.fil.ion.ucl.ac.uk/spm/) be in your MATLAB search
%  path. It has only been tested on SPM8/SPM12b and may not function
%  correctly with earlier versions of SPM. 
% _________________________________________________________________________
%  INPUTS
%	ol: filename for statistical image to overlay
%	ul: filename for anatomical image to use as underlay
%
% _________________________________________________________________________
%  EXAMPLES
%   >> bspmview('spmT_0001.img', 'T1.nii')  
%   >> bspmview('spmT_0001.img')   % uses default underlay
%	>> bspmview                    % opens dialogue for selecting overlay
%   
% _________________________________________________________________________
%  CREDITS
%	This software heavily relies on functions contained within the SPM
%	software, and in essence a translation of their functionality into a
%	simpler and more user-friendly format. In addition, this software was
%	inspired by and in some cases uses code segments from two other
%	statistical image viewers: XJVIEW.m by Xu Cui, Jian Li, and Xiaowei
%	Song (http://www.alivelearn.net/xjview8/developers/), and FIVE.m by
%	Aaron P. Schultz (http://mrtools.mgh.harvard.edu/index.php/Main_Page).
%

% ---------------------- Copyright (C) 2014 Bob Spunt ----------------------
%	Created:  2014-09-27
%	Email:    spunt@caltech.edu
% _________________________________________________________________________

% | CHECK INPUTS
% | =======================================================================
if nargin < 1 
    ol = uigetvol('Select an Image File for Overlay', 0);
    if isempty(ol), disp('Must select an overlay!'); return; end
else
    if iscell(ol), ol = char(ol); end
end
if nargin < 2
    ul=fullfile(fileparts(which('spm.m')), 'canonical', 'single_subj_T1.nii'); 
else
    if iscell(ul), ul = char(ul); end
end

% | GUI FIGURE
% | =======================================================================
try
    fonts   = default_fonts; 
    pos     = default_positions; 
    color   = default_colors; 
    S.hFig    = figure(...
    'Units', 'pixels', ...
    'Position',pos.gui,...
    'Resize','off',...
    'Color',color.bg,...
    'ColorMap',gray(64),...
    'NumberTitle','off',...
    'DockControls','off',...
    'MenuBar','none',...
    'Name','bspmVIEW',...
    'CloseRequestFcn', @cb_closegui, ...
    'DefaultTextColor',color.fg,...
    'DefaultTextInterpreter','none',...
    'DefaultTextFontName','Arial',...
    'DefaultTextFontSize',12,...
    'DefaultAxesColor',color.border,...
    'DefaultAxesXColor',color.border,...
    'DefaultAxesYColor',color.border,...
    'DefaultAxesZColor',color.border,...
    'DefaultAxesFontName','Arial',...
    'DefaultPatchFaceColor',color.fg,...
    'DefaultPatchEdgeColor',color.fg,...
    'DefaultSurfaceEdgeColor',color.fg,...
    'DefaultLineColor',color.border,...
    'DefaultUicontrolFontName',fonts.name,...
    'DefaultUicontrolFontSize',fonts.sz3,...
    'DefaultUicontrolInterruptible','on',...
    'Visible','off',...
    'Toolbar','none');
    set(S.hFig, 'ResizeFcn', @cb_resizegui); 
    uicontrol('Parent', S.hFig, 'Units', 'Normal', 'Style', 'Text', ...
    'pos', [0 0 1 .001], 'backg', color.blues(8,:));
    uicontrol('Parent', S.hFig, 'Units', 'Normal', 'Style', 'Text', ...
    'pos', [0 .001 .001 1], 'backg', color.blues(10,:));
    uicontrol('Parent', S.hFig, 'Units', 'Normal', 'Style', 'Text', ...
    'pos', [.999 .001 .001 .999], 'backg', color.blues(10,:));
catch lasterr
    rethrow(lasterr); 
end

% | INITIALIZE SPM REGISTRY & ORTHVIEWS
% | =======================================================================
% try
% | REGISTRY OBJECT (HREG)
S.hReg = uipanel('Parent',S.hFig,'Units','Pixels','Position',pos.pane.axes,...
        'BorderType', 'none', 'BackgroundColor',color.bg);
set(S.hReg, 'units', 'norm');
global st prevsect
prevsect    = ul;
% | CREATE GLOBAL VARIABLE ST, PREVSECT
bspm_orthviews('Reset');
st          = struct( ...
            'fig',          S.hFig,...
            'figax',        S.hReg,...
            'n',            0,...
            'bb',           [],...
            'callback',     {';'}, ...
            'Space',        eye(4),...
            'centre',       [],...
            'xhairs',       1,...
            'plugins',      {''},...
            'hld',          1,...
            'mode',         [],...
            'color',        color,...
            'pos',          pos,...
            'fonts',        fonts,...
            'direct',       'both',...
            'snap',         []);
st.vols     = cell(24,1);
st.ol       = load_overlay(ol, .001, 5);
spm_XYZreg('InitReg',S.hReg,st.ol.M,st.ol.DIM,[0;0;0]); % initialize registry object
st.ho = bspm_orthviews('Image', ul, [.025 .025 .95 .95]);
bspm_orthviews('MaxBB');
bspm_orthviews('Register', S.hReg);
setposition_axes; 
setthresh(st.ol.C0(3,:));
setxhaircolor;
setcolormap; 
put_figmenu; 
put_axesxyz; 
put_axesmenu;
put_lowerpane;
put_upperpane;
if nargout==1, S.handles = gethandles; end
% catch lasterr
%     rethrow(lasterr)
% end
% =========================================================================
% *
% * SUBFUNCTIONS
% *
% =========================================================================

% | GUI DEFAULTS
% =========================================================================
function color  = default_colors 
    color.bg        = [20/255 23/255 24/255];
    color.fg        = [248/255 248/255 248/255];
    color.border    = [023/255 024/255 020/255]*2;
    color.xhair     = [0.7020    0.8039    0.8902];
    color.panel     = [.01 .22 .34];
    color.blues = brewermap(40, 'Blues'); 
function fonts  = default_fonts
    fonts.name      = 'Arial'; 
    fonts.sz1       = 24;
    fonts.sz2       = 20; 
    fonts.sz3       = 16; 
    fonts.sz4       = 12; 
function pos    = default_positions 
    %% GENERAL
    screensize      = get(0, 'ScreenSize');
    pos.ss          = screensize(3:4);
    pos.gui         = [pos.ss(1:2)*.5 pos.ss(2)*.55 pos.ss(2)*.5];
    pos.aspratio    = pos.gui(3)/pos.gui(4);
    guiss = [pos.gui(3:4) pos.gui(3:4)]; 
    h1 = 5; 
    h2 = pos.gui(4) - h1; 
    h = h2-h1;
    w1 = 5; 
    w2 = pos.gui(3) - w1; 
    w = w2 - w1; 
    h1up = round(h*.925);
    hup = h-h1up; 
    %% PANELS   
    pos.pane.upper  = [w1 h1up w hup]; 
    pos.pane.axes   = [w1 h1 w h1up-h1]; 
    

    pos.peakvalue      = [.050 .050 .250 .50];
    pos.xyz            = [.325 .050 .325 .50];
    pos.clustersize    = [.675 .050 .250 .50];
    
    pos.maxval         = [.050 .050 .250 .50];
    pos.cmap           = [.325 .050 .325 .50];

    pos.pos            = [.100 .050 .250 .85];
    pos.neg            = [.400 .050 .250 .85];
    pos.posneg         = [.700 .050 .250 .85];
    
    pos.k              = [.025 .050 .200 .55];
    pos.tval           = [.250 .050 .225 .55];
    pos.pval           = [.500 .050 .275 .55];
    pos.df             = [.800 .050 .175 .55];
    
    pos.pslider        = [.050 .025 .900 .150];
function prop   = default_properties(varargin)
global st
prop.darkbg     = {'visible','on', 'clip', 'off', 'backg', st.color.bg, 'foreg', st.color.fg};
prop.lightbg    = {'visible','on', 'clip', 'off', 'backg', st.color.fg, 'foreg', [0 0 0]};
if ~isempty(varargin), prop.darkbg = [varargin{:} prop.darkbg]; prop.lightbg = [varargin{:} prop.lightbg]; end
prop.panel      = [prop.darkbg {'bordertype', 'none', 'titlepos', 'centertop', 'fontw', 'bold'}]; 
prop.edit       = [prop.lightbg {'style', 'edit', 'horiz', 'center'}];
prop.text       = [prop.darkbg {'style', 'text', 'horiz', 'center'}]; 
prop.popup      = [prop.lightbg {'style', 'popup'}]; 
prop.slider     = [prop.darkbg {'style', 'slide', 'min', 1.0000e-20, 'max', 1, 'sliderstep', [1 5], 'value', st.ol.P}];
prop.push       = [prop.darkbg {'style', 'push', 'horiz', 'center'}]; 
prop.radio      = [prop.darkbg {'style', 'radio', 'horiz', 'center'}];
prop.toggle     = [prop.darkbg {'style', 'toggle'}]; 
prop.checkbox   = [prop.darkbg {'style', 'check'}]; 
prop.listbox    = [prop.darkbg {'style', 'list'}]; 

% | GUI COMPONENTS
% =========================================================================
function put_upperpane(varargin)
    global st
    cnamepos = [.01 .01 .98 .90]; 
    prop = default_properties('units', 'pixels', 'fontu', 'norm', 'fonts', .60); 
    panelh  = uipanel('parent',st.fig, prop.panel{:}, 'pos', st.pos.pane.upper, 'tag', 'upperpanel');
    prop = default_properties('units', 'norm', 'fontu', 'norm', 'fonts', .60); 
    uicontrol('parent', panelh, prop.text{:}, 'pos', cnamepos, 'tag', 'ContrastName', 'string', st.ol.descrip); 
function put_lowerpane(varargin)
global st

% | Default properties
prop = default_properties('units', 'norm', 'fontn', 'arial', 'fonts', 19);  

% | Positioning
[h,axpos] = gethandles_axes;
lowpos = axpos(1,:);
lowpos(1) = axpos(3, 1); 
lowpos(3) = 1 - lowpos(1) - axpos(1,2);

% | Create the total panel
panelh = uipanel('parent', st.figax, prop.panel{:}, 'pos',lowpos, 'tag', 'lowerpanel'); 

% | Create each subpanel 
tpanepos = [.025 .500 .950 .225];
dpanepos = [.025 .025 .950 .200];
ipanepos = [.025 .750 .950 .225];
cpanepos = [.025 .250 .950 .225];
S.infopane   = uipanel('parent', panelh, prop.panel{:},'pos',ipanepos, 'title', 'Current Voxel');
S.threshpane = uipanel('parent', panelh, prop.panel{:},'pos',tpanepos,'title','Thresholding Options'); 
S.directpane = uipanel('parent', panelh, prop.panel{:},'pos',dpanepos,'title','Effect to Display'); 
S.cmappane  = uipanel('parent', panelh, prop.panel{:},'pos',cpanepos,'title','Color Map'); 

% | Uicontrols for Thresh Panel
Tpos = {st.pos.k st.pos.tval st.pos.pval st.pos.df};
Tstr = {'Extent' 'Thresh' 'P-Value' 'DF'};
Tdefvalues = [st.ol.K st.ol.U st.ol.P st.ol.DF];
Tstrform = {'%d' '%2.3f' '%d' '%d'}; 
Ecallback   = {@cb_updateoverlay, @cb_updateoverlay, @cb_updateoverlay, @cb_updateoverlay};
prop = default_properties('parent', S.threshpane, 'units', 'norm','fonts', 17);
prop2 = default_properties('parent', S.threshpane,'units', 'norm', 'fonts', 17);
for i = 1:length(Tstr)
    txpos = Tpos{i};
    txpos(2) = sum(txpos([2 4]))+.05;
    txpos(4) = .60*(1-txpos(4));  
    S.tx(i) = uicontrol(prop2.text{:}, 'pos', txpos, 'string', Tstr{i});
    S.ed(i) = uicontrol(prop.edit{:}, 'pos', Tpos{i}, 'string', sprintf(Tstrform{i}, Tdefvalues(i)), 'Tag', Tstr{i}, 'callback', Ecallback{i}); 
end

% | Uicontrols for Effects to Display Panel
prop = default_properties('parent', S.directpane, 'units', 'norm','fonts', 17); 
S.pos1 = uicontrol(prop.radio{:}, 'pos', st.pos.pos, 'tag', 'direct', 'str', 'positive', 'callback', @cb_directmenu);
S.pos2 = uicontrol(prop.radio{:}, 'pos', st.pos.neg, 'tag', 'direct', 'str', 'negative', 'callback', @cb_directmenu); 
S.pos3 = uicontrol(prop.radio{:}, 'pos', st.pos.posneg, 'tag', 'direct', 'str', 'pos/neg', 'value', 1, 'enable', 'inactive', 'callback', @cb_directmenu); 

% | Uicontrols for Current Voxel Panel
prop = default_properties('parent', S.infopane, 'units', 'norm','fonts', 17);
prop2 = default_properties('parent', S.infopane, 'units', 'norm','fonts', 17);
S.peakvalue = uicontrol(prop.edit{:}, 'pos', st.pos.peakvalue, 'tag', 'voxval', 'enable', 'inactive');
txpos = st.pos.peakvalue; txpos(2) = sum(txpos([2 4]))+.05; txpos(4) = .60*(1-txpos(4));  
S.peakvaluetx = uicontrol(prop2.text{:}, 'pos',  txpos, 'str', 'Value');
S.xyz = uicontrol(prop.edit{:}, 'pos', st.pos.xyz, 'tag', 'xyz', 'callback', @cb_changexyz);
txpos = st.pos.xyz; txpos(2) = sum(txpos([2 4]))+.05; txpos(4) = .60*(1-txpos(4));  
S.xyztx = uicontrol(prop2.text{:}, 'pos', txpos, 'str', 'Coordinate');
S.xyz = uicontrol(prop.edit{:}, 'pos', st.pos.clustersize, 'tag', 'clustersize', 'enable', 'inactive');
txpos = st.pos.clustersize; txpos(2) = sum(txpos([2 4]))+.05; txpos(4) = .60*(1-txpos(4));  
S.xyztx = uicontrol(prop2.text{:}, 'pos', txpos, 'str', 'Cluster Size');
setvoxelinfo; 

prop = default_properties('parent', S.cmappane, 'units', 'norm','fonts', 17);
S.maxval = uicontrol(prop.edit{:}, 'pos', st.pos.maxval, 'str', sprintf('%2.3f',max(st.ol.Z)), 'tag', 'maxval', 'callback', @cb_maxval);
txpos = st.pos.maxval; txpos(2) = sum(txpos([2 4]))+.05; txpos(4) = .60*(1-txpos(4));  
S.peakvaluetx = uicontrol(prop.text{:}, 'pos',  txpos, 'str', 'Color Max');
function put_figmenu
    global st
    %% Main Menu
    S.menu1 = uimenu('Parent', st.fig, 'Label', 'bspmVIEW');
    S.opencode  = uimenu(S.menu1, 'Label','Open GUI M-File', 'Callback', @cb_opencode); 
    S.exit      = uimenu(S.menu1, 'Label', 'Exit', 'Callback', {@cb_closegui, st.fig});
    %% Load Menu
    S.load = uimenu(st.fig,'Label','Load');
    S.loadol = uimenu(S.load,'Label','Overlay Image','CallBack', @cb_loadol);
    S.loadul = uimenu(S.load,'Label','Underlay Image', 'Separator', 'on', 'CallBack', @cb_loadul);
    %% Save Menu
    S.save              = uimenu(st.fig,'Label','Save');
    S.saveintensity     = uimenu(S.save,'Label','Save as intensity image','CallBack', @cb_saveimg);
    S.savemask          = uimenu(S.save,'Label','Save as mask', 'Separator', 'on', 'CallBack', @cb_saveimg);
    %% Options Menu
    S.options   = uimenu(st.fig,'Label','Options');
    S.crosshair = uimenu(S.options,'Label','Show Crosshairs','Checked', 'on', 'CallBack', @cb_crosshair);
function put_axesmenu
    [h,axpos] = gethandles_axes;
    cmenu = uicontextmenu;
    ctxhair = uimenu(cmenu, 'Label', 'Toggle Crosshairs', 'checked', 'on', 'callback', @cb_crosshair); 
    ctmax   = uimenu(cmenu, 'Label', 'Go to global max', 'callback', @cb_minmax, 'separator', 'on');
    ctmin   = uimenu(cmenu, 'Label', 'Go to global min', 'callback', @cb_minmax);
    ctclustmax  = uimenu(cmenu, 'Label', 'Go to cluster max', 'callback', @cb_clustminmax);
    ctclustmin  = uimenu(cmenu, 'Label', 'Go to cluster min', 'callback', @cb_clustminmax);
    ctsavemap   = uimenu(cmenu, 'Label', 'Save as intensity image', 'callback', @cb_saveimg, 'separator', 'on');
    ctsavemask  = uimenu(cmenu, 'Label', 'Save as mask image', 'callback', @cb_saveimg);
    ctsavemap   = uimenu(cmenu, 'Label', 'Save current cluster', 'callback', @cb_saveclust, 'separator', 'on');
    ctsavemask  = uimenu(cmenu, 'Label', 'Save current cluster as mask', 'callback', @cb_saveclust);
    ctsavergb  = uimenu(cmenu, 'Label', 'Do screencapture', 'callback', @cb_savergb, 'separator', 'on');
    for a = 1:3
        set(h.ax(a), 'uicontextmenu', cmenu); 
    end
function put_axesxyz
    global st
    h = gethandles_axes;
    xyz = round(bspm_XYZreg('GetCoords',st.registry.hReg));
    xyzstr = num2str([-99; xyz]); 
    xyzstr(1,:) = [];
    set(h.ax, 'YAxislocation', 'right'); 
    axidx = [3 2 1]; 
    for a = 1:length(axidx)
        yh = get(h.ax(axidx(a)), 'YLabel');
        st.vols{1}.ax{axidx(a)}.xyz = yh;
        if a==1
            set(yh, 'units', 'norm', 'fontunits', 'norm', 'fontsize', .075, ...
                'pos', [0 1 0], 'horiz', 'left', 'fontname', 'arial', ...
                'color', [1 1 1], 'string', xyzstr(a,:), 'rot', 0);
            set(yh, 'fontunits', 'points'); 
            fs = get(yh, 'fontsize');
        else
            set(yh, 'units', 'norm', 'fontsize', fs, ...
                'pos', [0 1 0], 'horiz', 'left', 'fontname', 'arial', ...
                'color', [1 1 1], 'string', xyzstr(a,:), 'rot', 0);
        end
    end

% | GUI CALLBACKS
% =========================================================================
function cb_updateoverlay(varargin)
    global st
    T0 = getthresh;
    T = T0; 
    tag = get(varargin{1}, 'tag');
    di = strcmpi({'positive' 'negative' 'both'}, T.direct); 
    switch tag
        case {'Thresh'}
            T.pval = bob_t2p(T.thresh, T.df);
        case {'P-Value'}
            T.thresh = bob_p2t(T.pval, T.df); 
        case {'DF'}
            T.thresh = bob_p2t(T.pval, T.df); 
            T.pval = bob_t2p(T.thresh, T.df);
        case {'Extent'}
            if sum(st.ol.C0(di,st.ol.C0(di,:)>=T.extent))==0
                headsup('No clusters survived. Defaulting to largest cluster at this voxelwise threshold.');
                T.extent = max(st.ol.C0(di,:));
            end
    end
    [st.ol.C0, st.ol.C0IDX] = getclustidx(st.ol.Y, T.thresh, T.extent);
    C = st.ol.C0(di,:); 
    if sum(C(C>=T.extent))==0
        T0.thresh = st.ol.U; 
        setthreshinfo(T0);
        headsup('No voxels survived. Try a different threshold.'); 
        return
    end
    setthresh(C, find(di)); 
    setthreshinfo(T);  
function cb_loadol(varargin)
    global st
    fname = uigetvol('Select an Image File for Overlay', 0);
    if isempty(fname), disp('Must select an overlay!'); return; end
    T = getthresh; 
    st.ol = load_overlay(fname, T.pval, T.extent);
    setthresh(st.ol.C0(3,:));
    setcolormap; 
    setposition_axes;
    setcontrastname; 
function cb_loadul(varargin)
    global st
    ul = uigetvol('Select an Image File for Underlay', 0);
    bspm_orthviews('Delete', st.ho); 
    st.ho = bspm_orthviews('Image', ul, [.025 .025 .95 .95]);
    bspm_orthviews('MaxBB');
    bspm_orthviews('AddBlobs', st.ho, st.ol.XYZ, st.ol.Z, st.ol.M);
    bspm_orthviews('Register', st.registry.hReg);
    setposition_axes; 
    put_axesxyz;
    put_axesmenu;
function cb_clustminmax(varargin)
global st
str = get(findobj(st.fig, 'tag', 'clustersize'), 'string'); 
if strcmp(str, 'n/a'), return; end
[xyz, voxidx] = getnearestvoxel;
clidx = spm_clusters(st.ol.XYZ);
clidx = clidx==(clidx(voxidx)); 
tmpXYZmm = st.ol.XYZmm(:,clidx); 
if regexp(get(varargin{1}, 'label'), 'cluster max')
    centre = tmpXYZmm(:,st.ol.Z(clidx)==max(st.ol.Z(clidx)));
elseif regexp(get(varargin{1}, 'label'), 'cluster min')
    centre = tmpXYZmm(:,st.ol.Z(clidx)==min(st.ol.Z(clidx)));
end
bspm_orthviews('reposition', centre);
function cb_minmax(varargin)
global st
lab = get(varargin{1}, 'label');
if regexp(lab, 'global max')
    centre = st.ol.XYZmm(:,st.ol.Z==max(st.ol.Z));
elseif regexp(lab, 'global min')
    centre = st.ol.XYZmm(:,st.ol.Z==min(st.ol.Z)); 
end
bspm_orthviews('reposition', centre); 
function cb_maxval(varargin)
val = str2num(get(varargin{1}, 'string')); 
bspm_orthviews('SetBlobsMax', 1, 1, val)
function cb_changexyz(varargin)
xyz = str2num(get(varargin{1}, 'string')); 
bspm_orthviews('reposition', xyz'); 
function cb_directmenu(varargin)
    global st
    str     = get(varargin{1}, 'string');
    allh = findobj(st.fig, 'Tag', 'direct'); 
    allhstr = get(allh, 'String');
    set(allh(strcmp(allhstr, str)), 'Value', 1, 'Enable', 'inactive'); 
    set(allh(~strcmp(allhstr, str)), 'Value', 0, 'Enable', 'on');
    T = getthresh;
    di = strcmpi({'positive' 'negative' 'both'}, T.direct);
    [st.ol.C0, st.ol.C0IDX] = getclustidx(st.ol.Y, T.thresh, T.extent);
    C = st.ol.C0(di,:);
    if sum(C>0)==0 
        headsup('Nothing survives at this threshold. Showing unthresholded image.');
        T.thresh = 0; 
        T.pval = bob_t2p(T.thresh, T.df);
        T.extent = 1; 
        [st.ol.C0, st.ol.C0IDX] = getclustidx(st.ol.Y, T.thresh, T.extent);
        C = st.ol.C0(di,:);
        setthreshinfo(T); 
    end
    setthresh(C, find(di));    
function cb_opencode(varargin)
    open(mfilename('fullpath'));
function cb_closegui(varargin)
    if length(varargin)==3
        delete(varargin{3}); 
    else
        delete(varargin{1}); 
    end
function cb_crosshair(varargin)
    state = get(varargin{1},'Checked');
    if strcmpi(state,'on');
        spm_orthviews('Xhairs','off')
        set(varargin{1},'Checked','off');
    end
    if strcmpi(state,'off');
        spm_orthviews('Xhairs','on')
        set(varargin{1},'Checked','on');
    end
function cb_resizegui(varargin)
    pos     = default_positions;
    cpos    = get(varargin{1}, 'pos');
    set(varargin{1}, 'pos', [cpos(1:2) cpos(4)*pos.aspratio cpos(4)]);
    put_axesxyz;  
function cb_saveimg(varargin)
    global st
    lab = get(varargin{1}, 'label');
    T = getthresh; 
    di = strcmpi({'positive' 'negative' 'both'}, T.direct); 
    clustidx = st.ol.C0(di,:);
    opt = [1 -1 1]; 
    outimg = st.ol.Y*opt(di);
    outhdr = st.ol.hdr; 
    outimg(clustidx==0) = NaN;
    putmsg = 'Save intensity image as'; 
    outhdr.descrip = 'Thresholded Intensity Image'; 
    [p,n] = fileparts(outhdr.fname); 
    deffn = sprintf('%s/Thresh_%s.nii', p, n);  
    if regexp(lab, 'Save as mask')
        outimg(outimg>0) = 1; 
        outhdr.descrip = 'Thresholded Mask Image'; 
        putmsg = 'Save mask image as'; 
        deffn = sprintf('%s/Mask_%s.nii', p, n);  
    end
    fn = uiputvol(deffn, putmsg);
    if isempty(fn), disp('User cancelled.'); return; end
    outhdr.fname = fn; 
    spm_write_vol(outhdr, outimg);
    fprintf('\nImage saved to %s\n', fn);         
function cb_saveclust(varargin)
    global st
    str = get(findobj(st.fig, 'tag', 'clustersize'), 'string'); 
    if strcmp(str, 'n/a'), return; end
    [xyz, voxidx] = bspm_XYZreg('NearestXYZ', round(st.centre), st.ol.XYZmm0);
    lab = get(varargin{1}, 'label');
    T = getthresh; 
    di = strcmpi({'positive' 'negative' 'both'}, T.direct);
    clidx = st.ol.C0IDX(di,:);
    clidx = clidx==(clidx(voxidx)); 
    opt = [1 -1 1]; 
    outimg = st.ol.Y*opt(di);
    outhdr = st.ol.hdr; 
    outimg(clidx==0) = NaN;
    putmsg = 'Save current cluster'; 
    outhdr.descrip = 'Thresholded Cluster Image'; 
    [p,n] = fileparts(outhdr.fname); 
    deffn = sprintf('%s/Cluster_%s_x=%d_y=%d_z=%d_%svoxels.nii', p, n, xyz, str);  
    if regexp(lab, 'Save current cluster as mask')
        outimg(outimg>0) = 1; 
        outhdr.descrip = 'Thresholded Mask Image'; 
        putmsg = 'Save mask image as'; 
        deffn = sprintf('%s/ClusterMask_%s_x=%d_y=%d_z=%d_%svoxels.nii', p, n, xyz, str);   
    end
    fn = uiputvol(deffn, putmsg);
    if isempty(fn), disp('User cancelled.'); return; end
    outhdr.fname = fn; 
    spm_write_vol(outhdr, outimg);
    fprintf('\nCluster image saved to %s\n', fn);     
function cb_savergb(varargin)
    %% Handles for axes
    % 1 - transverse
    % 2 - coronal
    % 3 - sagittal 
    % st.vols{1}.ax{1}.ax   - axes
    % st.vols{1}.ax{1}.d    - image
    % st.vols{1}.ax{1}.lx   - crosshair (x)
    % st.vols{1}.ax{1}.ly   - crosshair (y)
    global st
    setbackgcolor;
    im = screencapture(st.fig);
    setbackgcolor(st.color.bg)
    [imname, pname] = uiputfile({'*.png; *.jpg; *.pdf', 'Image'; '*.*', 'All Files (*.*)'}, 'Specify output directory and name', construct_filename);
    if isempty(imname), disp('User cancelled.'); return; end
    imwrite(im, fullfile(pname, imname)); 
    fprintf('\nImage saved to %s\n', fullfile(pname, imname));   
    
% | SETTERS
% =========================================================================
function setcontrastname
global st
connamh = findobj(st.fig, 'Tag', 'ContrastName'); 
set(connamh, 'String', st.ol.descrip); 
function setposition_axes
    global st
    %% Handles for axes
    % 1 - transverse
    % 2 - coronal
    % 3 - sagittal 
    % st.vols{1}.ax{1}.ax   - axes
    % st.vols{1}.ax{1}.d    - image
    % st.vols{1}.ax{1}.lx   - crosshair (x)
    % st.vols{1}.ax{1}.ly   - crosshair (y)
    h = gethandles_axes;
    axpos = cell2mat(get(h.ax, 'pos'));
    CBPIXSIZE = 80; 
    axpos(1:2, 1)   = 0; 
    axpos(1, 2)     = 0;
    axpos(3, 1)     = sum(axpos(2,[1 3]))+.005; 
    axpos(2:3, 2)   = sum(axpos(1,[2 4]))+.005;
    pz  = axpos(1,:);
    py  = axpos(2,:);
    px  = axpos(3,:);
    zrat = pz(3)/pz(4);
    yrat = py(3)/py(4);
    xrat = px(3)/px(4);
    VL = sum(py([2 4])); 
    while VL < 1
        px(4) = px(4) + .001; 
        px(3) = px(4)*xrat; 
        py(4) = px(4); 
        py(3) = py(4)*yrat; 
        pz(3) = py(3); 
        pz(4) = pz(3)/zrat; 
        px(1) = sum(py([1 3]))+.005;
        py(2) = sum(pz([2 4]))+.005; 
        VL = sum(py([2 4]));
    end
    axpos = [pz; py; px]; 
    for a = 1:3, set(h.ax(a), 'position', axpos(a,:)); end
    set(h.ax, 'units', 'pixels'); 
    axpos = cell2mat(get(h.ax, 'pos'));
    HL = round(sum(axpos(3, [1 3])) + CBPIXSIZE); 
    figsize = get(st.fig, 'pos'); 
    figsize(3) = HL; 
    set(st.fig, 'pos', figsize)
    for a = 1:3, set(h.ax(a), 'position', axpos(a,:)); end
    set(h.ax, 'units', 'norm');
    bspm_orthviews('Redraw');
function setthreshinfo(T)
    global st
    Tval = [T.extent T.thresh T.pval T.df]; 
    Tstr = {'Extent' 'Thresh' 'P-Value' 'DF'};
    Tstrform = {'%d' '%2.3f' '%d' '%d'}; 
    for i = 1:length(Tstr)
        set(findobj(st.fig, 'Tag', Tstr{i}), 'String', sprintf(Tstrform{i}, Tval(i)));
    end
function setthresh(C, di)
    global st
    if nargin==1, di = 3; end
    idx = find(C > 0); 
    st.ol.XYZ   = st.ol.XYZ0(:,idx);
    st.ol.XYZmm = st.ol.XYZmm0(:,idx);
    st.ol.C     = C(idx); 
    st.ol.Z     = st.ol.Y(idx); 
    if di~=3, st.ol.Z = abs(st.ol.Y(idx)); end
    bspm_orthviews('RemoveBlobs', st.ho);
    bspm_orthviews('AddBlobs', st.ho, st.ol.XYZ, st.ol.Z, st.ol.M);
    bspm_orthviews('Register', st.registry.hReg);
    setposition_axes;
    setcontrastname;
    if di==3, setcolormap(jet(64)); else setcolormap(hot(64)); end
    bspm_orthviews('Reposition');
function [voxval, clsize] = setvoxelinfo
    global st
    [nxyz,voxidx,d] = bspm_XYZreg('NearestXYZ', round(st.centre), st.ol.XYZmm); 
    if d > min(st.ol.VOX)
        voxval = 'n/a'; 
        clsize = 'n/a';
    else
        voxval = sprintf('%2.3f', st.ol.Z(voxidx));
        clsize = sprintf('%d', st.ol.C(voxidx));
    end
    set(findobj(st.fig, 'tag', 'xyz'), 'string', sprintf('%d, %d, %d', round(st.centre)));
    set(findobj(st.fig, 'tag', 'voxval'), 'string', voxval); 
    set(findobj(st.fig, 'tag', 'clustersize'), 'string', clsize); 
function setbackgcolor(newcolor)
global st
if nargin==0, newcolor = [0 0 0]; end
prop = {'backg' 'ycolor' 'xcolor' 'zcolor'}; 
for i = 1:length(prop)
    set(findobj(st.fig, prop{i}, st.color.bg), prop{i}, newcolor); 
end
h = gethandles_axes;
set(h.ax, 'ycolor', newcolor, 'xcolor', newcolor); 
function setxhaircolor(varargin)
    global st
    h = gethandles_axes;
    set(h.lx, 'color', st.color.xhair); 
    set(h.ly, 'color', st.color.xhair);

% | GETTERS
% =========================================================================
function h = gethandles(varargin)
    global st
    h.axial = st.vols{1}.ax{1}.ax;
    h.coronal = st.vols{1}.ax{2}.ax;
    h.sagittal = st.vols{1}.ax{3}.ax;
    h.colorbar = st.vols{1}.blobs{1}.cbar;    
function [clustsize, clustidx] = getclustidx(rawol, u, k)

    % raw data to XYZ
    DIM         = size(rawol); 
    [X,Y,Z]     = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
    XYZ         = [X(:)';Y(:)';Z(:)'];
    pos  = zeros(1, size(XYZ, 2)); 
    neg  = pos; 
    clustidx = zeros(3, size(XYZ, 2));
    
    % positive
    supra = (rawol(:)>=u)';    
    if sum(supra)
        tmp      = spm_clusters(XYZ(:, supra));
        clbin      = repmat(1:max(tmp), length(tmp), 1)==repmat(tmp', 1, max(tmp));
        pos(supra) = sum(repmat(sum(clbin), size(tmp, 2), 1) .* clbin, 2)';
        clustidx(1,supra) = tmp;
    end
    pos(pos < k)    = 0; 
    
    % negative
    rawol = rawol*-1; 
    supra = (rawol(:)>=u)';    
    if sum(supra)
        tmp      = spm_clusters(XYZ(:, supra));
        clbin      = repmat(1:max(tmp), length(tmp), 1)==repmat(tmp', 1, max(tmp));
        neg(supra) = sum(repmat(sum(clbin), size(tmp, 2), 1) .* clbin, 2)';
        clustidx(2,supra) = tmp;
    end
    neg(neg < k) = 0;
     
    % both
    clustsize = [pos; neg]; 
    clustsize(3,:) = sum(clustsize);
    clustidx(3,:) = sum(clustidx); 
function [h, axpos] = gethandles_axes(varargin)
    global st
    axpos = zeros(3,4);
    if isfield(st.vols{1}, 'blobs');
        h.cb = st.vols{1}.blobs{1}.cbar; 
    end
    for a = 1:3
        tmp = st.vols{1}.ax{a};
        h.ax(a) = tmp.ax; 
        h.d(a)  = tmp.d;
        h.lx(a) = tmp.lx; 
        h.ly(a) = tmp.ly;
        axpos(a,:) = get(h.ax(a), 'position');
    end
function T = getthresh
    global st
    T.extent = str2num(get(findobj(st.fig, 'Tag', 'Extent'), 'String')); 
    T.thresh = str2num(get(findobj(st.fig, 'Tag', 'Thresh'), 'String'));
    T.pval = str2num(get(findobj(st.fig, 'Tag', 'P-Value'), 'String'));
    T.df = str2num(get(findobj(st.fig, 'Tag', 'DF'), 'String'));
    tmph = findobj(st.fig, 'Tag', 'direct'); 
    opt = get(tmph, 'String');
    T.direct = opt(find(cell2mat(get(tmph, 'Value'))));
    if strcmp(T.direct, 'pos/neg'), T.direct = 'both'; end
function [xyz, voxidx, dist] = getnearestvoxel 
    global st
    [xyz, voxidx, dist] = bspm_XYZreg('NearestXYZ', round(st.centre), st.ol.XYZmm);
    
% | COLORBAR STUFF
% =========================================================================
function setcolormap(newmap, interval)
    global st
    if nargin < 1, newmap = jet(64); end
    if nargin < 2, interval = [st.vols{1}.blobs{1}.min st.vols{1}.blobs{1}.max]; end
    cbh = st.vols{1}.blobs{1}.cbar; 
    cmap = [gray(64); newmap];
    set(findobj(cbh, 'type', 'image'), 'CData', (65:128)', 'CdataMapping', 'direct');
    set(st.fig,'Colormap', cmap);
    bspm_orthviews('SetBlobsMax', 1, 1, max(st.ol.Z))
    set(findobj(st.fig, 'tag', 'maxval'), 'str',  sprintf('%2.3f',max(st.ol.Z))); 
function addcolourbar(vh,bh)
    global st
    axpos = zeros(3, 4);
    for a = 1:3
        axpos(a,:) = get(st.vols{vh}.ax{a}.ax, 'position');
    end
    cbpos = axpos(3,:); 
    cbpos(4) = cbpos(4)*.9; 
    cbpos(2) = cbpos(2) + (axpos(3,4)-cbpos(4))/2; 
    cbpos(1) = sum(cbpos([1 3])); 
    cbpos(3) = (1 - cbpos(1))/2; 
    cbpos(3) = min([cbpos(3) .30]); 
    cbpos(1) = cbpos(1) + (cbpos(3)/4); 
    yl = [st.vols{vh}.blobs{bh}.min st.vols{vh}.blobs{bh}.max];
    if range(yl) < 1
        yltick = [min(yl) max(yl)]; 
    else
        yltick = [ceil(min(yl)) floor(max(yl))];
    end
    st.vols{vh}.blobs{bh}.cbar = axes('Parent', st.figax, 'ycolor', st.color.fg, ...
        'position', cbpos, 'YAxisLocation', 'right', 'fontsize', 12, ...
        'ytick', yltick, 'tag', 'colorbar', ...
        'Box','off', 'YDir','normal', 'XTickLabel',[], 'XTick',[]); 
    set(st.vols{vh}.blobs{bh}.cbar, 'fontweight', 'bold', 'fontsize', st.fonts.sz3, 'fontname', st.fonts.name); 
    if isfield(st.vols{vh}.blobs{bh},'name')
        ylabel(st.vols{vh}.blobs{bh}.name,'parent',st.vols{vh}.blobs{bh}.cbar);
    end    
function cmap = getcmap(acmapname)
% get colormap of name acmapname
if ~isempty(acmapname)
    cmap = evalin('base',acmapname,'[]');
    if isempty(cmap) % not a matrix, is .mat file?
        acmat = spm_file(acmapname, 'ext','.mat');
        if exist(acmat, 'file')
            s    = struct2cell(load(acmat));
            cmap = s{1};
        end
    end
end
if size(cmap, 2)~=3
    warning('Colormap was not an N by 3 matrix')
    cmap = [];
end
function redraw_colourbar(vh,bh,interval,cdata)
global st
axpos = zeros(3, 4);
for a = 1:3
    axpos(a,:) = get(st.vols{vh}.ax{a}.ax, 'position');
end
cbpos = axpos(3,:); 
cbpos(4) = cbpos(4)*.9; 
cbpos(2) = cbpos(2) + (axpos(3,4)-cbpos(4))/2; 
cbpos(1) = sum(cbpos([1 3])); 
cbpos(3) = (1 - cbpos(1))/2; 
cbpos(1) = cbpos(1) + (cbpos(3)/4);
% only scale cdata if we have out-of-range truecolour values
if ndims(cdata)==3 && max(cdata(:))>1
    cdata=cdata./max(cdata(:));
end
% yl = [st.vols{vh}.blobs{bh}.min st.vols{vh}.blobs{bh}.max]; 
yl = interval;
if range(yl) < 1
    yltick = [min(yl) max(yl)]; 
else
    yltick = [ceil(min(yl)) floor(max(yl))];
end
image([0 1],interval,cdata,'Parent',st.vols{vh}.blobs{bh}.cbar);
set(st.vols{vh}.blobs{bh}.cbar, 'ycolor', st.color.fg, ...
    'position', cbpos, 'YAxisLocation', 'right', ...
    'ytick', yltick, ...
    'Box','off', 'YDir','normal', 'XTickLabel',[], 'XTick',[]); 
set(st.vols{vh}.blobs{bh}.cbar, 'fontweight', 'bold', 'fontsize', st.fonts.sz3, 'fontname', st.fonts.name); 
if isfield(st.vols{vh}.blobs{bh},'name')
    ylabel(st.vols{vh}.blobs{bh}.name,'parent',st.vols{vh}.blobs{bh}.cbar);
end
    
% | IMAGE MANIPULATION UTILITIES
% =========================================================================
function OL = load_overlay(fname, pval, k)
    global st
    if nargin<3, k = 5; end
    if nargin<2, pval = .001; end
    oh = spm_vol(fname); 
    od = spm_read_vols(oh);
    od(isnan(od)) = 0; 
    %% DEGREES OF FREEDOM
    try
        tmp = oh.descrip;
        idx1 = regexp(tmp,'[','ONCE');
        idx2 = regexp(tmp,']','ONCE');
        df = str2num(tmp(idx1+1:idx2-1));
        u = bob_p2t(pval, df);  
    catch
        headsup('Degrees of freedom could not be found in image header! Showing unthresholded image.')
        u = 0.01;
        k = 1;
        df = Inf;
        pval = Inf; 
    end
    if isempty(u)
        u = 0.01;
        k = 1; 
        df = Inf; 
        pval = Inf; 
    end
    [C, I] = getclustidx(od, u, k);
    if ~any(C(:))
        headsup('No suprathreshold voxels! Showing unthresholded image.')
        u = 0; 
        pval = bob_t2p(u, df);
        k = 1; 
        [C, I] = getclustidx(od, u, k); 
    end
    M           = oh.mat;         %-voxels to mm matrix
    DIM         = oh.dim';
    VOX         = abs(diag(M(:,1:3))); 
    [X,Y,Z]     = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
    XYZ        = [X(:)';Y(:)';Z(:)'];
    RCP         = XYZ; 
    RCP(4,:)    = 1;
    XYZmm      = M(1:3,:)*RCP;
    OL          = struct( ...
                'fname',    fname,...
                'descrip',  oh.descrip, ...
                'hdr',      oh, ...
                'DF',       df, ...
                'U',        u, ...
                'P',        pval, ...
                'K',        k, ...
                'Y',        od, ...
                'M',        M,...
                'DIM',      DIM,...
                'VOX',      VOX, ...
                'C0',        C, ...
                'C0IDX',       I, ...
                'XYZmm0',    XYZmm,...
                'XYZ0',      XYZ);   
function t = bob_p2t(alpha, df)
% BOB_P2T Get t-value from p-value + df
%
%   USAGE: t = bob_p2t(alpha, df)
%       
%   OUTPUT
%       t = crtical t-value
%
%   ARGUMENTS
%       alpha = p-value
%       df = degrees of freedom
%
% =========================================
if nargin<2, disp('USAGE: bob_p2t(p, df)'); return, end
t = tinv(1-alpha, df);
function p = bob_t2p(t, df)
% BOB_T2P Get p-value from t-value + df
%
%   USAGE: p = bob_t2p(t, df)
%       
%   OUTPUT
%       p = p-value
%
%   ARGUMENTS
%       t = t-value
%       df = degrees of freedom
%
% =========================================
if nargin<2, disp('USAGE: bob_t2p(p, df)'); return, end
p = tcdf(t, df);
p = 1 - p;
function [extent, info] = cluster_correct(im,u,alpha,range)
% BOB_SPM_CLUSTER_CORRECT Computer extent for cluster-level correction
%
% USAGE: [k info] = bob_spm_cluster_correct(im,u,alpha,range)
%
%
% THIS IS A MODIFICATION OF A FUNCTION BY DRS. THOMAS NICHOLS AND MARKO
% WILKE, CorrClusTh.m. ORIGINAL DOCUMENTATION PASTED BELOW:
%
% Find the corrected cluster size threshold for a given alpha
% function [k,Pc] =CorrClusTh(SPM,u,alpha,guess)
% SPM   - SPM data structure
% u     - Cluster defining threshold
%         If less than zero, u is taken to be uncorrected P-value
% alpha - FWE-corrected level (defaults to 0.05)
% guess - Set to NaN to use a Newton-Rhapson search (default)
%         Or provide a explicit list (e.g. 1:1000) of cluster sizes to
%         search over.
%         If guess is a (non-NaN) scalar nothing happens, except the the
%         corrected P-value of guess is printed. 
%
% Finds the corrected cluster size (spatial extent) threshold for a given
% cluster defining threshold u and FWE-corrected level alpha. 
%
%_________________________________________________________________________
% $Id: CorrClusTh.m,v 1.12 2008/06/10 19:03:13 nichols Exp $ Thomas Nichols, Marko Wilke
if nargin < 1
    disp('USAGE: [k info] = bob_spm_cluster_correct(im,u,alpha,range)'); 
    return; 
end
if nargin < 2, u = .001; end
if nargin < 3, alpha = .05; end
if nargin < 4, range = 5:200; end
if iscell(im), im = char(im); end

%% Get Design Variable %%
[impath, imname] = fileparts(im);
if exist([impath filesep 'I.mat'],'file') 
    matfile = [impath filesep 'I.mat']; 
    maskfile = [impath filesep 'mask.nii'];
elseif exist([impath filesep 'SPM.mat'],'file') 
    matfile = [impath filesep 'SPM.mat'];
else
    disp('Could not find an SPM.mat or I.mat variable, exiting.'); extent = []; info = []; return
end

%% Defaults %%
epsP = 1e-6;   % Corrected P-value convergence criterion (fraction of alpha)
du   = 1e-6;   % Step-size for Newton-Rhapson
maxi = 100;    % Maximum interations for refined search
STAT = 'T';    % Test statistic

%% Determime SPM or GLMFLEX %%
if strfind(matfile,'SPM.mat'), flexflag = 0; else flexflag = 1; end

%% Load and Compute Params %%
if flexflag % GLMFLEX
    II = load(matfile);
    try
        mask.hdr = spm_vol([II.I.OutputDir filesep 'mask.nii']);
    catch
        [p mf] = fileparts(im);
        mask.hdr = spm_vol([p filesep 'mask.nii']);
    end
    mask.data = spm_read_vols(mask.hdr);
    img.hdr = spm_vol(im);
    img.data = spm_read_vols(img.hdr);
    tmp = img.hdr.descrip; i1 = find(tmp=='['); i2 = find(tmp==']');
    df = str2num(tmp(i1(1)+1:i2(1)-1));
    df = [1 df];    
    n = 1;
    FWHM = II.I.FWHM{1};
    R = spm_resels_vol(mask.hdr,FWHM)';
    SS = sum(mask.data(:)==1);
    M = II.I.v.mat;
    VOX  = sqrt(diag(M(1:3,1:3)'*M(1:3,1:3)))';
    FWHMmm= FWHM.*VOX; % FWHM {mm}
    v2r  = 1/prod(FWHM(~isinf(FWHM)));% voxels to resels

else % SPM
    
    SPM = load(matfile);
    SPM = SPM.SPM;
    df   = [1 SPM.xX.erdf];
    STAT = 'T';
    n    = 1;
    R    = SPM.xVol.R;
    SS    = SPM.xVol.S;
    M    = SPM.xVol.M;
    VOX  = sqrt(diag(M(1:3,1:3)'*M(1:3,1:3)))';
    FWHM = SPM.xVol.FWHM;
    FWHMmm= FWHM.*VOX; 				% FWHM {mm}
    v2r  = 1/prod(FWHM(~isinf(FWHM))); %-voxels to resels
    
end
if ~nargout
    sf_ShowVolInfo(R,SS,VOX,FWHM,FWHMmm)
end
epsP = alpha*epsP;
Status = 'OK';
if u <= 1; u = spm_u(u,df,STAT); end

if length(range)==1 & ~isnan(range)
  
  %
  % Dummy case... just report P-value
  %

  k  = range;
  Pc = spm_P(1,k*v2r,u,df,STAT,R,n,SS);
  
  Status = 'JustPvalue';

elseif (spm_P(1,1*v2r,u,df,STAT,R,n,SS)<alpha)

  %
  % Crazy setting, where 1 voxel cluster is significant
  %

  k = 1;
  Pc = spm_P(1,1*v2r,u,df,STAT,R,n,SS);
  Status = 'TooRough';

elseif isnan(range)

  %
  % Automated search
  % 

  % Initial (lower bound) guess is the expected number of voxels per cluster
  [P Pn Em En EN] = spm_P(1,0,u,df,STAT,R,n,SS);
  kr = En; % Working in resel units
  rad = (kr)^(1/3); % Parameterize proportional to cluster diameter

  %
  % Crude linear search bound answer
  %
  Pcl  = 1;   % Lower bound on P
  radu = rad; % Upper bound on rad
  Pcu  = 0;   % Upper bound on P
  radl = Inf; % Lower bound on rad
  while (Pcl > alpha)
    Pcu  = Pcl;
    radl = radu; % Save previous result
    radu = radu*1.1;
    Pcl  = spm_P(1,radu^3   ,u,df,STAT,R,n,SS);
  end

  %
  % Newton-Rhapson refined search
  %
  d = 1;		    
  os = NaN;     % Old sign
  ms = (radu-radl)/10;  % Max step
  du = ms/100;
  % Linear interpolation for initial guess
  rad = radl*(alpha-Pcl)/(Pcu-Pcl)+radu*(Pcu-alpha)/(Pcu-Pcl);
  iter = 1;
  while abs(d) > epsP
    Pc  = spm_P(1,rad^3   ,u,df,STAT,R,n,SS);
    Pc1 = spm_P(1,(rad+du)^3,u,df,STAT,R,n,SS);
    d   = (alpha-Pc)/((Pc1-Pc)/du);
    os = sign(d);  % save old sign
    % Truncate search if step is too big
    if abs(d)>ms, 
      d = sign(d)*ms;
    end
    % Keep inside the given range
    if (rad+d)>radu, d = (radu-rad)/2; end
    if (rad+d)<radl, d = (rad-radl)/2; end
    % update
    rad = rad + d;
    iter = iter+1;
    if (iter>=maxi), 
      Status = 'TooManyIter';
      break; 
    end
  end
  % Convert back
  kr = rad^3;
  k = ceil(kr/v2r);
  Pc  = spm_P(1,k*v2r,u,df,STAT,R,n,SS);

%
% Brute force!
%
else
  Pc = 1;
  for k = range
    Pc = spm_P(1,k*v2r,u,df,STAT,R,n,SS);
    %fprintf('k=%d Pc=%g\n',k,Pc);
    if Pc <= alpha, 
      break; 
    end
  end;
  if (Pc > alpha)
    Status = 'OutOfRange';
  end
end
if ~nargout
    switch (Status)
     case {'JustPvalue'}
      fprintf(['  For a cluster-defining threshold of %0.4f a cluster size threshold of\n'...
           '  %d has corrected P-value %g\n\n'],...
          u,k,Pc);
     case {'OK'}
      fprintf(['  For a cluster-defining threshold of %0.4f the level %0.3f corrected\n'...
           '  cluster size threshold is %d and has size (corrected P-value) %g\n\n'],...
          u,alpha,k,Pc);
     case 'TooRough'
      fprintf(['\n  WARNING: Single voxel cluster is significant!\n\n',...
               '  For a cluster-defining threshold of %0.4f a k=1 voxel cluster\n'...
           '  size threshold has size (corrected P-value) %g\n\n'],...
          u,Pc); 
     case 'TooManyIter'
      fprintf(['\n  WARNING: Automated search failed to converge\n' ...
           '  Try systematic search.\n\n']); 
     case 'OutOfRange'  
      fprintf(['\n  WARNING: Within the range of cluster sizes searched (%g...%g)\n',...
             '  a corrected P-value <= alpha was not found (smallest P: %g)\n\n'],...
          range(1),range(end),Pc); 
      fprintf([  '  Try increasing the range or an automatic search.\n\n']); 
     otherwise
      error('Unknown status code');
    end
end
extent = k;
info.image = im;
info.extent = k;
info.alpha = alpha;
info.u = u;
info.Pc = Pc;
function success = bob_spm_save_rois_cluster(in, heightThresh, sizeThresh, mask)
% BOB_SPM_SAVE_CLUSTER
%
% USAGE: success = bob_spm_save_rois_cluster(in, heightThresh, sizeThresh, mask)
%
%   INPUTS
%       in:             image filename (full path if not in current dir)
%       heightThresh:   intensity threshold for defining clusters
%       sizeThresh:     extent threshold for defining clusters
%       mask:           optional mask filename (full path if not in 
%                       current directory)
%
%   OUTPUTS
%       success:        returns a 0 if no clusters could be identified
%                       after thresholding, 1 otherwise
%       
% ========================================================================%
if nargin<3, disp('USAGE: success = bob_spm_save_rois_cluster(in, heightThresh, sizeThresh, mask)'); return; end
if nargin<4, mask = []; end
    
% make sure image names are character arrays
% ------------------------------------------------------
if iscell(in), in = char(in); end;
maskflag = 1;
if isempty(mask) || length(mask)<1, maskflag = 0; mask = 'no mask';
elseif iscell(mask), mask = char(mask); end;

% write current images to command window
% ------------------------------------------------------
bob_display_message('Looking for Cluster');
fprintf('\nSource Image: %s\nMask Image: %s\nnHeight Threshold: %d\nSizeThreshold: %d\n', in, mask, heightThresh, sizeThresh);

% load images
% ------------------------------------------------------
in_hdr = spm_vol(in);
in = spm_read_vols(in_hdr);
imdims = size(in);
% if necessary, calculate critical t
if ismember(heightThresh,[.10 .05 .01 .005 .001 .0005 .0001]);
    tmp = in_hdr.descrip;
    idx1 = regexp(tmp,'[','ONCE');
    idx2 = regexp(tmp,']','ONCE');
    df = str2num(tmp(idx1+1:idx2-1));
    heightThresh = bob_p2t(heightThresh, df);
end
if maskflag
    mask = bob_reslice(mask,in_hdr.fname,1,1);
else
    mask = in>0;
end


% apply mask and height threshold
% ------------------------------------------------------
in(mask==0) = NaN;
in(in<heightThresh) = NaN;

% grab voxels
% ------------------------------------------------------
[X Y Z] = ind2sub(size(in), find(in > 0));
voxels = sortrows([X Y Z])';

% get cluster indices of voxels
% ------------------------------------------------------
cl_index = spm_clusters(voxels);
if isempty(cl_index)
	disp('No clusters found.')
    success = 0;
	return
end

% find index of voxels of sufficient size
% ------------------------------------------------------
cidx = unique(cl_index);
count = 0;
base_roi = zeros(imdims);
for i = cidx
    cluster_vox = voxels(:,cl_index==cidx(i));
    cluster_vox = cluster_vox';
    if length(cluster_vox)>=sizeThresh
        count = count + 1;
        cc = base_roi;
        for ii = 1:size(cluster_vox,1)
            cc(cluster_vox(ii,1),cluster_vox(ii,2),cluster_vox(ii,3)) = 1;
        end
        all_roi(:,:,:,count) = cc;
        all_size(count) = size(cluster_vox,1);
    end
end
if count==0
	disp('No clusters meet extent threshold.')
    success = 0;
	return
end
all_roi = double(all_roi);

% write the roi as an image
% ------------------------------------------------------
for i = 1:size(all_roi,4)
    roi_hdr = in_hdr;
    path = fileparts(in_hdr.fname);
    fulloutname = [path filesep 'ROI_Cluster' num2str(i) '_k=' num2str(all_size(i)) '.nii'];
    roi_hdr.fname = fulloutname;
    spm_write_vol(roi_hdr, all_roi(:,:,:,i));
    fprintf('\nROI of %d voxels written to:\n%s\n', all_size(i), [path filesep fulloutname]);
end
success = 1;
function [out, outmat] = bob_reslice(in, ref, int, nowrite)
% BOB_RESLICE 
%
% USAGE: [out M] = bob_reslice(in, ref, int, nowrite)
%
% ARGUMENTS
%   in: path to image to reslice
%   ref: path to reference image (image to reslice to)
%   int: interpolation method, 0=Nearest Neighbor, 1=Trilinear(default)
%   nowrite: option to not write new volume (default = 0)
%
% OUTPUT
%   out: the resliced image volume
%
% Most of the code is adapted from rest_Reslice in REST toolbox:
% Written by YAN Chao-Gan 090302 for DPARSF. Referenced from spm_reslice.
% State Key Laboratory of Cognitive Neuroscience and Learning 
% Beijing Normal University, China, 100875
% --------------------------------------------------------------------------
if nargin<4, nowrite = 0; end
if nargin<3, int = 1; end
if nargin<2, display('USAGE: out = bob_reslice(in, ref, int, nowrite)'); return; end
if iscell(in); in = char(in); end
if iscell(ref); ref = char(ref); end

% read in reference image
RefHead = spm_vol(ref); 
RefData = spm_read_vols(RefHead);
mat=RefHead.mat;
dim=RefHead.dim;

% read in image to reslice
SourceHead = spm_vol(in);
SourceData = spm_read_vols(SourceHead);

% do the reslicing
[x1,x2,x3] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
d     = [int*[1 1 1]' [1 1 0]'];
C = spm_bsplinc(SourceHead, d);
v = zeros(dim);
M = inv(SourceHead.mat)*mat; % M = inv(mat\SourceHead.mat) in spm_reslice.m
y1   = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
y2   = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
y3   = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
out    = spm_bsplins(C, y1,y2,y3, d);

%Revised by YAN Chao-Gan 121214. Apply a mask from the source image: don't extend values to outside brain.
tiny = 5e-2; % From spm_vol_utils.c
Mask = true(size(y1));
Mask = Mask & (y1 >= (1-tiny) & y1 <= (SourceHead.dim(1)+tiny));
Mask = Mask & (y2 >= (1-tiny) & y2 <= (SourceHead.dim(2)+tiny));
Mask = Mask & (y3 >= (1-tiny) & y3 <= (SourceHead.dim(3)+tiny));
out(~Mask) = 0;
outmat = mat;
if ~nowrite
    OutHead=SourceHead;
    OutHead.mat      = mat;
    OutHead.dim(1:3) = dim;
    [p n e] = fileparts(SourceHead.fname);
    newname = sprintf('%s_%dx%dx%d%s',n,dim,e);
    OutHead.fname = [p filesep newname];
    spm_write_vol(OutHead,out);
end

% | MISC UTILITIES
% =========================================================================
function vol = uigetvol(message, multitag)
% UIGETVOL Dialogue for selecting image volume file
%
%   USAGE: vol = uigetvol(message, multitag)
%       
%       message = to display to user
%       multitag = (default = 0) tag to allow selecting multiple images
%
% EX: img = uigetvol('Select Image to Process'); 
%
if nargin < 2, multitag = 0; end
if nargin < 1, message = 'Select Image File'; end
if ~multitag
    [imname, pname] = uigetfile({'*.img; *.nii', 'Image File'; '*.*', 'All Files (*.*)'}, message);
else
    [imname, pname] = uigetfile({'*.img; *.nii', 'Image File'; '*.*', 'All Files (*.*)'}, message, 'MultiSelect', 'on');
end
if isequal(imname,0) || isequal(pname,0)
    vol = [];
else
    vol = fullfile(pname, strcat(imname));
end
function vol = uiputvol(defname, prompt)
    if nargin < 1, defname = 'myimage.nii'; end
    if nargin < 2, prompt = 'Save image as'; end
    [imname, pname] = uiputfile({'*.img; *.nii', 'Image File'; '*.*', 'All Files (*.*)'}, prompt, defname);
    if isequal(imname,0) || isequal(pname,0)
        vol = [];
    else
        vol = fullfile(pname, imname); 
    end
function out = adjustbrightness(in)
    lim = .5;
    dat.min = min(in(in>0)); 
    dat.max = max(in(in>0));
    dat.dim = size(in);
    out = double(in)./255; 
    out(out>0) = out(out>0) + (lim-nanmean(nanmean(out(out>0))))*(1 - out(out>0)); 
    out(out>0) = scaledata(out(out>0), [dat.min dat.max]);
function fn = construct_filename
global st
[p,n]   = fileparts(st.ol.hdr.fname);
idx     = regexp(st.ol.descrip, ': ');
if ~isempty(idx)
    n = strtrim(st.ol.descrip(idx+1:end));
    n = regexprep(n, ' ', '_'); 
end
fn = sprintf('%s/%s_x=%d_y=%d_z=%d.png', p, n, round(st.centre));        
function handles = headsup(msg, wait4resp)
% HEADSUP Present message to user and wait for a response
%
%  USAGE: handles = headsup(msg, *wait4resp)    *optional input
% __________________________________________________________________________
%  INPUTS
%   msg: character array to present to user 
%

% ---------------------- Copyright (C) 2014 Bob Spunt ----------------------
%	Created:  2014-09-30
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 1, disp('USAGE: handles = headsup(msg, *wait4resp)'); return; end
if nargin < 2, wait4resp = 1; end
if ~iscell(msg), msg = char(msg); end
handles(1) = figure(...
    'Units', 'norm', ...
    'Position',[.425 .45 .15 .10],...
    'Resize','off',...
    'Color', [0.8941    0.1020    0.1098]*.60, ...
    'NumberTitle','off',...
    'DockControls','off',...
    'MenuBar','none',...
    'Name','Heads Up',...
    'Visible','on',...
    'Toolbar','none');
handles(2) = uicontrol('parent', handles(1), 'units', 'norm', 'style',  'text', 'backg', [0.8941    0.1020    0.1098]*.60,'foreg', [248/255 248/255 248/255], 'horiz', 'center', ...
    'pos', [.1 .40 .8 .40], 'fontname', 'arial', 'fontw', 'bold', 'fontsize', 16, 'string', msg, 'visible', 'on'); 
handles(3) = uicontrol('parent', handles(1), 'units', 'norm', 'style', 'push', 'foreg', [0 0 0], 'horiz', 'center', ...
    'pos', [.4 .10 .2 .30], 'fontname', 'arial', 'fontw', 'bold', 'fontsize', 16, 'string', 'OK', 'visible', 'on', 'callback', {@cb_ok, handles});
if wait4resp, uiwait(handles(1)); end
function cb_ok(varargin)
delete(varargin{:}); % Bye-bye figure
function [h, hCap] = layControl(parent,varargin)
%LAYCONTROL creates uicontrols within a figure's region or uipanel.
%
%% Description:
% Function LCONTROL creates diferent uicontrols un a figure's
% region or a uipanel, it starts filling the region with consecutive calls
% of the function using a left-to-right and top-to-bottom algorithm to
% ensure that all the uicontrols created are contained in the region or
% panel.
%
%% Use examples:
%
%           h = layControl(hParent,'align','Center','tall',20, ...
%                'width',100, 'separation',2, 'region',[.5 0.5  0.5 0.5]);
%
%           h = layControl(hparent,'panel',[0.5 0.5 0.35 0.35],...
%               'title','My panel', 'BackGroundColor',[.5 .5 .5]);
%
%   [h, hCap] = layControl(hParent,'nextLine','caption','Type here', ...
%               'BackGroundColor','b',...
%               'Style','Edit','String','Edit here');
%
%           h = layControl(hParent,'goBack',3,'Style','togglebutton');
%
%           h = layControl(hParent,'style','popupmenu',...
%               'string',{'One','Two'},...
%               'Callback',@changeSelection);
%
% *IMPORTANT INFORMATION:*
%    when creating a panel the function also creates a 'region' for 
%    creating new objects in there, but if you want thos objects to be 
%    aligned properly and actualy BE in the panel, you'll have to pass
%    the panel's handle as the 'parent' parameter to the function.
%
%    The diference between a 'caption' clause and a 'style','text' one is
%    that the 'caption' argument ensures that another object of at least
%    the same width of the caption can be placed in the same line, this
%    object can be any standard uicontrol, the handles of captions are
%    returned in the second output argument, while the 'Style','text' does
%    not ensures this and the handle is returned in the first output
%    argument.
%    
%
%
%% Input arguments:
%
% Following arguments MUST be specified first (at least once):
%
%      align: Especifies the default alignment of 'text' style uicontrol
%       tall: Default vertical size of the objects 
%      width: Default horizontal size of the objects
% separation: Separation betewen objects
%   nextLine: Jumps to the next line
%     goBack: Go back 'n' number of lines whithout reseting the 'x' coord.
%   dontWalk: Plots next object in the same 'x' coord.
% 
% Following arguments must be followed by the position:
%
%     region: creates a 'virtual' region and no objetc is created 
%      panel: creates a uipanel object in the region specified
%    
% Following arguments specify the object being created:
%
%      Style: the uicontrol 'Style' property e.g. 'edit' or 'popupmenu'
%    caption: creates a 'text' uicontrol but using the default alingment.
%
%% Output arguments:
%    h : the handle of object being created, 0 (zero) for argument 'region'
% hCap : the handle of captions being created.
%
%% Aditional information
%     Author: Roberto Cifuentes
%    Created: 2014-04-27
%Last Update: 2014-05-06

    persistent xMax yMax xMin  ultX ultY tall width separation align 
    if isempty(align)
        align = 'left' ;
    end
    dontWalk = false;
    k = 0 ;h = 0 ; hCap = 0 ; lastCap = false;
    while k < numel(varargin)
        k = k + 1 ;
        switch lower(varargin{k})
            case 'dontwalk' % Flag for not adding or restoring the default x coord
                dontWalk = true;
            case 'tall' % Default height for objects               
                k = k+1 ;
                tall = varargin{k}; 
            case 'width' % default width of objects
                k = k+1 ;
                width = varargin{k} ; 
            case 'separation' % Default spacer for objects
                k = k+1;
                separation = varargin{k} ;
            case 'nextline' % moves to next line with/without restoring the x coord.
                   if ~dontWalk
                        ultX = xMin + separation ;
                    end
                    ultY = ultY - tall - separation ;
            case 'goback' % Moves back N lines without restoring the x coord
                k = k+1 ;
                ultY = ultY+ (tall +separation)*varargin{k} ;
            case 'align' % default alignment for 'Caption's
                k = k+1 ;
                align =  varargin{k} ;               
            case 'region' % Creates a region normali the parent is a figure!!
                separatorFromPanelsTop = 25 ;% modify depending on figure's menu bar configuration
                k = k+1 ;
                curPos = varargin{k};
                regPos = get(parent,'Position') ;
                
                xMin = curPos(1) * regPos(3)+ separation;
                xMax = xMin + curPos(3)*regPos(3) - separation;   
                yMax = curPos(2)*regPos(4) + curPos(4)* regPos(4) - separatorFromPanelsTop ;
                
                ultX = xMin + separation ; 
                ultY = yMax - separation ;                
            case 'panel'
                separatorFromPanelsTop = 40 ;% modify depending on figure's menu bar configuration
                k = k+1 ;
                curPos = varargin{k};
                h(end+1,1) = uipanel(parent,'Position',curPos,'Title','');                
                regPos = get(parent,'Position') ;
                
                xMin = separation;
                xMax = curPos(3)*regPos(3) - separation;   
                yMax = curPos(4)*regPos(4) - separatorFromPanelsTop ;
                
                ultX = xMin + separation ; 
                ultY = yMax - separation ;                

            case 'caption'
                lastCap = true;
                k = k+1 ;
                if ultX+separation+2*width >= xMax
                    if ~dontWalk
                        ultX = xMin + separation ;
                    end
                    ultY = ultY - tall - separation ;
                end
                posicion = [ ultX+separation, ultY, width, tall];
                if ~dontWalk
                    ultX = posicion(1) + width ;
                end
                hCap(end+1,1) = uicontrol(parent,'Style','text', ...
                    'String',varargin{k}, ...
                    'HorizontalAlignment',align, ...
                    'Position', posicion ...
                    ) ;
                dontWalk = false ;
            case 'style'
                lastCap = false;
                k = k+1 ;
        
                if ultX+separation+width >= xMax
                   if ~dontWalk
                        ultX = xMin + separation ;
                   end                    
                    ultY = ultY - tall - separation ;
                end
                posicion = [ultX+separation ultY width tall];
                if ~dontWalk
                    ultX = posicion(1) + width ;
                end
                h(end+1,1) = uicontrol(parent,'Style',varargin{k}, ...
                    'Position',posicion);
                switch lower(varargin{k})
                    case {'edit','listbox','popupmenu'}
                        set(h(end,1),'BackGroundColor','w');
                end
                 dontWalk = false;
            otherwise
               if lastCap
                   set(hCap(end,1),varargin{k},varargin{k+1});
               else
                   set(h(end,1),varargin{k},varargin{k+1});
               end
                k = k+1 ;
        end

%         fprintf('\n xMin: %6.2f  ultX: %6.2f ultY: %6.2f ',xMin,ultX,ultY);
    end
    
    if numel(h)>1
        h = h(2:end,1) ;
    end
    if numel(hCap)>1
        hCap = hCap(2:end,1);
    end

% | MAXIMUM INTENSITY PROJECTION (MIP; FROM SPM8)
% =========================================================================
function addmip(varargin)
%% MIP

% Add MIP to Panel
hMIPax  = axes('Parent',hReg,'Position',[.05 .05 .9 .9],'Visible','off');

%% MIP
hMIPax  = bspm_mip_ui(Z,XYZmm,M,DIM,hMIPax,units);

%% Cross Register the MIP with Registry, Push Coords
bspm_XYZreg('XReg',hReg, hMIPax,'bspm_mip_ui');
function varargout = bspm_mip_ui(varargin)
%-Condition arguments
%==========================================================================
if nargin==0
    error('Insufficient arguments')
elseif ~ischar(varargin{1})
    varargout={bspm_mip_ui('Display',varargin{1:end})}; return
end
%-Axis offsets for 3d MIPs:
%==========================================================================
%-MIP pane dimensions, origin offsets and #pixels per mm
%-See spm_project.c for derivation
mipmat = char(spm_get_defaults('stats.results.mipmat'));
load(mipmat, 'DXYZ', 'CXYZ', 'scale');
Po(1)  =                  CXYZ(2) -2;
Po(2)  = DXYZ(3)+DXYZ(1) -CXYZ(1) +2;
Po(3)  = DXYZ(3)         -CXYZ(3) +2;
Po(4)  = DXYZ(2)         +CXYZ(1) -2;
%==========================================================================
switch lower(varargin{1}), case 'display'
%==========================================================================
    % hMIPax = bspm_mip_ui('Display',Z,XYZ,M,DIM,F,units)
    
    if nargin<5
        F      = gcf;
        hMIPax = [];
    else
        F = varargin{6};
        if ischar(F), F=spm_figure('FindWin',F); end
        if ~ishandle(F), error('Invalid handle'), end
        switch get(F,'Type'), case 'figure'
            hMIPax = [];
            case 'axes'
                hMIPax = F;
                F      = get(hMIPax,'Parent');
            otherwise
                error('F not a figure or axis handle')
        end
    end
    if nargin<4, error('Insufficient arguments'), end
    Z       = varargin{2};
    XYZ     = varargin{3};
    M       = varargin{4};
    DIM     = varargin{5};
    try
        units = varargin{7};
    catch
        units = {'mm' 'mm' 'mm'};
    end

    xyz = bspm_XYZreg('RoundCoords',[0;0;0],M,DIM);

    %-Display (MIP) transformation matrix
    %----------------------------------------------------------------------
    Md      = eye(4);
    Ms      = diag([scale(1:3) 1]);
    
    %-Display MIP
    %----------------------------------------------------------------------
    set(F,'Units','normalized')
    pXYZ = Ms*Md*[XYZ;ones(1,size(XYZ,2))];
    bspm_mip(Z,pXYZ(1:3,:),Ms*Md*M,units);
    hMIPim = get(gca,'Children');

    %-Print coordinates
    %----------------------------------------------------------------------
    hMIPxyz = text(0,max(get(hMIPax,'YLim'))/2,...
        {'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]',xyz(1:3))},...
        'Interpreter','TeX','FontName',spm_platform('font','times'),...
        'Color',[1,1,1]*.7,...
        'HorizontalAlignment','Center',...
        'VerticalAlignment','Bottom',...
        'Rotation',90,...
        'Tag','hMIPxyz',...
        'UserData',xyz);
    
    %-Create point markers
    %----------------------------------------------------------------------
    xyz = Ms*Md*[xyz(:);1];
    xyz = xyz(1:3);
    hX1r  = text(Po(1)+xyz(2),Po(2)+xyz(1),'<',...
        'Color','r','Fontsize',20, 'FontWeight', 'bold',...
        'Tag','hX1r',...
        'ButtonDownFcn','bspm_mip_ui(''MoveStart'')');
    hX2r  = text(Po(1)+xyz(2),Po(3)-xyz(3),'<',...
        'Color','r','Fontsize',20, 'FontWeight', 'bold',...
        'Tag','hX2r',...
        'ButtonDownFcn','bspm_mip_ui(''MoveStart'')');
    hX3r  = text(Po(4)+xyz(1),Po(3)-xyz(3),'<',...
        'Color','r','Fontsize',20, 'FontWeight', 'bold',...
        'Tag','hX3r',...
        'ButtonDownFcn','bspm_mip_ui(''MoveStart'')');
    hXr   = [hX1r,hX2r,hX3r];


    if DIM(3) == 1
        %-2 dimensional data
        %------------------------------------------------------------------
        set(hXr(3),'Visible','off');
        set(hXr(2),'Visible','off');

    end

    %-Create UIContextMenu for marker jumping
    %-----------------------------------------------------------------------
    h = uicontextmenu('Tag','MIPconmen','UserData',hMIPax);
    if isempty(XYZ), str='off'; else str='on'; end
    uimenu(h,'Separator','off','Label','goto nearest suprathreshold voxel',...
        'CallBack',['bspm_mip_ui(''Jump'',',...
        'get(get(gcbo,''Parent''),''UserData''),''nrvox'');'],...
        'Interruptible','off','BusyAction','Cancel','Enable',str);
    uimenu(h,'Separator','off','Label','goto nearest local maximum',...
        'CallBack',['bspm_mip_ui(''Jump'',',...
        'get(get(gcbo,''Parent''),''UserData''),''nrmax'');'],...
        'Interruptible','off','BusyAction','Cancel','Enable',str);
    uimenu(h,'Separator','off','Label','goto global maximum',...
        'CallBack',['bspm_mip_ui(''Jump'',',...
        'get(get(gcbo,''Parent''),''UserData''),''glmax'');'],...
        'Interruptible','off','BusyAction','Cancel','Enable',str);
    uimenu(h,'Separator','on','Label','save MIP as...',...
        'CallBack',['bspm_mip_ui(''Save'', ',...
        'get(get(gcbo,''Parent''),''UserData''));'],...
        'Interruptible','off','BusyAction','Cancel');
    h1 = uimenu(h,'Separator','on','Label','Extract betas');
    uimenu(h1,'Label','This voxel',...
        'CallBack','beta=bspm_mip_ui(''Extract'', ''voxel'')',...
        'Interruptible','off','BusyAction','Cancel');
    uimenu(h1,'Label','This cluster',...
        'CallBack','beta=bspm_mip_ui(''Extract'', ''cluster'')',...
        'Interruptible','off','BusyAction','Cancel');

    set(hMIPim,'UIContextMenu',h)

    %-Save handles and data
    %----------------------------------------------------------------------
    set(hMIPax,'Tag','hMIPax','UserData',...
        struct(...
        'hReg',     [],...
        'hMIPxyz',  hMIPxyz,...
        'XYZ',      XYZ,...
        'Z',        Z,...
        'M',        M,...
        'Md',       Md,...
        'Ms',       Ms,...
        'DIM',      DIM,...
        'units',    {units},...
        'hXr',      hXr))

    varargout = {hMIPax};


    %======================================================================
    case 'getcoords'
    %======================================================================
        % xyz = bspm_mip_ui('GetCoords',h)
        if nargin<2, h=bspm_mip_ui('FindMIPax'); else h=varargin{2}; end
        varargout = {get(getfield(get(h,'UserData'),'hMIPxyz')  ,'UserData')};



    %======================================================================
    case 'setcoords'
    %======================================================================
        % [xyz,d] = bspm_mip_ui('SetCoords',xyz,h,hC)
        if nargin<4, hC=0; else hC=varargin{4}; end
        if nargin<3, h=bspm_mip_ui('FindMIPax'); else h=varargin{3}; end
        if nargin<2, error('Set coords to what?'), else xyz=varargin{2}; end

        MD  = get(h,'UserData');

        %-Check validity of coords only when called without a caller handle
        %------------------------------------------------------------------
        if hC<=0
            [xyz,d] = bspm_XYZreg('RoundCoords',xyz,MD.M,MD.DIM);
            if d>0 && nargout<2, warning(sprintf(...
                    '%s: Co-ords rounded to nearest voxel center: Discrepancy %.2f',...
                    mfilename,d));
            end
        else
            d = [];
        end

        %-Move marker points, update internal cache in hMIPxyz
        %------------------------------------------------------------------
        bspm_mip_ui('PosnMarkerPoints',MD.Md(1:3,:)*[xyz;1],h);
        set(MD.hMIPxyz,'UserData',reshape(xyz(1:3),3,1))
        set(MD.hMIPxyz,'String',{'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]',xyz(1:3))})

        %-Tell the registry, if we've not been called by the registry...
        %------------------------------------------------------------------
        if ~isempty(MD.hReg) && MD.hReg~=hC, bspm_XYZreg('SetCoords',xyz,MD.hReg,h); end

        %-Return arguments
        %------------------------------------------------------------------
        varargout = {xyz,d};



    %======================================================================
    case 'posnmarkerpoints'
    %======================================================================
        % bspm_mip_ui('PosnMarkerPoints',xyz,h)
        if nargin<3, h=bspm_mip_ui('FindMIPax'); else h=varargin{3}; end
        if nargin<2, xyz = bspm_mip_ui('GetCoords',h); else xyz = varargin{2}; end
        MD = get(h,'UserData');
        
        %-Get handles of marker points from UserData of hMIPax
        %------------------------------------------------------------------
        hX = MD.hXr;

        %-Set marker points
        %------------------------------------------------------------------
        set(hX,'Units','Data')
        if length(hX)==1
            vx  = sqrt(sum(MD.M(1:3,1:3).^2));
            tmp = MD.M\[xyz ; 1];
            tmp = tmp(1:2).*vx(1:2)';
            set(hX,'Position',[tmp(1), tmp(2), 1])
        else
            pxyz = MD.Ms*[xyz(1:3);1];
            set(hX(1),'Position',[ Po(1) + pxyz(2), Po(2) + pxyz(1), 0])
            set(hX(2),'Position',[ Po(1) + pxyz(2), Po(3) - pxyz(3), 0])
            set(hX(3),'Position',[ Po(4) + pxyz(1), Po(3) - pxyz(3), 0])
        end


    %======================================================================
    case 'jump'
    %======================================================================
        % [xyz,d] = bspm_mip_ui('Jump',h,loc)
        if nargin<3, loc='nrvox'; else loc=varargin{3}; end
        if nargin<2, h=bspm_mip_ui('FindMIPax'); else h=varargin{2}; end

        %-Get current location & MipData
        %------------------------------------------------------------------
        oxyz = bspm_mip_ui('GetCoords',h);
        MD   = get(h,'UserData');


        %-Compute location to jump to
        %------------------------------------------------------------------
        if isempty(MD.XYZ), loc='dntmv'; end
        switch lower(loc), case 'dntmv'
            spm('alert!','No suprathreshold voxels to jump to!',mfilename,0);
            varargout = {oxyz, 0};
            return
            case 'nrvox'
                str       = 'nearest suprathreshold voxel';
                [xyz,i,d] = bspm_XYZreg('NearestXYZ',oxyz,MD.XYZ);
            case 'nrmax'
                str       = 'nearest local maximum';
                iM        = inv(MD.M);
                XYZvox    = iM(1:3,:)*[MD.XYZ; ones(1,size(MD.XYZ,2))];
                [null,null,XYZvox] = spm_max(MD.Z,XYZvox);
                XYZ       = MD.M(1:3,:)*[XYZvox; ones(1,size(XYZvox,2))];
                [xyz,i,d] = bspm_XYZreg('NearestXYZ',oxyz,XYZ);
            case 'glmax'
                str       = 'global maximum';
                [null, i] = max(MD.Z); i = i(1);
                xyz       = MD.XYZ(:,i);
                d         = sqrt(sum((oxyz-xyz).^2));
            case 'nrchan'
                str       = 'nearest suprathreshold channel';
                if ~isfield(MD, 'hChanPlot'), bspm_mip_ui('Channels',h); end
                MD        = get(h,'UserData');
                [xyz,i,d] = bspm_XYZreg('NearestXYZ',[oxyz(1); oxyz(2); 0],MD.Channels.pos);
                xyz(3)    = oxyz(3);
                str       = [str sprintf(' (%s)',MD.Channels.name{i})];
            otherwise
                warning('Unknown jumpmode')
                varargout = {xyz,0};
                return
        end

        %-Write jump report, jump, and return arguments
        %------------------------------------------------------------------
        fprintf(['\n\t%s:\tJumped %0.2fmm from [%3.0f, %3.0f, %3.0f],\n\t\t\t',...
            'to %s at [%3.0f, %3.0f, %3.0f]\n'],...
            mfilename, d, oxyz, str, xyz)

        bspm_mip_ui('SetCoords',xyz,h,h);
        varargout = {xyz, d};


    %======================================================================
    case 'findmipax'
    %======================================================================
        % hMIPax = bspm_mip_ui('FindMIPax',h)
        % Checks / finds hMIPax handles
        %-**** h is handle of hMIPax, or figure containing MIP (default gcf)
        if nargin<2, h=get(0,'CurrentFigure'); else h=varargin{2}; end
        if ischar(h), h=spm_figure('FindWin',h); end
        if ~ishandle(h), error('invalid handle'), end
        if ~strcmp(get(h,'Tag'),'hMIPax'), h=findobj(h,'Tag','hMIPax'); end
        if isempty(h), error('MIP axes not found'), end
        if length(h)>1, error('Multiple MIPs in this figure'), end
        varargout = {h};



    %======================================================================
    case 'movestart'
    %======================================================================
        % bspm_mip_ui('MoveStart')
        [cO,cF] = gcbo;
        hMIPax  = get(cO,'Parent');
        MD      = get(hMIPax,'UserData');

        %-Store useful quantities in UserData of gcbo, the object to be dragged
        %------------------------------------------------------------------
        set(hMIPax,'Units','Pixels')
        set(cO,'UserData',struct(...
            'hReg',     MD.hReg,...
            'xyz',      bspm_mip_ui('GetCoords',hMIPax),...
            'MIPaxPos', get(hMIPax,'Position')*[1,0;0,1;0,0;0,0],...
            'hMIPxyz',  MD.hMIPxyz,...
            'M',        MD.M,...
            'Md',       MD.Md,...
            'Ms',       MD.Ms,...
            'DIM',      MD.DIM,...
            'hX',       MD.hXr))

        %-Initiate dragging
        %------------------------------------------------------------------
        if strcmp(get(cF,'SelectionType'),'normal') || isempty(MD.XYZ)
            %-Set Figure callbacks for drop but no drag (DragType 0)
            %--------------------------------------------------------------
            set(MD.hMIPxyz,'Visible','on','String',...
                {'{\bfSPM}{\itmip}','\itPoint & drop...'})
            set(cF,'WindowButtonUpFcn',    'bspm_mip_ui(''Move'',0)',...
                'Interruptible','off')
            set(cF,'Pointer','CrossHair')
            
       %-Set Figure callbacks for drag'n'drop (DragType 1)
       %-------------------------------------------------------------------
        elseif strcmp(get(cF,'SelectionType'),'extend')
            set(MD.hMIPxyz,'Visible','on','String',...
                {'{\bfSPM}{\itmip}','\itDynamic drag & drop...'})
            set(cF,'WindowButtonMotionFcn','bspm_mip_ui(''Move'',1)',...
                'Interruptible','off')
            set(cF,'WindowButtonUpFcn',    'bspm_mip_ui(''MoveEnd'')',...
                'Interruptible','off')
            set(cF,'Pointer','Fleur')
            
        %-Set Figure callbacks for drag'n'drop with co-ord updating (DragType 2)
        %------------------------------------------------------------------
        elseif strcmp(get(cF,'SelectionType'),'alt')
            set(MD.hMIPxyz,'Visible','on','String',...
                {'{\bfSPM}{\itmip}','\itMagnetic drag & drop...'})
            set(cF,'WindowButtonMotionFcn','bspm_mip_ui(''Move'',2)',...
                'Interruptible','off')
            set(cF,'WindowButtonUpFcn',    'bspm_mip_ui(''MoveEnd'')',...
                'Interruptible','off')
            set(cF,'Pointer','Fleur')
        end



    %======================================================================
    case 'move'
    %======================================================================
        % bspm_mip_ui('Move',DragType)
        if nargin<2, DragType = 2; else DragType = varargin{2}; end
        cF = gcbf;
        cO = gco(cF);

        %-Get useful data from UserData of gcbo, the object to be dragged
        %------------------------------------------------------------------
        MS  = get(cO,'UserData');

        %-Work out where we are moving to - Use HandleGraphics to give position
        %------------------------------------------------------------------
        set(cF,'Units','pixels')
        d = get(cF,'CurrentPoint') - MS.MIPaxPos;
        set(cO,'Units','pixels')
        set(cO,'Position',d)
        set(cO,'Units','data')
        d = get(cO,'Position');

        %-Work out xyz, depending on which view is being manipulated
        %------------------------------------------------------------------
        sMarker = get(cO,'Tag');
        pxyz = MS.Ms*[MS.xyz;1];
        if strcmp(sMarker,'hX1r')
            xyz = [d(2) - Po(2); d(1) - Po(1); pxyz(3)];
        elseif strcmp(sMarker,'hX2r')
            xyz = [pxyz(1); d(1) - Po(1); Po(3) - d(2)];
        elseif strcmp(sMarker,'hX3r')
            xyz = [d(1) - Po(4); pxyz(2); Po(3) - d(2)];
        else
            error('Can''t work out which marker point')
        end
        xyz = inv(MS.Ms*MS.Md) * [xyz(:);1]; xyz = xyz(1:3);
        
        %-Round coordinates according to DragType & set in hMIPxyz's UserData
        %------------------------------------------------------------------
        if DragType==0
            xyz    = bspm_XYZreg('RoundCoords',xyz,MS.M,MS.DIM);
        elseif DragType==1
            xyz    = bspm_XYZreg('RoundCoords',xyz,MS.M,MS.DIM);
            hMIPax = get(cO,'Parent');
            MD     = get(hMIPax,'UserData');
            i      = bspm_XYZreg('FindXYZ',xyz,MD.XYZ);
        elseif DragType==2
            hMIPax = get(cO,'Parent');
            MD     = get(hMIPax,'UserData');
            [xyz,i,d] = bspm_XYZreg('NearestXYZ',xyz,MD.XYZ);
        end
        set(MS.hMIPxyz,'UserData',xyz)

        %-Move marker points
        %------------------------------------------------------------------
        set(MS.hX,'Units','Data')
        xyz2 = MS.Ms * MS.Md * [xyz(:);1];
        xyz2 = xyz2(1:3);
        set(MS.hX(1),'Position',[ Po(1) + xyz2(2), Po(2) + xyz2(1), 0])
        set(MS.hX(2),'Position',[ Po(1) + xyz2(2), Po(3) - xyz2(3), 0])
        set(MS.hX(3),'Position',[ Po(4) + xyz2(1), Po(3) - xyz2(3), 0])


        %-Update dynamic co-ordinate strings (if appropriate DragType)
        %------------------------------------------------------------------
        if DragType==0
            bspm_mip_ui('MoveEnd')
        elseif DragType==1
            if isempty(i)
                str = {'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]',xyz)};
            else
                str = {'{\bfSPM}{\itmip}: ',...
                    sprintf('[%g, %g, %g]: %.4f',xyz,MD.Z(i))};
            end
            set(MD.hMIPxyz,'String',str)
        elseif DragType==2
            set(MD.hMIPxyz,'String',...
                {'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]: %.4f',xyz,MD.Z(i))})
        else
            error('Illegal DragType')
        end



    %======================================================================
    case 'moveend'
    %======================================================================
        % bspm_mip_ui('MoveEnd')
        cF = gcbf;
        cO = gco(cF);
        hMIPax  = get(cO,'Parent');
        MS      = get(cO,'UserData');

        %-Reset WindowButton functions, pointer & SPMmip info-string
        %------------------------------------------------------------------
        set(gcbf,'WindowButtonMotionFcn',' ')
        set(gcbf,'WindowButtonUpFcn',' ')
        set(gcbf,'Pointer','arrow')
        set(MS.hMIPxyz,'String',...
            {'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]',get(MS.hMIPxyz,'UserData'))})

        %-Set coordinates after drag'n'drop, tell registry
        %------------------------------------------------------------------
        % don't need to set internal coordinates 'cos 'move' does that
        if ~isempty(MS.hReg)
            bspm_XYZreg('SetCoords',get(MS.hMIPxyz,'UserData'),MS.hReg,hMIPax);
        end

    %======================================================================
    case 'channels'
    %======================================================================
        % bspm_mip_ui('Channels',h) 
        % M/EEG specific to display channels on 1-slice MIP

        if nargin<2, h=bspm_mip_ui('FindMIPax'); else h=varargin{2}; end
        MD  = get(h,'UserData');

        if ~isfield(MD, 'hChanPlot')
            % first time call
            D = spm_eeg_load;
            if isempty(D), return; end

            DIM = get(findobj('Tag','hFxyz'), 'UserData');

            [mod, Cind] = spm_eeg_modality_ui(D, 1, 1);
            otherind = setdiff(1:nchannels(D), Cind);
            if ~isempty(otherind)
                D = chantype(D, otherind, 'Other');
            end
            [Cel x, y] = spm_eeg_locate_channels(D, DIM.DIM(1), Cind);
            Cel = DIM.M * [Cel'; ones(2,size(Cel,1))];
            Cel = Cel(1:2,:)';
            pos = [Cel'; zeros(1,size(Cel,1))];
            Cel(:,1) = Cel(:,1) + Po(2);
            Cel(:,2) = Cel(:,2) + Po(1);
            hold(h,'on'), hChanPlot = plot(h,Cel(:, 2), Cel(:, 1), 'b*');

            hChanText = cell(1,size(Cel,1));
            name = D.chanlabels(Cind);
            figure(spm_figure('FindWin'));
            for i = 1:size(Cel, 1)
                hChanText{i} = text(Cel(i, 2)+0.5, Cel(i, 1), name{i}, 'Color', 'b');
            end

            MD.hChanPlot     = hChanPlot;
            MD.hChanText     = hChanText;
            MD.Channels.pos  = pos;
            MD.Channels.name = D.chanlabels(Cind);
            set(h, 'UserData', MD);
        else
            if strcmp(get(MD.hChanPlot, 'Visible'), 'on');
                % switch it off
                set(MD.hChanPlot, 'Visible', 'off');
                for i = 1:length(MD.hChanText)
                    set(MD.hChanText{i}, 'Visible', 'off');
                end
            else
                % switch it on
                set(MD.hChanPlot, 'Visible', 'on');
                for i = 1:length(MD.hChanText)
                    set(MD.hChanText{i}, 'Visible', 'on');
                end
            end

        end

        
    %======================================================================
    case 'save'
    %======================================================================
        % bspm_mip_ui('Save',h) 
        if nargin<2, h=bspm_mip_ui('FindMIPax'); else h=varargin{2}; end
        %MD  = get(h,'UserData');
        %pXYZ = MD.Ms*MD.Md*[MD.XYZ;ones(1,size(MD.XYZ,2))];
        %mip  = spm_mip(MD.Z,pXYZ(1:3,:),MD.Ms*MD.Md*MD.M,MD.units);
        mip = get(findobj(h,'Type','image'),'CData');
        [f, p] = uiputfile('*.png', 'Save as', 'mip.png');
        if ~isequal(f,0) && ~isequal(p,0)
            if ndims(mip) == 3
                imwrite(mip,fullfile(p,f),'png');
            else
                imwrite(mip,gray(64),fullfile(p,f),'png');
            end
            fprintf('Saving SPM MIP in %s\n',fullfile(p,f));
        end
    
    
    %======================================================================
    case 'extract'
    %======================================================================
        % beta = bspm_mip_ui('Extract',action)
        if nargin<2, action='voxel'; else action=varargin{2}; end
        xSPM = evalin('base','xSPM');
        SPM  = evalin('base','SPM');
        
        XYZmm = spm_results_ui('GetCoords');
        [XYZmm,i] = bspm_XYZreg('NearestXYZ',XYZmm,xSPM.XYZmm);
        spm_results_ui('SetCoords',xSPM.XYZmm(:,i));
        
        switch lower(action)
            case 'voxel'
                % current voxel
                XYZ = SPM.xVol.iM(1:3,:)*[XYZmm;1];
                
            case 'cluster'
                % current cluster
                A   = spm_clusters(xSPM.XYZ);
                j   = find(A == A(i));
                XYZ = xSPM.XYZ(:,j);
                
            otherwise
                error('Unknown action.');
        end
        
        beta = spm_get_data(SPM.Vbeta,XYZ);
        
        varargout = {beta};
        
        
    %======================================================================
    otherwise
    %======================================================================
        error('Unknown action string')

end
function mip = bspm_mip(Z,XYZ,M,units)
% SPM Maximum Intensity Projection
% FORMAT mip = spm_mip(Z,XYZ,M,units)
% Z       - vector point list of SPM values for MIP
% XYZ     - matrix of coordinates of points (mip coordinates)
% M       - voxels - > mip matrix or size of voxels (mm)
% units   - defining space     [default {'mm' 'mm' 'mm'}]
%
% mip     - maximum intensity projection
%           if no output, the mip is displayed in current figure.
%__________________________________________________________________________
%
% If the data are 2 dimensional [DIM(3) = 1] the projection is simply an
% image, otherwise:
%
% spm_mip creates and displays a maximum intensity projection of a point
% list of voxel values (Z) and their location (XYZ) in three orthogonal
% views of the brain. It is assumed voxel locations conform to the space
% defined in the atlas of Talairach and Tournoux (1988); unless the third
% dimension is time.
%
% This routine loads a mip outline from MIP.mat. This is an image with
% contours and grids defining the space of Talairach & Tournoux (1988).
% mip95 corresponds to the Talairach atlas, mip96 to the MNI templates.
% The outline and grid are superimposed at intensity 0.4.
%
% A customised mip outline can be used instead of the default.
%
% A default colormap of 64 levels is assumed. The pointlist image is
% scaled to fit in the interval [1/9,1]*64 for display. Flat images
% are scaled to 1*64.
%
% If M is not specified, it is assumed the XYZ locations are 
% in Talairach mm.
%__________________________________________________________________________
% Copyright (C) 1996-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_mip.m 5245 2013-02-06 17:28:06Z guillaume $

%-Get units and grid scaling
%--------------------------------------------------------------------------
try, units; catch, units = {'mm' 'mm' 'mm'}; end
try, M;     catch, M = 1; end
Grid = 0.4;

%-Transpose locations if necessary
%--------------------------------------------------------------------------
if size(XYZ,1) ~= 3, XYZ = XYZ';         end
if size(Z,1)   ~= 1, Z   = Z';           end
if size(M,1)   == 1, M   = speye(4,4)*M; end

%-Scale & offset point list values to fit in [1/(1+Scal),1]
%--------------------------------------------------------------------------
Z    = Z - min(Z);
mx   = max(Z);
Scal = 8;
if isempty(mx)
    Z = [];
elseif isfinite(mx) && mx
    Z = (1 + Scal*Z/mx)/(Scal + 1);
else
    Z = ones(1,length(Z));
end
%-Display format
%==========================================================================
% load various grids, DXYZ, CXYZ, scale (see spm_mip_ui and spm_project)
load(char(spm_get_defaults('stats.results.mipmat')));

%-Single slice case
%--------------------------------------------------------------------------
if isempty(units{3}) && ~strcmp(units{2},'mm')

    %-2d case: Time-Frequency or Frequency-Frequency
    %----------------------------------------------------------------------
    mip = 4*grid_trans;
      
elseif isempty(units{3})
    %-2d case
    %----------------------------------------------------------------------
    mip = 4*grid_trans + mask_trans; 
elseif strcmp(units{3},'ms') || strcmp(units{3},'Hz')
    
    %-3d case: Space-time
    %----------------------------------------------------------------------
    mip = 4*grid_time + mask_trans;

else
    %-3d case: Space
    %----------------------------------------------------------------------
    mip = 4*grid_all + mask_all;
end
%-Create maximum intensity projection
%--------------------------------------------------------------------------
mip  = mip/max(mip(:));
c    = [0 0 0 ;
        0 0 1 ;
        0 1 0 ;
        0 1 1 ;
        1 0 0 ;
        1 0 1 ; 
        1 1 0 ; 
        1 1 1 ] - 0.5;
c    = c*M(1:3,1:3);
dim  = [(max(c) - min(c)) size(mip)];
d    = spm_project(Z,round(XYZ),dim,DXYZ,CXYZ);
mip  = max(d,Grid*mip);
mip  = rot90((1 - mip)*64);
%-And display it
%--------------------------------------------------------------------------
if ~nargout
    image(mip); axis tight; axis off;
end

% | SPM_OPTHVIEWS (MODIFIED FROM ORIGINAL SPM8 CODE)
% =========================================================================
function varargout = bspm_orthviews(action,varargin)
% John Ashburner et al% Display orthogonal views of a set of images
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
%         st.registry.hReg \_ See bspm_XYZreg for documentation
%         st.registry.hMe  /
%
% For each of the displayed images, there is a non-empty entry in the
% vols cell array.  Handles returned by "spm_orthviews('Image',.....)"
% indicate the position in the cell array of the newly created ortho-view.
% Operations on each ortho-view require the handle to be passed.
%
% When a new image is displayed, the cell entry contains the information
% returned by spm_vol (type help spm_vol for more info).  In addition,
% there are a few other fields, some of which are documented here:
%
%         premul  - a matrix to premultiply the .mat field by.  Useful
%                   for re-orienting images.
%         window  - either 'auto' or an intensity range to display the
%                   image with.
%         mapping - Mapping of image intensities to grey values. Currently
%                   one of 'linear', 'histeq', loghisteq',
%                   'quadhisteq'. Default is 'linear'.
%                   Histogram equalisation depends on the image toolbox
%                   and is only available if there is a license available
%                   for it.
%         ax      - a cell array containing an element for the three
%                   views.  The fields of each element are handles for
%                   the axis, image and crosshairs.
%
%         blobs   - optional.  Is there for using to superimpose blobs.
%                   vol     - 3D array of image data
%                   mat     - a mapping from vox-to-mm (see spm_vol, or
%                             help on image formats).
%                   max     - maximum intensity for scaling to.  If it
%                             does not exist, then images are auto-scaled.
%
%                   There are two colouring modes: full colour, and split
%                   colour.  When using full colour, there should be a
%                   'colour' field for each cell element.  When using
%                   split colourscale, there is a handle for the colorbar
%                   axis.
%
%                   colour  - if it exists it contains the
%                             red,green,blue that the blobs should be
%                             displayed in.
%                   cbar    - handle for colorbar (for split colourscale).
global st

persistent zoomlist reslist

if isempty(st), reset_st; end

if ~nargin, action = ''; end

if ~any(strcmpi(action,{'reposition','pos'}))
    spm('Pointer','Watch');
end
    
switch lower(action)
    case 'image'
        H = specify_image(varargin{1});
        
        if ~isempty(H)
            if numel(varargin)>=2
                st.vols{H}.area = varargin{2};
            else
                st.vols{H}.area = [0 0 1 1];
            end
            if isempty(st.bb), st.bb = maxbb; end
            resolution;
            bbox;
            cm_pos;
        end
        varargout{1} = H;
        mmcentre     = mean(st.Space*[maxbb';1 1],2)';
        st.centre    = mmcentre(1:3);
        redraw_all

    case 'caption'
        if ~isnumeric(varargin{1})
            varargin{1} = cellstr(varargin{1});
            xlh = NaN(numel(varargin{1}),1);
            for i=1:numel(varargin{1})
                h = bspm_orthviews('Caption',i,varargin{1}{i},varargin{3:end});
                if ~isempty(h), xlh(i) = h; end
            end
            varargout{1} = xlh;
            return;
        end
        
        vh = valid_handles(varargin{1});
        nh = numel(vh);
        
        xlh = nan(nh, 1);
        for i = 1:nh
            xlh(i) = get(st.vols{vh(i)}.ax{3}.ax, 'XLabel');
            if iscell(varargin{2})
                if i <= length(varargin{2})
                    set(xlh(i), 'String', varargin{2}{i});
                end
            else
                set(xlh(i), 'String', varargin{2});
            end
            for np = 4:2:nargin
                property = varargin{np-1};
                value = varargin{np};
                set(xlh(i), property, value);
            end
        end
        varargout{1} = xlh;
        
    case 'bb'
        if ~isempty(varargin) && all(size(varargin{1})==[2 3]), st.bb = varargin{1}; end
        bbox;
        redraw_all;
        
    case 'redraw'
        redraw_all;
        callback;
        if isfield(st,'registry')
            bspm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
        end
        
    case 'reload_mats'
        if nargin > 1
            handles = valid_handles(varargin{1});
        else
            handles = valid_handles;
        end
        for i = handles
            fnm = spm_file(st.vols{i}.fname, 'number', st.vols{i}.n);
            st.vols{i}.mat = spm_get_space(fnm);
        end
        % redraw_all (done in bspm_orthviews('reorient','context_quit'))
        
    case 'reposition'
        
        if isempty(varargin), tmp = findcent;
        else tmp = varargin{1}; end
        if numel(tmp) == 3
            h = valid_handles(st.snap);
            if ~isempty(h)
                tmp = st.vols{h(1)}.mat * ...
                    round(st.vols{h(1)}.mat\[tmp(:); 1]);
            end
            st.centre = tmp(1:3);
        end
        redraw_all;
        callback;
        if isfield(st,'registry')
            xyz = bspm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
        end
        cm_pos;
        xyzstr = num2str([-99; round(xyz)]); 
        xyzstr(1,:) = [];
        axidx = [3 2 1];
        setvoxelinfo; 
        for a = 1:length(axidx)
            yh = st.vols{1}.ax{axidx(a)}.xyz;
            set(yh, 'string', xyzstr(a,:));
        end
           
    case 'setcoords'
        st.centre = varargin{1};
        st.centre = st.centre(:);
        redraw_all;
        callback;
        cm_pos;
        
    case 'space'
        if numel(varargin) < 1
            st.Space = eye(4);
            st.bb = maxbb;
            resolution;
            bbox;
            redraw_all;
        else
            space(varargin{:});
            resolution;
            bbox;
            redraw_all;
        end
        
    case 'maxbb'
        st.bb = maxbb;
        bbox;
        redraw_all;
        
    case 'resolution'
        resolution(varargin{:});
        bbox;
        redraw_all;
        
    case 'window'
        if numel(varargin)<2
            win = 'auto';
        elseif numel(varargin{2})==2
            win = varargin{2};
        end
        for i=valid_handles(varargin{1})
            st.vols{i}.window = win;
        end
        redraw(varargin{1});
        
    case 'delete'
        my_delete(varargin{1});
        
    case 'move'
        move(varargin{1},varargin{2});
        % redraw_all;
        
    case 'reset'
        my_reset;
        
    case 'pos'
        if isempty(varargin)
            H = st.centre(:);
        else
            H = pos(varargin{1});
        end
        varargout{1} = H;
        
    case 'interp'
        st.hld = varargin{1};
        redraw_all;
        
    case 'xhairs'
        xhairs(varargin{:});
        
    case 'register'
        register(varargin{1});
        
    case 'addblobs'
        addblobs(varargin{:});
        % redraw(varargin{1});
        
    case 'setblobsmax'
        st.vols{varargin{1}}.blobs{varargin{2}}.max = varargin{3};
        bspm_orthviews('redraw')
        
    case 'addcolouredblobs'
        addcolouredblobs(varargin{:});
        % redraw(varargin{1});
        
    case 'addimage'
        addimage(varargin{1}, varargin{2});
        % redraw(varargin{1});
        
    case 'addcolouredimage'
        addcolouredimage(varargin{1}, varargin{2},varargin{3});
        % redraw(varargin{1});
        
    case 'addtruecolourimage'
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
        
    case 'addcolourbar'
        addcolourbar(varargin{1}, varargin{2});
        
    case {'removeblobs','rmblobs'}
        rmblobs(varargin{1});
        redraw(varargin{1});
        
    case 'addcontext'
        if nargin == 1
            handles = 1:max_img;
        else
            handles = varargin{1};
        end
        addcontexts(handles);
        
    case {'removecontext','rmcontext'}
        if nargin == 1
            handles = 1:max_img;
        else
            handles = varargin{1};
        end
        rmcontexts(handles);
        
    case 'context_menu'
        c_menu(varargin{:});
        
    case 'valid_handles'
        if nargin == 1
            handles = 1:max_img;
        else
            handles = varargin{1};
        end
        varargout{1} = valid_handles(handles);

    case 'zoom'
        zoom_op(varargin{:});
        
    case 'zoommenu'
        if isempty(zoomlist)
            zoomlist = [NaN 0 5    10  20 40 80 Inf];
            reslist  = [1   1 .125 .25 .5 .5 1  1  ];
        end
        if nargin >= 3
            if all(cellfun(@isnumeric,varargin(1:2))) && ...
                    numel(varargin{1})==numel(varargin{2})
                zoomlist = varargin{1}(:);
                reslist  = varargin{2}(:);
            else
                warning('bspm_orthviews:zoom',...
                        'Invalid zoom or resolution list.')
            end
        end
        if nargout > 0
            varargout{1} = zoomlist;
        end
        if nargout > 1
            varargout{2} = reslist;
        end
        
    otherwise
        addonaction = strcmpi(st.plugins,action);
        if any(addonaction)
            feval(['spm_ov_' st.plugins{addonaction}],varargin{:});
        end
end

spm('Pointer','Arrow');
function H  = specify_image(img)
global st
H = [];
if isstruct(img)
    V = img(1);
else
    try
        V = spm_vol(img);
    catch
        fprintf('Can not use image "%s"\n', img);
        return;
    end
end
if numel(V)>1, V=V(1); end

ii = 1;
while ~isempty(st.vols{ii}), ii = ii + 1; end
DeleteFcn = ['bspm_orthviews(''Delete'',' num2str(ii) ');'];
V.ax = cell(3,1);
for i=1:3
    ax = axes('Visible','off', 'Parent', st.figax, ...
        'YDir','normal', 'DeleteFcn',DeleteFcn, 'ButtonDownFcn',@repos_start);
    d  = image(0, 'Tag','Transverse', 'Parent',ax, 'DeleteFcn',DeleteFcn);
    set(ax, 'Ydir','normal', 'ButtonDownFcn', @repos_start);
    lx = line(0,0, 'Parent',ax, 'DeleteFcn',DeleteFcn, 'Color',[0 0 1]);
    ly = line(0,0, 'Parent',ax, 'DeleteFcn',DeleteFcn, 'Color',[0 0 1]);
    if ~st.xhairs
        set(lx, 'Visible','off');
        set(ly, 'Visible','off');
    end
    V.ax{i} = struct('ax',ax,'d',d,'lx',lx,'ly',ly);
end
V.premul    = eye(4);
V.window    = 'auto';
V.mapping   = 'linear';
st.vols{ii} = V;
H = ii;
function bb = maxbb
global st
mn = [Inf Inf Inf];
mx = -mn;
for i=valid_handles
    premul = st.Space \ st.vols{i}.premul;
    bb = spm_get_bbox(st.vols{i}, 'fv', premul);
    mx = max([bb ; mx]);
    mn = min([bb ; mn]);
end
bb = [mn ; mx];
function H  = pos(handle)
global st
H = [];
for i=valid_handles(handle)
    is = inv(st.vols{i}.premul*st.vols{i}.mat);
    H = is(1:3,1:3)*st.centre(:) + is(1:3,4);
end
function mx = maxval(vol)
if isstruct(vol)
    mx = -Inf;
    for i=1:vol.dim(3)
        tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
        imx = max(tmp(isfinite(tmp)));
        if ~isempty(imx), mx = max(mx,imx); end
    end
else
    mx = max(vol(isfinite(vol)));
end
function mn = minval(vol)
if isstruct(vol)
    mn = Inf;
    for i=1:vol.dim(3)
        tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
        imn = min(tmp(isfinite(tmp)));
        if ~isempty(imn), mn = min(mn,imn); end
    end
else
    mn = min(vol(isfinite(vol)));
end
function m  = max_img
m = 24;
function centre = findcent
global st
obj    = get(st.fig,'CurrentObject');
centre = [];
cent   = [];
cp     = [];
for i=valid_handles
    for j=1:3
        if ~isempty(obj)
            if (st.vols{i}.ax{j}.ax == obj),
                cp = get(obj,'CurrentPoint');
            end
        end
        if ~isempty(cp)
            cp   = cp(1,1:2);
            is   = inv(st.Space);
            cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);
            switch j
                case 1
                    cent([1 2])=[cp(1)+st.bb(1,1)-1 cp(2)+st.bb(1,2)-1];
                case 2
                    cent([1 3])=[cp(1)+st.bb(1,1)-1 cp(2)+st.bb(1,3)-1];
                case 3
                    if st.mode ==0
                        cent([3 2])=[cp(1)+st.bb(1,3)-1 cp(2)+st.bb(1,2)-1];
                    else
                        cent([2 3])=[st.bb(2,2)+1-cp(1) cp(2)+st.bb(1,3)-1];
                    end
            end
            break;
        end
    end
    if ~isempty(cent), break; end
end
if ~isempty(cent), centre = st.Space(1:3,1:3)*cent(:) + st.Space(1:3,4); end
function handles = valid_handles(handles)
global st
if ~nargin, handles = 1:max_img; end
if isempty(st) || ~isfield(st,'vols')
    handles = [];
else
    handles = handles(:)';
    handles = handles(handles<=max_img & handles>=1 & ~rem(handles,1));
    for h=handles
        if isempty(st.vols{h}), handles(handles==h)=[]; end
    end
end
function img = scaletocmap(inpimg,mn,mx,cmap,miscol)
if nargin < 5, miscol=1; end
cml = size(cmap,1);
scf = (cml-1)/(mx-mn);
img = round((inpimg-mn)*scf)+1;
img(img<1)   = 1;
img(img>cml) = cml;
img(~isfinite(img)) = miscol;
function item_parent = addcontext(volhandle)
global st
% create context menu
set(0,'CurrentFigure',st.fig);
% contextmenu
item_parent = uicontextmenu;

% contextsubmenu 0
item00 = uimenu(item_parent, 'Label','unknown image', 'UserData','filename');
bspm_orthviews('context_menu','image_info',item00,volhandle);
item0a = uimenu(item_parent, 'UserData','pos_mm', 'Separator','on', ...
    'Callback','bspm_orthviews(''context_menu'',''repos_mm'');');
item0b = uimenu(item_parent, 'UserData','pos_vx', ...
    'Callback','bspm_orthviews(''context_menu'',''repos_vx'');');
item0c = uimenu(item_parent, 'UserData','v_value');

% contextsubmenu 1
item1    = uimenu(item_parent,'Label','Zoom', 'Separator','on');
[zl, rl] = bspm_orthviews('ZoomMenu');
for cz = numel(zl):-1:1
    if isinf(zl(cz))
        czlabel = 'Full Volume';
    elseif isnan(zl(cz))
        czlabel = 'BBox, this image > ...';
    elseif zl(cz) == 0
        czlabel = 'BBox, this image nonzero';
    else
        czlabel = sprintf('%dx%d mm', 2*zl(cz), 2*zl(cz));
    end
    item1_x = uimenu(item1, 'Label',czlabel,...
        'Callback', sprintf(...
        'bspm_orthviews(''context_menu'',''zoom'',%d,%d)',zl(cz),rl(cz)));
    if isinf(zl(cz)) % default display is Full Volume
        set(item1_x, 'Checked','on');
    end
end

% contextsubmenu 2
checked   = {'off','off'};
checked{st.xhairs+1} = 'on';
item2     = uimenu(item_parent,'Label','Crosshairs','Callback','bspm_orthviews(''context_menu'',''Xhair'');','Checked',checked{2});

% contextsubmenu 3
if st.Space == eye(4)
    checked = {'off', 'on'};
else
    checked = {'on', 'off'};
end
item3     = uimenu(item_parent,'Label','Orientation');
item3_1   = uimenu(item3,      'Label','World space', 'Callback','bspm_orthviews(''context_menu'',''orientation'',3);','Checked',checked{2});
item3_2   = uimenu(item3,      'Label','Voxel space (1st image)', 'Callback','bspm_orthviews(''context_menu'',''orientation'',2);','Checked',checked{1});
item3_3   = uimenu(item3,      'Label','Voxel space (this image)', 'Callback','bspm_orthviews(''context_menu'',''orientation'',1);','Checked','off');

% contextsubmenu 3
if isempty(st.snap)
    checked = {'off', 'on'};
else
    checked = {'on', 'off'};
end
item3     = uimenu(item_parent,'Label','Snap to Grid');
item3_1   = uimenu(item3,      'Label','Don''t snap', 'Callback','bspm_orthviews(''context_menu'',''snap'',3);','Checked',checked{2});
item3_2   = uimenu(item3,      'Label','Snap to 1st image', 'Callback','bspm_orthviews(''context_menu'',''snap'',2);','Checked',checked{1});
item3_3   = uimenu(item3,      'Label','Snap to this image', 'Callback','bspm_orthviews(''context_menu'',''snap'',1);','Checked','off');

% contextsubmenu 4
if st.hld == 0
    checked = {'off', 'off', 'on'};
elseif st.hld > 0
    checked = {'off', 'on', 'off'};
else
    checked = {'on', 'off', 'off'};
end
item4     = uimenu(item_parent,'Label','Interpolation');
item4_1   = uimenu(item4,      'Label','NN',    'Callback','bspm_orthviews(''context_menu'',''interpolation'',3);', 'Checked',checked{3});
item4_2   = uimenu(item4,      'Label','Trilin', 'Callback','bspm_orthviews(''context_menu'',''interpolation'',2);','Checked',checked{2});
item4_3   = uimenu(item4,      'Label','Sinc',  'Callback','bspm_orthviews(''context_menu'',''interpolation'',1);','Checked',checked{1});

% contextsubmenu 5
% item5     = uimenu(item_parent,'Label','Position', 'Callback','bspm_orthviews(''context_menu'',''position'');');

% contextsubmenu 6
item6       = uimenu(item_parent,'Label','Image','Separator','on');
item6_1     = uimenu(item6,      'Label','Window');
item6_1_1   = uimenu(item6_1,    'Label','local');
item6_1_1_1 = uimenu(item6_1_1,  'Label','auto', 'Callback','bspm_orthviews(''context_menu'',''window'',2);');
item6_1_1_2 = uimenu(item6_1_1,  'Label','manual', 'Callback','bspm_orthviews(''context_menu'',''window'',1);');
item6_1_1_3 = uimenu(item6_1_1,  'Label','percentiles', 'Callback','bspm_orthviews(''context_menu'',''window'',3);');
item6_1_2   = uimenu(item6_1,    'Label','global');
item6_1_2_1 = uimenu(item6_1_2,  'Label','auto', 'Callback','bspm_orthviews(''context_menu'',''window_gl'',2);');
item6_1_2_2 = uimenu(item6_1_2,  'Label','manual', 'Callback','bspm_orthviews(''context_menu'',''window_gl'',1);');
if license('test','image_toolbox') == 1
    offon = {'off', 'on'};
    checked = offon(strcmp(st.vols{volhandle}.mapping, ...
        {'linear', 'histeq', 'loghisteq', 'quadhisteq'})+1);
    item6_2     = uimenu(item6,      'Label','Intensity mapping');
    item6_2_1   = uimenu(item6_2,    'Label','local');
    item6_2_1_1 = uimenu(item6_2_1,  'Label','Linear', 'Checked',checked{1}, ...
        'Callback','bspm_orthviews(''context_menu'',''mapping'',''linear'');');
    item6_2_1_2 = uimenu(item6_2_1,  'Label','Equalised histogram', 'Checked',checked{2}, ...
        'Callback','bspm_orthviews(''context_menu'',''mapping'',''histeq'');');
    item6_2_1_3 = uimenu(item6_2_1,  'Label','Equalised log-histogram', 'Checked',checked{3}, ...
        'Callback','bspm_orthviews(''context_menu'',''mapping'',''loghisteq'');');
    item6_2_1_4 = uimenu(item6_2_1,  'Label','Equalised squared-histogram', 'Checked',checked{4}, ...
        'Callback','bspm_orthviews(''context_menu'',''mapping'',''quadhisteq'');');
    item6_2_2   = uimenu(item6_2,    'Label','global');
    item6_2_2_1 = uimenu(item6_2_2,  'Label','Linear', 'Checked',checked{1}, ...
        'Callback','bspm_orthviews(''context_menu'',''mapping_gl'',''linear'');');
    item6_2_2_2 = uimenu(item6_2_2,  'Label','Equalised histogram', 'Checked',checked{2}, ...
        'Callback','bspm_orthviews(''context_menu'',''mapping_gl'',''histeq'');');
    item6_2_2_3 = uimenu(item6_2_2,  'Label','Equalised log-histogram', 'Checked',checked{3}, ...
        'Callback','bspm_orthviews(''context_menu'',''mapping_gl'',''loghisteq'');');
    item6_2_2_4 = uimenu(item6_2_2,  'Label','Equalised squared-histogram', 'Checked',checked{4}, ...
        'Callback','bspm_orthviews(''context_menu'',''mapping_gl'',''quadhisteq'');');
end

% contextsubmenu 7
item7     = uimenu(item_parent,'Label','Overlay');
item7_1   = uimenu(item7,      'Label','Add blobs');
item7_1_1 = uimenu(item7_1,    'Label','local',  'Callback','bspm_orthviews(''context_menu'',''add_blobs'',2);');
item7_1_2 = uimenu(item7_1,    'Label','global', 'Callback','bspm_orthviews(''context_menu'',''add_blobs'',1);');
item7_2   = uimenu(item7,      'Label','Add image');
item7_2_1 = uimenu(item7_2,    'Label','local',  'Callback','bspm_orthviews(''context_menu'',''add_image'',2);');
item7_2_2 = uimenu(item7_2,    'Label','global', 'Callback','bspm_orthviews(''context_menu'',''add_image'',1);');
item7_3   = uimenu(item7,      'Label','Add coloured blobs','Separator','on');
item7_3_1 = uimenu(item7_3,    'Label','local',  'Callback','bspm_orthviews(''context_menu'',''add_c_blobs'',2);');
item7_3_2 = uimenu(item7_3,    'Label','global', 'Callback','bspm_orthviews(''context_menu'',''add_c_blobs'',1);');
item7_4   = uimenu(item7,      'Label','Add coloured image');
item7_4_1 = uimenu(item7_4,    'Label','local',  'Callback','bspm_orthviews(''context_menu'',''add_c_image'',2);');
item7_4_2 = uimenu(item7_4,    'Label','global', 'Callback','bspm_orthviews(''context_menu'',''add_c_image'',1);');
item7_5   = uimenu(item7,      'Label','Remove blobs',        'Visible','off','Separator','on');
item7_6   = uimenu(item7,      'Label','Remove coloured blobs','Visible','off');
item7_6_1 = uimenu(item7_6,    'Label','local', 'Visible','on');
item7_6_2 = uimenu(item7_6,    'Label','global','Visible','on');
item7_7   = uimenu(item7,      'Label','Set blobs max', 'Visible','off');

for i=1:3
    set(st.vols{volhandle}.ax{i}.ax,'UIcontextmenu',item_parent);
    st.vols{volhandle}.ax{i}.cm = item_parent;
end

% process any plugins
for k = 1:numel(st.plugins)
    feval(['spm_ov_', st.plugins{k}],'context_menu',volhandle,item_parent);
    if k==1
        h = get(item_parent,'Children');
        set(h(1),'Separator','on'); 
    end
end
function cm_handles = get_cm_handles
global st
cm_handles = [];
for i = valid_handles
    cm_handles = [cm_handles st.vols{i}.ax{1}.cm];
end
function addblobs(handle, xyz, t, mat, name)
global st
if nargin < 5
    name = '';
end
for i=valid_handles(handle)
    if ~isempty(xyz)
        rcp      = round(xyz);
        dim      = max(rcp,[],2)';
        off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
        vol      = zeros(dim)+NaN;
        vol(off) = t;
        vol      = reshape(vol,dim);
        st.vols{i}.blobs=cell(1,1);
        mx = max([eps max(t)]);
        mn = min([0 min(t)]);
        st.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'max',mx, 'min',mn,'name',name);
        addcolourbar(handle,1);
    end
end
function addimage(handle, fname)
global st
for i=valid_handles(handle)
    if isstruct(fname)
        vol = fname(1);
    else
        vol = spm_vol(fname);
    end
    mat = vol.mat;
    st.vols{i}.blobs=cell(1,1);
    mx = max([eps maxval(vol)]);
    mn = min([0 minval(vol)]);
    st.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'max',mx,'min',mn);
    addcolourbar(handle,1);
end
function addcolouredblobs(handle, xyz, t, mat, colour, name)
if nargin < 6
    name = '';
end
global st
for i=valid_handles(handle)
    if ~isempty(xyz)
        rcp      = round(xyz);
        dim      = max(rcp,[],2)';
        off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
        vol      = zeros(dim)+NaN;
        vol(off) = t;
        vol      = reshape(vol,dim);
        if ~isfield(st.vols{i},'blobs')
            st.vols{i}.blobs=cell(1,1);
            bset = 1;
        else
            bset = numel(st.vols{i}.blobs)+1;
        end
        mx = max([eps maxval(vol)]);
        mn = min([0 minval(vol)]);
        st.vols{i}.blobs{bset} = struct('vol',vol, 'mat',mat, ...
            'max',mx, 'min',mn, 'colour',colour, 'name',name);
    end
end
function addcolouredimage(handle, fname,colour)
global st
for i=valid_handles(handle)
    if isstruct(fname)
        vol = fname(1);
    else
        vol = spm_vol(fname);
    end
    mat = vol.mat;
    if ~isfield(st.vols{i},'blobs')
        st.vols{i}.blobs=cell(1,1);
        bset = 1;
    else
        bset = numel(st.vols{i}.blobs)+1;
    end
    mx = max([eps maxval(vol)]);
    mn = min([0 minval(vol)]);
    st.vols{i}.blobs{bset} = struct('vol',vol, 'mat',mat, ...
        'max',mx, 'min',mn, 'colour',colour);
end
function addtruecolourimage(handle,fname,colourmap,prop,mx,mn)
% adds true colour image to current displayed image
global st
for i=valid_handles(handle)
    if isstruct(fname)
        vol = fname(1);
    else
        vol = spm_vol(fname);
    end
    mat = vol.mat;
    if ~isfield(st.vols{i},'blobs')
        st.vols{i}.blobs=cell(1,1);
        bset = 1;
    else
        bset = numel(st.vols{i}.blobs)+1;
    end
    c = struct('cmap', colourmap,'prop',prop);
    st.vols{i}.blobs{bset} = struct('vol',vol, 'mat',mat, ...
        'max',mx, 'min',mn, 'colour',c);
    addcolourbar(handle,bset);
end
function addcontexts(handles)
for ii = valid_handles(handles)
    addcontext(ii);
end
bspm_orthviews('reposition',bspm_orthviews('pos'));
function rmblobs(handle)
global st
for i=valid_handles(handle)
    if isfield(st.vols{i},'blobs')
        for j=1:numel(st.vols{i}.blobs)
            if isfield(st.vols{i}.blobs{j},'cbar') && ishandle(st.vols{i}.blobs{j}.cbar),
                delete(st.vols{i}.blobs{j}.cbar);
            end
        end
        st.vols{i} = rmfield(st.vols{i},'blobs');
    end
end
function rmcontexts(handles)
global st
for ii = valid_handles(handles)
    for i=1:3
        set(st.vols{ii}.ax{i}.ax,'UIcontextmenu',[]);
        try, st.vols{ii}.ax{i} = rmfield(st.vols{ii}.ax{i},'cm'); end
    end
end
function register(hreg)
global st
%tmp = uicontrol('Position',[0 0 1 1],'Visible','off','Parent',st.fig);
h   = valid_handles;
if ~isempty(h)
    tmp = st.vols{h(1)}.ax{1}.ax;
    st.registry = struct('hReg',hreg,'hMe', tmp);
    bspm_XYZreg('Add2Reg',st.registry.hReg,st.registry.hMe, 'bspm_orthviews');
else
    warning('Nothing to register with');
end
st.centre = bspm_XYZreg('GetCoords',st.registry.hReg);
st.centre = st.centre(:);
function callback
global st
if ~iscell(st.callback), st.callback = { st.callback }; end
for i=1:numel(st.callback)
    if isa(st.callback{i},'function_handle')
        feval(st.callback{i});
    else
        eval(st.callback{i});
    end
end
function xhairs(state)
global st
if ~nargin, if st.xhairs, state = 'off'; else state = 'on'; end; end
st.xhairs = 0;
opt = 'on';
if ~strcmpi(state,'on')
    opt = 'off';
else
    st.xhairs = 1;
end
for i=valid_handles
    for j=1:3
        set(st.vols{i}.ax{j}.lx,'Visible',opt);
        set(st.vols{i}.ax{j}.ly,'Visible',opt);
    end
end
function cm_pos
global st
for i = 1:numel(valid_handles)
    if isfield(st.vols{i}.ax{1},'cm')
        set(findobj(st.vols{i}.ax{1}.cm,'UserData','pos_mm'),...
            'Label',sprintf('mm:  %.1f %.1f %.1f',bspm_orthviews('pos')));
        pos = bspm_orthviews('pos',i);
        set(findobj(st.vols{i}.ax{1}.cm,'UserData','pos_vx'),...
            'Label',sprintf('vx:  %.1f %.1f %.1f',pos));
        try
            Y = spm_sample_vol(st.vols{i},pos(1),pos(2),pos(3),st.hld);
        catch
            Y = NaN;
            fprintf('Cannot access file "%s".\n', st.vols{i}.fname);
        end
        set(findobj(st.vols{i}.ax{1}.cm,'UserData','v_value'),...
            'Label',sprintf('Y = %g',Y));
    end
end
function my_reset
global st
if ~isempty(st) && isfield(st,'registry') && ishandle(st.registry.hMe)
    delete(st.registry.hMe); st = rmfield(st,'registry');
end
my_delete(1:max_img);
reset_st;
function my_delete(handle)
global st
% remove blobs (and colourbars, if any)
rmblobs(handle);
% remove displayed axes
for i=valid_handles(handle)
    kids = get(st.fig,'Children');
    for j=1:3
        try
            if any(kids == st.vols{i}.ax{j}.ax)
                set(get(st.vols{i}.ax{j}.ax,'Children'),'DeleteFcn','');
                delete(st.vols{i}.ax{j}.ax);
            end
        end
    end
    st.vols{i} = [];
end
function resolution(res)
global st
if ~nargin, res = 1; end % Default minimum resolution 1mm
for i=valid_handles
    % adapt resolution to smallest voxel size of displayed images
    res  = min([res,sqrt(sum((st.vols{i}.mat(1:3,1:3)).^2))]);
end
res      = res/mean(svd(st.Space(1:3,1:3)));
Mat      = diag([res res res 1]);
st.Space = st.Space*Mat;
st.bb    = st.bb/res;
function move(handle,pos)
global st
for i=valid_handles(handle)
    st.vols{i}.area = pos;
end
bbox;
function space(handle,M,dim)
global st
if ~isempty(st.vols{handle})
    if nargin < 2
        M = st.vols{handle}.mat;
        dim = st.vols{handle}.dim(1:3);
    end
    Mat   = st.vols{handle}.premul(1:3,1:3)*M(1:3,1:3);
    vox   = sqrt(sum(Mat.^2));
    if det(Mat(1:3,1:3))<0, vox(1) = -vox(1); end
    Mat   = diag([vox 1]);
    Space = (M)/Mat;
    bb    = [1 1 1; dim];
    bb    = [bb [1;1]];
    bb    = bb*Mat';
    bb    = bb(:,1:3);
    bb    = sort(bb);
    st.Space = Space;
    st.bb = bb;
end
function zoom_op(fov,res)
global st
if nargin < 1, fov = Inf; end
if nargin < 2, res = Inf; end

if isinf(fov)
    st.bb = maxbb;
elseif isnan(fov) || fov == 0
    current_handle = valid_handles;
    if numel(current_handle) > 1 % called from check reg context menu
        current_handle = get_current_handle;
    end
    if fov == 0
        % zoom to bounding box of current image ~= 0
        thr = 'nz';
    else
        % zoom to bounding box of current image > chosen threshold
        thr = spm_input('Threshold (Y > ...)', '+1', 'r', '0', 1);
    end
    premul = st.Space \ st.vols{current_handle}.premul;
    st.bb = spm_get_bbox(st.vols{current_handle}, thr, premul);
else
    vx    = sqrt(sum(st.Space(1:3,1:3).^2));
    vx    = vx.^(-1);
    pos   = bspm_orthviews('pos');
    pos   = st.Space\[pos ; 1];
    pos   = pos(1:3)';
    st.bb = [pos-fov*vx; pos+fov*vx];
end
resolution(res);
bbox;
redraw_all;
if isfield(st.vols{1},'sdip')
    spm_eeg_inv_vbecd_disp('RedrawDip');
end
function repos_start(varargin)
    % don't use right mouse button to start reposition
    if ~strcmpi(get(gcbf,'SelectionType'),'alt')
        set(gcbf,'windowbuttonmotionfcn',@repos_move, 'windowbuttonupfcn',@repos_end);
        bspm_orthviews('reposition');
    end
function repos_move(varargin)
    bspm_orthviews('reposition');
function repos_end(varargin)
    set(gcbf,'windowbuttonmotionfcn','', 'windowbuttonupfcn','');
function bbox
global st
Dims = diff(st.bb)'+1;

TD = Dims([1 2])';
CD = Dims([1 3])';
if st.mode == 0, SD = Dims([3 2])'; else SD = Dims([2 3])'; end

un    = get(st.fig,'Units');set(st.fig,'Units','Pixels');
sz    = get(st.fig,'Position');set(st.fig,'Units',un);
sz    = sz(3:4);
sz(2) = sz(2)-40;

for i=valid_handles
    area   = st.vols{i}.area(:);
    area   = [area(1)*sz(1) area(2)*sz(2) area(3)*sz(1) area(4)*sz(2)];
    if st.mode == 0
        sx = area(3)/(Dims(1)+Dims(3))/1.02;
    else
        sx = area(3)/(Dims(1)+Dims(2))/1.02;
    end
    sy     = area(4)/(Dims(2)+Dims(3))/1.02;
    s      = min([sx sy]);
    
    offy   = (area(4)-(Dims(2)+Dims(3))*1.02*s)/2 + area(2);
    sky    = s*(Dims(2)+Dims(3))*0.02;
    if st.mode == 0
        offx = (area(3)-(Dims(1)+Dims(3))*1.02*s)/2 + area(1);
        skx  = s*(Dims(1)+Dims(3))*0.02;
    else
        offx = (area(3)-(Dims(1)+Dims(2))*1.02*s)/2 + area(1);
        skx  = s*(Dims(1)+Dims(2))*0.02;
    end
    
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
    if st.mode == 0
        set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
            'Position',[offx+s*Dims(1)+skx offy s*Dims(3) s*Dims(2)],...
            'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
            'Visible','on','XTick',[],'YTick',[]);
    else
        set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
            'Position',[offx+s*Dims(1)+skx offy+s*Dims(2)+sky s*Dims(2) s*Dims(3)],...
            'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
            'Visible','on','XTick',[],'YTick',[]);
    end
end
function redraw(arg1)
global st
bb   = st.bb;
Dims = round(diff(bb)'+1);
is   = inv(st.Space);
cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);

for i = valid_handles(arg1)
    M = st.Space\st.vols{i}.premul*st.vols{i}.mat;
    TM0 = [ 1 0 0 -bb(1,1)+1
            0 1 0 -bb(1,2)+1
            0 0 1 -cent(3)
            0 0 0 1];
    TM = inv(TM0*M);
    TD = Dims([1 2]);
    
    CM0 = [ 1 0 0 -bb(1,1)+1
            0 0 1 -bb(1,3)+1
            0 1 0 -cent(2)
            0 0 0 1];
    CM = inv(CM0*M);
    CD = Dims([1 3]);
    
    if st.mode ==0
        SM0 = [ 0 0 1 -bb(1,3)+1
                0 1 0 -bb(1,2)+1
                1 0 0 -cent(1)
                0 0 0 1];
        SM = inv(SM0*M); 
        SD = Dims([3 2]);
    else
        SM0 = [ 0 -1 0 +bb(2,2)+1
                0  0 1 -bb(1,3)+1
                1  0 0 -cent(1)
                0  0 0 1];
        SM = inv(SM0*M);
        SD = Dims([2 3]);
    end
    
    try
        imgt = spm_slice_vol(st.vols{i},TM,TD,st.hld)';
        imgc = spm_slice_vol(st.vols{i},CM,CD,st.hld)';
        imgs = spm_slice_vol(st.vols{i},SM,SD,st.hld)';
        imgc2 = adjustbrightness(imgc); 
        imgs2 = adjustbrightness(imgs); 
        imgt2 = adjustbrightness(imgt); 
        ok   = true;
    catch
        fprintf('Cannot access file "%s".\n', st.vols{i}.fname);
        fprintf('%s\n',getfield(lasterror,'message'));
        ok   = false;
    end
    
    if ok
        % get min/max threshold
        if strcmp(st.vols{i}.window,'auto')
            mn = -Inf;
            mx = Inf;
        else
            mn = min(st.vols{i}.window);
            mx = max(st.vols{i}.window);
        end
        % threshold images
        imgt = max(imgt,mn); imgt = min(imgt,mx);
        imgc = max(imgc,mn); imgc = min(imgc,mx);
        imgs = max(imgs,mn); imgs = min(imgs,mx);
        % compute intensity mapping, if histeq is available
        if license('test','image_toolbox') == 0
            st.vols{i}.mapping = 'linear';
        end
        switch st.vols{i}.mapping
            case 'linear'
            case 'histeq'
                % scale images to a range between 0 and 1
                imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                imgc = reshape(img(numel(imgt1)+(1:numel(imgc1))),size(imgc1));
                imgs = reshape(img(numel(imgt1)+numel(imgc1)+(1:numel(imgs1))),size(imgs1));
                mn = 0;
                mx = 1;
            case 'quadhisteq'
                % scale images to a range between 0 and 1
                imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                img  = histeq([imgt1(:).^2; imgc1(:).^2; imgs1(:).^2],1024);
                imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                imgc = reshape(img(numel(imgt1)+(1:numel(imgc1))),size(imgc1));
                imgs = reshape(img(numel(imgt1)+numel(imgc1)+(1:numel(imgs1))),size(imgs1));
                mn = 0;
                mx = 1;
            case 'loghisteq'
                sw = warning('off','MATLAB:log:logOfZero');
                imgt = log(imgt-min(imgt(:)));
                imgc = log(imgc-min(imgc(:)));
                imgs = log(imgs-min(imgs(:)));
                warning(sw);
                imgt(~isfinite(imgt)) = 0;
                imgc(~isfinite(imgc)) = 0;
                imgs(~isfinite(imgs)) = 0;
                % scale log images to a range between 0 and 1
                imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                imgc = reshape(img(numel(imgt1)+(1:numel(imgc1))),size(imgc1));
                imgs = reshape(img(numel(imgt1)+numel(imgc1)+(1:numel(imgs1))),size(imgs1));
                mn = 0;
                mx = 1;
        end
        % recompute min/max for display
        if strcmp(st.vols{i}.window,'auto')
            mx = -inf; mn = inf;
        end
        if ~isempty(imgt)
            tmp = imgt(isfinite(imgt));
            mx = max([mx max(max(tmp))]);
            mn = min([mn min(min(tmp))]);
        end
        if ~isempty(imgc)
            tmp = imgc(isfinite(imgc));
            mx = max([mx max(max(tmp))]);
            mn = min([mn min(min(tmp))]);
        end
        if ~isempty(imgs)
            tmp = imgs(isfinite(imgs));
            mx = max([mx max(max(tmp))]);
            mn = min([mn min(min(tmp))]);
        end
        if mx==mn, mx=mn+eps; end
        
        if isfield(st.vols{i},'blobs')
            if ~isfield(st.vols{i}.blobs{1},'colour')
                % Add blobs for display using the split colourmap
                scal = 64/(mx-mn);
                dcoff = -mn*scal;
                imgt = imgt*scal+dcoff;
                imgc = imgc*scal+dcoff;
                imgs = imgs*scal+dcoff;
                
                if isfield(st.vols{i}.blobs{1},'max')
                    mx = st.vols{i}.blobs{1}.max;
                else
                    mx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
                    st.vols{i}.blobs{1}.max = mx;
                end
                if isfield(st.vols{i}.blobs{1},'min')
                    mn = st.vols{i}.blobs{1}.min;
                else
                    mn = min([0 minval(st.vols{i}.blobs{1}.vol)]);
                    st.vols{i}.blobs{1}.min = mn;
                end
                
                vol  = st.vols{i}.blobs{1}.vol;
                M    = st.Space\st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
                tmpt = spm_slice_vol(vol,inv(TM0*M),TD,[0 NaN])';
                tmpc = spm_slice_vol(vol,inv(CM0*M),CD,[0 NaN])';
                tmps = spm_slice_vol(vol,inv(SM0*M),SD,[0 NaN])';
                
                %tmpt_z = find(tmpt==0);tmpt(tmpt_z) = NaN;
                %tmpc_z = find(tmpc==0);tmpc(tmpc_z) = NaN;
                %tmps_z = find(tmps==0);tmps(tmps_z) = NaN;
                
                sc   = 64/(mx-mn);
                off  = 65.51-mn*sc;
                msk  = find(isfinite(tmpt)); imgt(msk) = off+tmpt(msk)*sc;
                msk  = find(isfinite(tmpc)); imgc(msk) = off+tmpc(msk)*sc;
                msk  = find(isfinite(tmps)); imgs(msk) = off+tmps(msk)*sc;
                

                cmap = get(st.fig,'Colormap');
                if size(cmap,1)~=128
                    setcolormap(jet(64));
%                     spm_figure('Colormap','gray-hot')
                end
                figure(st.fig)
                redraw_colourbar(i,1,[mn mx],(1:64)'+64);
            elseif isstruct(st.vols{i}.blobs{1}.colour)
                % Add blobs for display using a defined colourmap
                
                % colourmaps
                gryc = (0:63)'*ones(1,3)/63;
                
                % scale grayscale image, not isfinite -> black
                gimgt = scaletocmap(imgt,mn,mx,gryc,65);
                gimgc = scaletocmap(imgc,mn,mx,gryc,65);
                gimgs = scaletocmap(imgs,mn,mx,gryc,65);
                gryc  = [gryc; 0 0 0];
                cactp = 0;
                
                for j=1:numel(st.vols{i}.blobs)
                    % colourmaps
                    actc = st.vols{i}.blobs{j}.colour.cmap;
                    actp = st.vols{i}.blobs{j}.colour.prop;
                    
                    % get min/max for blob image
                    if isfield(st.vols{i}.blobs{j},'max')
                        cmx = st.vols{i}.blobs{j}.max;
                    else
                        cmx = max([eps maxval(st.vols{i}.blobs{j}.vol)]);
                    end
                    if isfield(st.vols{i}.blobs{j},'min')
                        cmn = st.vols{i}.blobs{j}.min;
                    else
                        cmn = -cmx;
                    end
                    
                    % get blob data
                    vol  = st.vols{i}.blobs{j}.vol;
                    M    = st.Space\st.vols{i}.premul*st.vols{i}.blobs{j}.mat;
                    tmpt = spm_slice_vol(vol,inv(TM0*M),TD,[0 NaN])';
                    tmpc = spm_slice_vol(vol,inv(CM0*M),CD,[0 NaN])';
                    tmps = spm_slice_vol(vol,inv(SM0*M),SD,[0 NaN])';
                    
                    % actimg scaled round 0, black NaNs
                    topc = size(actc,1)+1;
                    tmpt = scaletocmap(tmpt,cmn,cmx,actc,topc);
                    tmpc = scaletocmap(tmpc,cmn,cmx,actc,topc);
                    tmps = scaletocmap(tmps,cmn,cmx,actc,topc);
                    actc = [actc; 0 0 0];
                    
                    % combine gray and blob data to truecolour
                    if isnan(actp)
                        if j==1, imgt = gryc(gimgt(:),:); end
                        imgt(tmpt~=size(actc,1),:) = actc(tmpt(tmpt~=size(actc,1)),:);
                        if j==1, imgc = gryc(gimgc(:),:); end
                        imgc(tmpc~=size(actc,1),:) = actc(tmpc(tmpc~=size(actc,1)),:);
                        if j==1, imgs = gryc(gimgs(:),:); end
                        imgs(tmps~=size(actc,1),:) = actc(tmps(tmps~=size(actc,1)),:);
                    else
                        cactp = cactp + actp;
                        if j==1, imgt = actc(tmpt(:),:)*actp; else imgt = imgt + actc(tmpt(:),:)*actp; end
                        if j==numel(st.vols{i}.blobs), imgt = imgt + gryc(gimgt(:),:)*(1-cactp); end
                        if j==1, imgc = actc(tmpc(:),:)*actp; else imgc = imgc + actc(tmpc(:),:)*actp; end
                        if j==numel(st.vols{i}.blobs), imgc = imgc + gryc(gimgc(:),:)*(1-cactp); end
                        if j==1, imgs = actc(tmps(:),:)*actp; else imgs = imgs + actc(tmps(:),:)*actp; end
                        if j==numel(st.vols{i}.blobs), imgs = imgs + gryc(gimgs(:),:)*(1-cactp); end
                    end
                    if j==numel(st.vols{i}.blobs)
                        imgt = reshape(imgt,[size(gimgt) 3]);
                        imgc = reshape(imgc,[size(gimgc) 3]);
                        imgs = reshape(imgs,[size(gimgs) 3]);
                    end
                    
                     % colourbar
                    csz   = size(st.vols{i}.blobs{j}.colour.cmap);
                    cdata = reshape(st.vols{i}.blobs{j}.colour.cmap, [csz(1) 1 csz(2)]);
                    redraw_colourbar(i,j,[cmn cmx],cdata);
                end
                
            else
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
                
                colour = zeros(numel(st.vols{i}.blobs),3);
                for j=1:numel(st.vols{i}.blobs) % get colours of all images first
                    if isfield(st.vols{i}.blobs{j},'colour')
                        colour(j,:) = reshape(st.vols{i}.blobs{j}.colour, [1 3]);
                    else
                        colour(j,:) = [1 0 0];
                    end
                end
                %colour = colour/max(sum(colour));
                
                for j=1:numel(st.vols{i}.blobs)
                    if isfield(st.vols{i}.blobs{j},'max')
                        mx = st.vols{i}.blobs{j}.max;
                    else
                        mx = max([eps max(st.vols{i}.blobs{j}.vol(:))]);
                        st.vols{i}.blobs{j}.max = mx;
                    end
                    if isfield(st.vols{i}.blobs{j},'min')
                        mn = st.vols{i}.blobs{j}.min;
                    else
                        mn = min([0 min(st.vols{i}.blobs{j}.vol(:))]);
                        st.vols{i}.blobs{j}.min = mn;
                    end
                    
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
                    cdata=permute(shiftdim((1/64:1/64:1)'* ...
                        colour(j,:),-1),[2 1 3]);
                    redraw_colourbar(i,j,[mn mx],cdata);
                end
                
                imgt = repmat(1-wt,[1 1 3]).*imgt+cimgt;
                imgc = repmat(1-wc,[1 1 3]).*imgc+cimgc;
                imgs = repmat(1-ws,[1 1 3]).*imgs+cimgs;
                
                imgt(imgt<0)=0; imgt(imgt>1)=1;
                imgc(imgc<0)=0; imgc(imgc>1)=1;
                imgs(imgs<0)=0; imgs(imgs>1)=1;
            end
        else
            scal = 64/(mx-mn);
            dcoff = -mn*scal;
            imgt = imgt*scal+dcoff;
            imgc = imgc*scal+dcoff;
            imgs = imgs*scal+dcoff;
        end

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
        if st.mode ==0
            set(st.vols{i}.ax{3}.lx,'HitTest','off',...
                'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
            set(st.vols{i}.ax{3}.ly,'HitTest','off',...
                'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(cent(3)-bb(1,3)+1));
        else
            set(st.vols{i}.ax{3}.lx,'HitTest','off',...
                'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
            set(st.vols{i}.ax{3}.ly,'HitTest','off',...
                'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(bb(2,2)+1-cent(2)));
        end
        
        if ~isempty(st.plugins) % process any addons
            for k = 1:numel(st.plugins)
                if isfield(st.vols{i},st.plugins{k})
                    feval(['spm_ov_', st.plugins{k}], ...
                        'redraw', i, TM0, TD, CM0, CD, SM0, SD);
                end
            end
        end
    end
end
drawnow;
function redraw_all
redraw(1:max_img);
function reset_st
global st
fig = spm_figure('FindWin','Graphics');
bb  = []; %[ [-78 78]' [-112 76]' [-50 85]' ];
st  = struct('n', 0, 'vols',{cell(max_img,1)}, 'bb',bb, 'Space',eye(4), ...
             'centre',[0 0 0], 'callback',';', 'xhairs',1, 'hld',1, ...
             'fig',fig, 'mode',1, 'plugins',{{}}, 'snap',[]);

xTB = spm('TBs');
if ~isempty(xTB)
    pluginbase = {spm('Dir') xTB.dir};
else
    pluginbase = {spm('Dir')};
end
for k = 1:numel(pluginbase)
    pluginpath = fullfile(pluginbase{k},'bspm_orthviews');
    if isdir(pluginpath)
        pluginfiles = dir(fullfile(pluginpath,'spm_ov_*.m'));
        if ~isempty(pluginfiles)
            if ~isdeployed, addpath(pluginpath); end
            for l = 1:numel(pluginfiles)
                pluginname = spm_file(pluginfiles(l).name,'basename');
                st.plugins{end+1} = strrep(pluginname, 'spm_ov_','');
            end
        end
    end
end
function c_menu(varargin)
global st

switch lower(varargin{1})
    case 'image_info'
        if nargin <3
            current_handle = get_current_handle;
        else
            current_handle = varargin{3};
        end
        if isfield(st.vols{current_handle},'fname')
            [p,n,e,v] = spm_fileparts(st.vols{current_handle}.fname);
            if isfield(st.vols{current_handle},'n')
                v = sprintf(',%d',st.vols{current_handle}.n);
            end
            set(varargin{2}, 'Label',[n e v]);
        end
        delete(get(varargin{2},'children'));
        if exist('p','var')
            item1 = uimenu(varargin{2}, 'Label', p);
        end
        if isfield(st.vols{current_handle},'descrip')
            item2 = uimenu(varargin{2}, 'Label',...
                st.vols{current_handle}.descrip);
        end
        dt = st.vols{current_handle}.dt(1);
        item3 = uimenu(varargin{2}, 'Label', sprintf('Data type: %s', spm_type(dt)));
        str   = 'Intensity: varied';
        if size(st.vols{current_handle}.pinfo,2) == 1
            if st.vols{current_handle}.pinfo(2)
                str = sprintf('Intensity: Y = %g X + %g',...
                    st.vols{current_handle}.pinfo(1:2)');
            else
                str = sprintf('Intensity: Y = %g X', st.vols{current_handle}.pinfo(1)');
            end
        end
        item4  = uimenu(varargin{2}, 'Label',str);
        item5  = uimenu(varargin{2}, 'Label', 'Image dimensions', 'Separator','on');
        item51 = uimenu(varargin{2}, 'Label',...
            sprintf('%dx%dx%d', st.vols{current_handle}.dim(1:3)));
        
        prms   = spm_imatrix(st.vols{current_handle}.mat);
        item6  = uimenu(varargin{2}, 'Label', 'Voxel size', 'Separator','on');
        item61 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', prms(7:9)));
        
        O      = st.vols{current_handle}.mat\[0 0 0 1]'; O=O(1:3)';
        item7  = uimenu(varargin{2}, 'Label', 'Origin', 'Separator','on');
        item71 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', O));
        
        R      = spm_matrix([0 0 0 prms(4:6)]);
        item8  = uimenu(varargin{2}, 'Label', 'Rotations', 'Separator','on');
        item81 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(1,1:3)));
        item82 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(2,1:3)));
        item83 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(3,1:3)));
        item9  = uimenu(varargin{2},...
            'Label','Specify other image...',...
            'Callback','bspm_orthviews(''context_menu'',''swap_img'');',...
            'Separator','on');
        
    case 'repos_mm'
        oldpos_mm = bspm_orthviews('pos');
        newpos_mm = spm_input('New Position (mm)','+1','r',sprintf('%.2f %.2f %.2f',oldpos_mm),3);
        bspm_orthviews('reposition',newpos_mm);
        
    case 'repos_vx'
        current_handle = get_current_handle;
        oldpos_vx = bspm_orthviews('pos', current_handle);
        newpos_vx = spm_input('New Position (voxels)','+1','r',sprintf('%.2f %.2f %.2f',oldpos_vx),3);
        newpos_mm = st.vols{current_handle}.mat*[newpos_vx;1];
        bspm_orthviews('reposition',newpos_mm(1:3));
        
    case 'zoom'
        zoom_all(varargin{2:end});
        bbox;
        redraw_all;
        
    case 'xhair'
        bspm_orthviews('Xhairs',varargin{2:end});
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles)
            z_handle = findobj(cm_handles(i),'label','Crosshairs');
            if st.xhairs
                set(z_handle,'Checked','on');
            else
                set(z_handle,'Checked','off');
            end
        end
        
    case 'orientation'
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles)
            z_handle = get(findobj(cm_handles(i),'label','Orientation'),'Children');
            set(z_handle,'Checked','off');
        end
        if varargin{2} == 3
            bspm_orthviews('Space');
            for i = 1:numel(cm_handles),
                z_handle = findobj(cm_handles(i),'label','World space');
                set(z_handle,'Checked','on');
            end
        elseif varargin{2} == 2,
            bspm_orthviews('Space',1);
            for i = 1:numel(cm_handles)
                z_handle = findobj(cm_handles(i),'label',...
                    'Voxel space (1st image)');
                set(z_handle,'Checked','on');
            end
        else
            bspm_orthviews('Space',get_current_handle);
            z_handle = findobj(st.vols{get_current_handle}.ax{1}.cm, ...
                'label','Voxel space (this image)');
            set(z_handle,'Checked','on');
            return;
        end
        
    case 'snap'
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles)
            z_handle = get(findobj(cm_handles(i),'label','Snap to Grid'),'Children');
            set(z_handle,'Checked','off');
        end
        if varargin{2} == 3
            st.snap = [];
        elseif varargin{2} == 2
            st.snap = 1;
        else
            st.snap = get_current_handle;
            z_handle = get(findobj(st.vols{get_current_handle}.ax{1}.cm,'label','Snap to Grid'),'Children');
            set(z_handle(1),'Checked','on');
            return;
        end
        for i = 1:numel(cm_handles)
            z_handle = get(findobj(cm_handles(i),'label','Snap to Grid'),'Children');
            set(z_handle(varargin{2}),'Checked','on');
        end
        
    case 'interpolation'
        tmp        = [-4 1 0];
        st.hld     = tmp(varargin{2});
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles)
            z_handle = get(findobj(cm_handles(i),'label','Interpolation'),'Children');
            set(z_handle,'Checked','off');
            set(z_handle(varargin{2}),'Checked','on');
        end
        redraw_all;
        
    case 'window'
        current_handle = get_current_handle;
        if varargin{2} == 2
            bspm_orthviews('window',current_handle);
        elseif varargin{2} == 3
            pc = spm_input('Percentiles', '+1', 'w', '3 97', 2, 100);
            wn = spm_summarise(st.vols{current_handle}, 'all', ...
                @(X) spm_percentile(X, pc));
            bspm_orthviews('window',current_handle,wn);
        else
            if isnumeric(st.vols{current_handle}.window)
                defstr = sprintf('%.2f %.2f', st.vols{current_handle}.window);
            else
                defstr = '';
            end
            [w,yp] = spm_input('Range','+1','e',defstr,[1 inf]);
            while numel(w) < 1 || numel(w) > 2
                uiwait(warndlg('Window must be one or two numbers','Wrong input size','modal'));
                [w,yp] = spm_input('Range',yp,'e',defstr,[1 inf]);
            end
            if numel(w) == 1
                w(2) = w(1)+eps;
            end
            bspm_orthviews('window',current_handle,w);
        end
        
    case 'window_gl'
        if varargin{2} == 2
            for i = 1:numel(get_cm_handles)
                st.vols{i}.window = 'auto';
            end
        else
            current_handle = get_current_handle;
            if isnumeric(st.vols{current_handle}.window)
                defstr = sprintf('%d %d', st.vols{current_handle}.window);
            else
                defstr = '';
            end
            [w,yp] = spm_input('Range','+1','e',defstr,[1 inf]);
            while numel(w) < 1 || numel(w) > 2
                uiwait(warndlg('Window must be one or two numbers','Wrong input size','modal'));
                [w,yp] = spm_input('Range',yp,'e',defstr,[1 inf]);
            end
            if numel(w) == 1
                w(2) = w(1)+eps;
            end
            for i = 1:numel(get_cm_handles)
                st.vols{i}.window = w;
            end
        end
        redraw_all;
        
    case 'mapping'
        checked = strcmp(varargin{2}, ...
            {'linear', 'histeq', 'loghisteq', 'quadhisteq'});
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
        end
        redraw_all;
        
    case 'mapping_gl'
        checked = strcmp(varargin{2}, ...
            {'linear', 'histeq', 'loghisteq', 'quadhisteq'});
        checked = checked(end:-1:1); % Handles are stored in inverse order
        cm_handles = get_cm_handles;
        for k = valid_handles
            st.vols{k}.mapping = varargin{2};
            z_handle = get(findobj(cm_handles(k), ...
                'label','Intensity mapping'),'Children');
            for l = 1:numel(z_handle)
                c_handle = get(z_handle(l), 'Children');
                set(c_handle, 'checked', 'off');
                set(c_handle(checked), 'checked', 'on');
            end
        end
        redraw_all;
        
    case 'swap_img'
        current_handle = get_current_handle;
        newimg = spm_select(1,'image','select new image');
        if ~isempty(newimg)
            new_info = spm_vol(newimg);
            fn = fieldnames(new_info);
            for k=1:numel(fn)
                st.vols{current_handle}.(fn{k}) = new_info.(fn{k});
            end
            bspm_orthviews('context_menu','image_info',get(gcbo, 'parent'));
            redraw_all;
        end
        
    case 'add_blobs'
        % Add blobs to the image - in split colortable
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        spm_input('!DeleteInputObj');
        [SPM,xSPM] = spm_getSPM;
        if ~isempty(SPM)
            for i = 1:numel(cm_handles)
                addblobs(cm_handles(i),xSPM.XYZ,xSPM.Z,xSPM.M);
                % Add options for removing blobs
                c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Remove blobs');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                item7_3_1 = uimenu(c_handle,'Label','local','Callback','bspm_orthviews(''context_menu'',''remove_blobs'',2);');
                if varargin{2} == 1,
                    item7_3_2 = uimenu(c_handle,'Label','global','Callback','bspm_orthviews(''context_menu'',''remove_blobs'',1);');
                end
                % Add options for setting maxima for blobs
                c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Set blobs max');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                uimenu(c_handle,'Label','local','Callback','bspm_orthviews(''context_menu'',''setblobsmax'',2);');
                if varargin{2} == 1
                    uimenu(c_handle,'Label','global','Callback','bspm_orthviews(''context_menu'',''setblobsmax'',1);');
                end
            end
            redraw_all;
        end
        
    case 'remove_blobs'
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        for i = 1:numel(cm_handles)
            rmblobs(cm_handles(i));
            % Remove options for removing blobs
            c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Remove blobs');
            delete(get(c_handle,'Children'));
            set(c_handle,'Visible','off');
            % Remove options for setting maxima for blobs
            c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Set blobs max');
            set(c_handle,'Visible','off');
        end
        redraw_all;
        
    case 'add_image'
        % Add blobs to the image - in split colortable
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        spm_input('!DeleteInputObj');
        fname = spm_select(1,'image','select image');
        if ~isempty(fname)
            for i = 1:numel(cm_handles)
                addimage(cm_handles(i),fname);
                % Add options for removing blobs
                c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Remove blobs');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                item7_3_1 = uimenu(c_handle,'Label','local','Callback','bspm_orthviews(''context_menu'',''remove_blobs'',2);');
                if varargin{2} == 1
                    item7_3_2 = uimenu(c_handle,'Label','global','Callback','bspm_orthviews(''context_menu'',''remove_blobs'',1);');
                end
                % Add options for setting maxima for blobs
                c_handle = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Set blobs max');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                uimenu(c_handle,'Label','local','Callback','bspm_orthviews(''context_menu'',''setblobsmax'',2);');
                if varargin{2} == 1
                    uimenu(c_handle,'Label','global','Callback','bspm_orthviews(''context_menu'',''setblobsmax'',1);');
                end
            end
            redraw_all;
        end
        
    case 'add_c_blobs'
        % Add blobs to the image - in full colour
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        spm_input('!DeleteInputObj');
        [SPM,xSPM] = spm_getSPM;
        if ~isempty(SPM)
            c = spm_input('Colour','+1','m',...
                'Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
            colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
            c_names = {'red';'yellow';'green';'cyan';'blue';'magenta'};
            hlabel = sprintf('%s (%s)',xSPM.title,c_names{c});
            for i = 1:numel(cm_handles)
                addcolouredblobs(cm_handles(i),xSPM.XYZ,xSPM.Z,xSPM.M,colours(c,:),xSPM.title);
                addcolourbar(cm_handles(i),numel(st.vols{cm_handles(i)}.blobs));
                c_handle    = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Remove coloured blobs');
                ch_c_handle = get(c_handle,'Children');
                set(c_handle,'Visible','on');
                %set(ch_c_handle,'Visible',on');
                item7_4_1   = uimenu(ch_c_handle(2),'Label',hlabel,'ForegroundColor',colours(c,:),...
                    'Callback','c = get(gcbo,''UserData'');bspm_orthviews(''context_menu'',''remove_c_blobs'',2,c);',...
                    'UserData',c);
                if varargin{2} == 1
                    item7_4_2 = uimenu(ch_c_handle(1),'Label',hlabel,'ForegroundColor',colours(c,:),...
                        'Callback','c = get(gcbo,''UserData'');bspm_orthviews(''context_menu'',''remove_c_blobs'',1,c);',...
                        'UserData',c);
                end
            end
            redraw_all;
        end
        
    case 'remove_c_blobs'
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
        for i = 1:numel(cm_handles)
            if isfield(st.vols{cm_handles(i)},'blobs')
                for j = 1:numel(st.vols{cm_handles(i)}.blobs)
                    if all(st.vols{cm_handles(i)}.blobs{j}.colour == colours(varargin{3},:));
                        if isfield(st.vols{cm_handles(i)}.blobs{j},'cbar')
                            delete(st.vols{cm_handles(i)}.blobs{j}.cbar);
                        end
                        st.vols{cm_handles(i)}.blobs(j) = [];
                        break;
                    end
                end
                rm_c_menu = findobj(st.vols{cm_handles(i)}.ax{1}.cm,'Label','Remove coloured blobs');
                delete(gcbo);
                if isempty(st.vols{cm_handles(i)}.blobs)
                    st.vols{cm_handles(i)} = rmfield(st.vols{cm_handles(i)},'blobs');
                    set(rm_c_menu, 'Visible', 'off');
                end
            end
        end
        redraw_all;
        
    case 'add_c_image'
        % Add truecolored image
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        spm_input('!DeleteInputObj');
        fname = spm_select([1 Inf],'image','select image(s)');
        for k = 1:size(fname,1)
            c = spm_input(sprintf('Image %d: Colour',k),'+1','m',...
                'Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
            colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
            c_names = {'red';'yellow';'green';'cyan';'blue';'magenta'};
            hlabel = sprintf('%s (%s)',fname(k,:),c_names{c});
            for i = 1:numel(cm_handles)
                addcolouredimage(cm_handles(i),fname(k,:),colours(c,:));
                addcolourbar(cm_handles(i),numel(st.vols{cm_handles(i)}.blobs));
                c_handle    = findobj(findobj(st.vols{cm_handles(i)}.ax{1}.cm,'label','Overlay'),'Label','Remove coloured blobs');
                ch_c_handle = get(c_handle,'Children');
                set(c_handle,'Visible','on');
                %set(ch_c_handle,'Visible',on');
                item7_4_1 = uimenu(ch_c_handle(2),'Label',hlabel,'ForegroundColor',colours(c,:),...
                    'Callback','c = get(gcbo,''UserData'');bspm_orthviews(''context_menu'',''remove_c_blobs'',2,c);','UserData',c);
                if varargin{2} == 1
                    item7_4_2 = uimenu(ch_c_handle(1),'Label',hlabel,'ForegroundColor',colours(c,:),...
                        'Callback','c = get(gcbo,''UserData'');bspm_orthviews(''context_menu'',''remove_c_blobs'',1,c);',...
                        'UserData',c);
                end
            end
            redraw_all;
        end
        
    case 'setblobsmax'
        if varargin{2} == 1
            % global
            cm_handles = valid_handles;
            mx = -inf;
            for i = 1:numel(cm_handles)
                if ~isfield(st.vols{cm_handles(i)}, 'blobs'), continue, end
                for j = 1:numel(st.vols{cm_handles(i)}.blobs)
                    mx = max(mx, st.vols{cm_handles(i)}.blobs{j}.max);
                end
            end
            mx = spm_input('Maximum value', '+1', 'r', mx, 1);
            for i = 1:numel(cm_handles)
                if ~isfield(st.vols{cm_handles(i)}, 'blobs'), continue, end
                for j = 1:numel(st.vols{cm_handles(i)}.blobs)
                    st.vols{cm_handles(i)}.blobs{j}.max = mx;
                end
            end
        else
            % local (should handle coloured blobs, but not implemented yet)
            cm_handle = get_current_handle;
            colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
            if ~isfield(st.vols{cm_handle}, 'blobs'), return, end
            for j = 1:numel(st.vols{cm_handle}.blobs)
                if nargin < 4 || ...
                        all(st.vols{cm_handle}.blobs{j}.colour == colours(varargin{3},:))
                    mx = st.vols{cm_handle}.blobs{j}.max;
                    mx = spm_input('Maximum value', '+1', 'r', mx, 1);
                    st.vols{cm_handle}.blobs{j}.max = mx;
                end
            end
        end
        redraw_all;
end
function current_handle = get_current_handle
cm_handle      = get(gca,'UIContextMenu');
cm_handles     = get_cm_handles;
current_handle = find(cm_handles==cm_handle);
function zoom_all(zoom,res)
cm_handles = get_cm_handles;
zoom_op(zoom,res);
for i = 1:numel(cm_handles)
    z_handle = get(findobj(cm_handles(i),'label','Zoom'),'Children');
    set(z_handle,'Checked','off');
    if isinf(zoom)
        set(findobj(z_handle,'Label','Full Volume'),'Checked','on');
    elseif zoom > 0
        set(findobj(z_handle,'Label',sprintf('%dx%d mm', 2*zoom, 2*zoom)),'Checked','on');
    end % leave all unchecked if either bounding box option was chosen
end

% | bspm_XYZreg (MODIFIED FROM ORIGINAL SPM8 CODE)
% =========================================================================
function varargout=bspm_XYZreg(varargin)
% Registry for GUI XYZ locations, and point list utility functions
%
%                           ----------------
%
% PointList & voxel centre utilities...
%
% FORMAT [xyz,d] = bspm_XYZreg('RoundCoords',xyz,M,D)
% FORMAT [xyz,d] = bspm_XYZreg('RoundCoords',xyz,V)
% Rounds specified xyz location to nearest voxel centre
% xyz - (Input) 3-vector of X, Y & Z locations, in "real" co-ordinates
% M   - 4x4 transformation matrix relating voxel to "real" co-ordinates
% D   - 3 vector of image X, Y & Z dimensions (DIM)
% V   - 9-vector of image and voxel sizes, and origin [DIM,VOX,ORIGIN]'
%       M derived as [ [diag(V(4:6)), -(V(7:9).*V(4:6))]; [zeros(1,3) ,1]]
%       DIM    - D
%       VOX    - Voxel dimensions in units of "real" co-ordinates
%       ORIGIN - Origin of "real" co-ordinates in voxel co-ordinates
% xyz - (Output) co-ordinates of nearest voxel centre in "real" co-ordinates
% d   - Euclidean distance between requested xyz & nearest voxel centre
%
% FORMAT i = bspm_XYZreg('FindXYZ',xyz,XYZ)
% finds position of specified voxel in XYZ pointlist
% xyz - 3-vector of co-ordinates
% XYZ - Pointlist: 3xn matrix of co-ordinates
% i   - Column(s) of XYZ equal to xyz
%
% FORMAT [xyz,i,d] = bspm_XYZreg('NearestXYZ',xyz,XYZ)
% find nearest voxel in pointlist to specified location
% xyz - (Input) 3-vector of co-ordinates
% XYZ - Pointlist: 3xn matrix of co-ordinates
% xyz - (Output) co-ordinates of nearest voxel in XYZ pointlist
%       (ties are broken in favour of the first location in the pointlist)
% i   - Column of XYZ containing co-ordinates of nearest pointlist location
% d   - Euclidean distance between requested xyz & nearest pointlist location
%
% FORMAT d = bspm_XYZreg('Edist',xyz,XYZ)
% Euclidean distances between co-ordinates xyz & points in XYZ pointlist
% xyz - 3-vector of co-ordinates
% XYZ - Pointlist: 3xn matrix of co-ordinates
% d   - n row-vector of Euclidean distances between xyz & points of XYZ
%
%                           ----------------
% Registry functions
%
% FORMAT [hReg,xyz] = bspm_XYZreg('InitReg',hReg,M,D,xyz)
% Initialise registry in graphics object
% hReg - Handle of HandleGraphics object to build registry in. Object must
%        be un'Tag'ged and have empty 'UserData'
% M    - 4x4 transformation matrix relating voxel to "real" co-ordinates, used
%        and stored for checking validity of co-ordinates
% D    - 3 vector of image X, Y & Z dimensions (DIM), used
%        and stored for checking validity of co-ordinates
% xyz  - (Input) Initial co-ordinates [Default [0;0;0]]
%        These are rounded to the nearest voxel centre
% hReg - (Output) confirmation of registry handle
% xyz  - (Output) Current registry co-ordinates, after rounding
%
% FORMAT bspm_XYZreg('UnInitReg',hReg)
% Clear registry information from graphics object
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object.
%        Object's 'Tag' & 'UserData' are cleared
%
% FORMAT xyz = bspm_XYZreg('GetCoords',hReg)
% Get current registry co-ordinates
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
% 
% FORMAT [xyz,d] = bspm_XYZreg('SetCoords',xyz,hReg,hC,Reg)
% Set co-ordinates in registry & update registered HGobjects/functions
% xyz  - (Input) desired co-ordinates
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
%        If hReg doesn't contain a registry, a warning is printed.
% hC   - Handle of caller object (to prevent circularities) [Default 0]
%        If caller object passes invalid registry handle, then bspm_XYZreg
%        attempts to blank the 'hReg' fiend of hC's 'UserData', printing
%        a warning notification.
% Reg  - Alternative nx2 cell array Registry of handles / functions
%        If specified, overrides use of registry held in hReg
%        [Default getfield(get(hReg,'UserData'),'Reg')]
% xyz  - (Output) Desired co-ordinates are rounded to nearest voxel if hC
%        is not specified, or is zero. Otherwise, caller is assummed to
%        have checked verity of desired xyz co-ordinates. Output xyz returns
%        co-ordinates actually set.
% d    - Euclidean distance between desired and set co-ordinates.
%
% FORMAT nReg = bspm_XYZreg('XReg',hReg,{h,Fcn}pairs)
% Cross registration object/function pairs with the registry, push xyz co-ords
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
% h    - Handle of HandleGraphics object to be registered
%        The 'UserData' of h must be a structure with an 'Reg' field, which
%        is set to hReg, the handle of the registry (back registration)
% Fcn  - Handling function for HandleGraphics object h
%        This function *must* accept XYZ updates via the call:
%                feval(Fcn,'SetCoords',xyz,h,hReg)
%        and should *not* call back the registry with the update!
%        {h,Fcn} are appended to the registry (forward registration)
% nReg - New registry cell array: Handles are checked for validity before
%        entry. Invalid handles are omitted, generating a warning.
%
% FORMAT nReg = bspm_XYZreg('Add2Reg',hReg,{h,Fcn}pairs)
% Add object/function pairs for XYZ updates to registry (forward registration)
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
% h    - Handle of HandleGraphics object to be registered
% Fcn  - Handling function for HandleGraphics object h
%        This function *must* accept XYZ updates via the call:
%                feval(Fcn,'SetCoords',xyz,h,hReg)
%        and should *not* call back the registry with the update!
%        {h,Fcn} are appended to the registry (forward registration)
% nReg - New registry cell array: Handles are checked for validity before
%        entry. Invalid handles are omitted, generating a warning.
%
% FORMAT bspm_XYZreg('SetReg',h,hReg)
% Set registry field of object's UserData (back registration)
% h    - Handle of HandleGraphics object to be registered
%        The 'UserData' of h must be a structure with an 'Reg' field, which
%        is set to hReg, the handle of the registry (back registration)
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
%
% FORMAT nReg = bspm_XYZreg('unXReg',hReg,hD1,hD2,hD3,...)
% Un-cross registration of HandleGraphics object hD
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
% hD?  - Handles of HandleGraphics object to be unregistered
%        The 'UserData' of hD must be a structure with a 'Reg' field, which
%        is set to empty (back un-registration)
% nReg - New registry cell array: Registry entries with handle entry hD are 
%        removed from the registry (forward un-registration)
%        Handles not in the registry generate a warning
%
% FORMAT nReg = bspm_XYZreg('Del2Reg',hReg,hD)
% Delete HandleGraphics object hD from registry (forward un-registration)
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
% hD?  - Handles of HandleGraphics object to be unregistered
% nReg - New registry cell array: Registry entries with handle entry hD are 
%        removed from the registry. Handles not in registry generate a warning
%
% FORMAT bspm_XYZreg('UnSetReg',h)
% Unset registry field of object's UserData (back un-registration)
% h - Handle of HandleGraphics object to be unregistered
%     The 'UserData' of hD must be a structure with a 'Reg' field, which
%     is set to empty (back un-registration)
%
% FORMAT bspm_XYZreg('CleanReg',hReg)
% Clean invalid handles from registry
% hReg - Handle of 'hReg' 'Tag'ged registry HandleGraphics object
%
% FORMAT Reg = bspm_XYZreg('VReg',Reg,Warn)
% Prune invalid handles from Registry cell array
% Reg  - (Input) nx2 cell array of {handle,function} pairs
% Warn - If specified, print warning if find invalid handles
% Reg  - (Output) mx2 cell array of valid {handle,function} pairs
%
% FORMAT hReg = bspm_XYZreg('FindReg',h)
% Find/check registry object
% h    - handle of Registry, or figure containing Registry (default gcf)
%        If ischar(h), then uses spm_figure('FindWin',h) to locate named figures
% hReg - handle of confirmed registry object
%        Errors if h is not a registry or a figure containing a unique registry
%        Registry object is identified by 'hReg' 'Tag'
%_______________________________________________________________________
%
% bspm_XYZreg provides a framework for modular inter-GUI communication of
% XYZ co-orginates, and various utility functions for pointlist handling
% and rounding in voxel co-ordinates.
%
%-----------------------------------------------------------------------
%                                                           THE REGISTRY
%
% The concept of the registry is of a central entity which "knows"
% about other GUI objects holding XYZ co-ordinates, and keeps them all
% in sync. Changes to the registry's XYZ co-ordinates are passed on to
% registered functions by the registry (forward registration).
% Individual objects which can change the XYZ co-ordinates should
% therefore update the registry with the new co-ordinates (back
% registration), so that the registry can tell all registered objects
% about the new location, and a framework is provided for this.
%
% The registry is held as the 'UserData of a HandleGraphics object,
% whose handle therefore identifies the registry. The registry object
% is 'Tag'ged 'hReg' for identification (though this 'Tag' is not used
% for locating the registry, so multiple registry incarnations are
% possible). The registry object's 'UserData' is a structure containing
% the current XYZ co-ordinates, the voxel-to-co-ordinates matrix M, the
% image dimensions D, and the Registry itself. The registry is a nx2
% cell array containing n handle/function pairs.
%
% The model is that all GUI objects requiring linking to a common XYZ
% location via the registry each be identified by a HandleGraphics
% handle. This handle can be the handle of the particular instantiation
% of the GUI control itself (as is the case with the MIP-GUI of
% spm_mip_ui where the axis handle is used to identify the MIP to use);
% the handle of another HandleGraphics object associated with the GUI
% control (as is the case with the XYZ editable widgets of
% spm_results_ui where the handle of the bounding frame uicontrol is
% used); or may be 0, the handle of the root object, which allows non
% GUI functions (such as a function that just prints information) to be
% added to the registry. The registry itself thus conforms to this
% model. Each object has an associated "handling function" (so this
% function is the registry's handling function). The registry itself
% consists of object-handle/handling-function pairs.
%
% If an object and it's handling function are entered in the registry,
% then the object is said to be "forward registered", because the
% registry will now forward all location updates to that object, via
% it's handling function. The assummed syntax is:
% feval(Fcn,'SetCoords',xyz,h,hReg), where Fcn is the handling function
% for the GUI control identified by handle h, xyz are the new
% co-ordinates, and hReg is the handle of the registry.
%
% An optional extension is "back registration", whereby the GUI
% controls inform the registry of the new location when they are
% updated. All that's required is that the objects call the registry's
% 'SetCoords' function: bspm_XYZreg('SetCoords',xyz,hReg,hC), where hReg
% is the registry object's handle, and hC is the handle associated with
% the calling GUI control. The specification of the caller GUI control
% allows the registry to avoid circularities: If the object is "forward
% registered" for updates, then the registry function doesn't try to
% update the object which just updated the registry! (Similarly, the
% handle of the registry object, hReg, is passed to the handling
% function during forward XYZ updating, so that the handling function's
% 'SetCoords' facility can be constructed to accept XYZ updates from
% various sources, and only inform the registry if not called by the
% registry, and hence avoid circularities.)
%
% A framework is provided for "back" registration. Really all that is
% required is that the GUI controls know of the registry object (via
% it's handle hReg), and call the registry's 'SetCoords' facility when
% necessary. This can be done in many ways, but a simple structure is
% provided, mirroring that of the registry's operation. This framework
% assummes that the GUI controls identification object's 'UserData' is
% a structure with a field named 'hReg', which stores the handle of the
% registry (if back registered), or is empty (if not back registered,
% i.e. standalone). bspm_XYZreg provides utility functions for
% setting/unsetting this field, and for "cross registering" - that is
% both forward and back registration in one command. Cross registering
% involves adding the handle/function pair to the registry, and setting
% the registry handle in the GUI control object's 'UserData' 'hReg'
% field. It's up to the handling function to read the registry handle
% from it's objects 'UserData' and act accordingly. A simple example of
% such a function is provided in bspm_XYZreg_Ex2.m, illustrated below.
%
% SubFunctions are provided for getting and setting the current
% co-ordinates; adding and deleting handle/function pairs from the
% registry (forward registration and un-registration), setting and
% removing registry handle information from the 'hReg' field of the
% 'UserData' of a HG object (backward registration & un-registration);
% cross registration and unregistration (including pushing of current
% co-ordinates); setting and getting the current XYZ location. See the
% FORMAT statements and the example below...
%
%                           ----------------
% Example
% %-Create a window:
% F = figure;
% %-Create an object to hold the registry
% hReg = uicontrol(F,'Style','Text','String','hReg',...
%   'Position',[100 200 100 025],...
%   'FontName','Times','FontSize',14,'FontWeight','Bold',...
%   'HorizontalAlignment','Center');
% %-Setup M & D
% V = [65;87;26;02;02;04;33;53;08];
% M = [ [diag(V(4:6)), -(V(7:9).*V(4:6))]; [zeros(1,3) ,1]];
% D = V(1:3);
% %-Initialise a registry in this object, with initial co-ordinates [0;0;0]
% bspm_XYZreg('InitReg',hReg,M,D,[0;0;0])
% % (ans returns [0;0;0] confirming current co-ordinates
% %-Set co-ordinates to [10;10;10]
% bspm_XYZreg('SetCoords',[10,10,10],hReg)
% % (warns of co-ordinate rounding to [10,10,12], & returns ans as [10;10;12])
%
% %-Forward register a command window xyz reporting function: bspm_XYZreg_Ex1.m
% bspm_XYZreg('Add2Reg',hReg,0,'bspm_XYZreg_Ex1')
% % (ans returns new registry, containing just this handle/function pair
% %-Set co-ordinates to [0;10;12]
% [xyz,d] = bspm_XYZreg('SetCoords',[0,10,12],hReg);
% % (bspm_XYZreg_Ex1 called, and prints co-ordinates and handles)
% %-Have a peek at the registry information
% RD = get(hReg,'UserData')
% RD.xyz    %-The current point according to the registry
% RD.Reg    %-The nx2 cell array of handle/function pairs
%
% %-Create an example GUI XYZ control, using bspm_XYZreg_Ex2.m
% hB = bspm_XYZreg_Ex2('Create',M,D,xyz);
% % (A figure window with a button appears, whose label shows the current xyz
% %-Press the button, and enter new co-ordinates [0;0;0] in the Cmd window...
% % (...the button's internal notion of the current location is changed, but
% % (the registry isn't informed:
% bspm_XYZreg('GetCoords',hReg)
% (...returns [0;10;12])
% %-"Back" register the button
% bspm_XYZreg('SetReg',hB,hReg)
% %-Check the back registration
% if ( hReg == getfield(get(hB,'UserData'),'hReg') ), disp('yes!'), end
% %-Now press the button, and enter [0;0;0] again...
% % (...this time the registry is told, and the registry tells bspm_XYZreg_Ex1,
% % (which prints out the new co-ordinates!
% %-Forward register the button to receive updates from the registry
% nReg = bspm_XYZreg('Add2Reg',hReg,hB,'bspm_XYZreg_Ex2')
% % (The new registry is returned as nReg, showing two entries
% %-Set new registry co-ordinates to [10;10;12]
% [xyz,d] = bspm_XYZreg('SetCoords',[10;10;12],hReg);
% % (...the button updates too!
%
% %-Illustration of robustness: Delete the button & use the registry
% delete(hB)
% [xyz,d] = bspm_XYZreg('SetCoords',[10;10;12],hReg);
% % (...the invalid handle hB in the registry is ignored)
% %-Peek at the registry
% getfield(get(hReg,'UserData'),'Reg')
% %-Delete hB from the registry by "cleaning"
% bspm_XYZreg('CleanReg',hReg)
% % (...it's gone
%
% %-Make a new button and cross register
% hB = bspm_XYZreg_Ex2('Create',M,D)
% % (button created with default co-ordinates of [0;0;0]
% nReg = bspm_XYZreg('XReg',hReg,hB,'bspm_XYZreg_Ex2')
% % (Note that the registry pushes the current co-ordinates to the button
% %-Use the button & bspm_XYZreg('SetCoords'... at will!
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes, Chloe Hutton
% $Id: bspm_XYZreg.m 3664 2010-01-07 16:08:51Z volkmar $



%=======================================================================
switch lower(varargin{1}), case 'roundcoords'
%=======================================================================
% [xyz,d] = bspm_XYZreg('RoundCoords',xyz,M,D)
% [xyz,d] = bspm_XYZreg('RoundCoords',xyz,V)
if nargin<3, error('Insufficient arguments'), end
if nargin<4
    V = varargin{3};
    M = [ [diag(V(4:6)), -(V(7:9).*V(4:6))]; [zeros(1,3) ,1]];
    D = V(1:3);
else
    M = varargin{3};
    D = varargin{4};
end
    
%-Round xyz to coordinates of actual voxel centre
%-Do rounding in voxel coordinates & ensure within image size
%-Watch out for infinities!
%-----------------------------------------------------------------------
xyz  = [varargin{2}(:); 1];
xyz(isinf(xyz)) = 1e10*sign(xyz(isinf(xyz)));
rcp  = round(inv(M)*xyz);
rcp  = max([min([rcp';[D',1]]);[1,1,1,1]])';
rxyz = M*rcp;

%-Work out Euclidean distance between points xyz & rounded xyz
d = sqrt(sum((xyz-rxyz).^2));

varargout = {rxyz(1:3),d};



%=======================================================================
case 'findxyz'
%=======================================================================
% i = bspm_XYZreg('FindXYZ',xyz,XYZ)
if nargin<3, error('Insufficient arguments'), end
XYZ = varargin{3};
xyz = varargin{2};
    
%-Find XYZ = xyz
%-----------------------------------------------------------------------
i = find(all([XYZ(1,:)==xyz(1);XYZ(2,:)==xyz(2);XYZ(3,:)==xyz(3)],1));

varargout = {i};



%=======================================================================
case 'nearestxyz'
%=======================================================================
% [xyz,i,d] = bspm_XYZreg('NearestXYZ',xyz,XYZ)
if nargin<3, error('Insufficient arguments'), end
    
%-Find in XYZ nearest point to coordinates xyz (Euclidean distance) 
%-----------------------------------------------------------------------
[d,i] = min(bspm_XYZreg('Edist',varargin{2},varargin{3}));
varargout = {varargin{3}(:,i),i,d};



%=======================================================================
case 'edist'
%=======================================================================
% d = bspm_XYZreg('Edist',xyz,XYZ)
if nargin<3, error('Insufficient arguments'), end
    
%-Calculate (Euclidean) distances from pointlist co-ords to xyz
%-----------------------------------------------------------------------
varargout = {sqrt(sum([ (varargin{3}(1,:) - varargin{2}(1));...
            (varargin{3}(2,:) - varargin{2}(2));...
            (varargin{3}(3,:) - varargin{2}(3)) ].^2))};



%=======================================================================
case 'initreg'      % Initialise registry in handle h
%=======================================================================
% [hReg,xyz] = bspm_XYZreg('InitReg',hReg,M,D,xyz)
if nargin<5, xyz=[0;0;0]; else, xyz=varargin{5}; end
if nargin<4, error('Insufficient arguments'), end
D    = varargin{4};
M    = varargin{3};
hReg = varargin{2};

%-Check availability of hReg object for building a registry in
%-----------------------------------------------------------------------
if ~isempty(get(hReg,'UserData')), error('Object already has UserData...'), end
if ~isempty(get(hReg,'Tag')), error('Object already ''Tag''ed...'), end

%-Check co-ordinates are in range
%-----------------------------------------------------------------------
[xyz,d] = bspm_XYZreg('RoundCoords',xyz,M,D);
if d>0 & nargout<2, warning(sprintf('%s: Co-ords rounded to neatest voxel center: Discrepancy %.2f',mfilename,d)), end

%-Set up registry
%-----------------------------------------------------------------------
RD = struct('xyz',xyz,'M',M,'D',D,'Reg',[]);
RD.Reg = {};
set(hReg,'Tag','hReg','UserData',RD)

%-Return current co-ordinates
%-----------------------------------------------------------------------
varargout = {hReg,xyz};



%=======================================================================
case 'uninitreg'    % UnInitialise registry in handle hReg
%=======================================================================
% bspm_XYZreg('UnInitReg',hReg)
hReg = varargin{2};
if ~strcmp(get(hReg,'Tag'),'hReg'), warning('Not an XYZ registry'), return, end
set(hReg,'Tag','','UserData',[])



%=======================================================================
case 'getcoords'    % Get current co-ordinates
%=======================================================================
% xyz = bspm_XYZreg('GetCoords',hReg)
if nargin<2, hReg=bspm_XYZreg('FindReg'); else, hReg=varargin{2}; end
if ~ishandle(hReg), error('Invalid object handle'), end
if ~strcmp(get(hReg,'Tag'),'hReg'), error('Not a registry'), end
varargout = {getfield(get(hReg,'UserData'),'xyz')};



%=======================================================================
case 'setcoords'    % Set co-ordinates & update registered functions
%=======================================================================
% [xyz,d] = bspm_XYZreg('SetCoords',xyz,hReg,hC,Reg)
% d returned empty if didn't check, warning printed if d not asked for & round
% Don't check if callerhandle specified (speed)
% If Registry cell array Reg is specified, then only these handles are updated
hC=0; mfn=''; if nargin>=4
    if ~ischar(varargin{4}), hC=varargin{4}; else mfn=varargin{4}; end
end
hReg = varargin{3};

%-Check validity of hReg registry handle
%-----------------------------------------------------------------------
%-Return if hReg empty, in case calling objects functions don't check isempty
if isempty(hReg), return, end
%-Check validity of hReg registry handle, correct calling objects if necc.
if ~ishandle(hReg)
    str = sprintf('%s: Invalid registry handle (%.4f)',mfilename,hReg);
    if hC>0
        %-Remove hReg from caller
        bspm_XYZreg('SetReg',hC,[])
        str = [str,sprintf('\n\t\t\t...removed from caller (%.4f)',hC)];
    end
    warning(str)
    return
end
xyz  = varargin{2};

RD      = get(hReg,'UserData');

%-Check validity of coords only when called without a caller handle
%-----------------------------------------------------------------------
if hC<=0
    [xyz,d] = bspm_XYZreg('RoundCoords',xyz,RD.M,RD.D);
    if d>0 & nargout<2, warning(sprintf(...
        '%s: Co-ords rounded to neatest voxel center: Discrepancy %.2f',...
        mfilename,d)), end
else
    d = 0;
end

%-Sort out valid handles, eliminate caller handle, update co-ords with
% registered handles via their functions
%-----------------------------------------------------------------------
if nargin<5
    RD.Reg = bspm_XYZreg('VReg',RD.Reg);
    Reg    = RD.Reg;
else
    Reg = bspm_XYZreg('VReg',varargin{5});
end
if hC>0 & length(Reg), Reg(find([Reg{:,1}]==varargin{4}),:) = []; end
for i = 1:size(Reg,1)
    feval(Reg{i,2},'SetCoords',xyz,Reg{i,1},hReg);
end

%-Update registry (if using hReg) with location & cleaned Reg cellarray
%-----------------------------------------------------------------------
if nargin<5
    RD.xyz  = xyz;
    set(hReg,'UserData',RD)
end

varargout = {xyz,d};

if ~strcmp(mfn,'spm_graph')
    sHdl=findobj(0,'Tag','SPMGraphSatelliteFig');
    axHdl=findobj(sHdl,'Type','axes','Tag','SPMGraphSatelliteAxes');
    %tag for true axis, as legend is of type axis, too
    for j=1:length(axHdl)
        autoinp=get(axHdl(j),'UserData');
        if ~isempty(autoinp), spm_graph([],[],hReg,axHdl(j)); end
    end
end


%=======================================================================
case 'xreg'     % Cross register object handles & functions
%=======================================================================
% nReg = bspm_XYZreg('XReg',hReg,{h,Fcn}pairs)
if nargin<4, error('Insufficient arguments'), end
hReg = varargin{2};

%-Quick check of registry handle
%-----------------------------------------------------------------------
if isempty(hReg),   warning('Empty registry handle'), return, end
if ~ishandle(hReg), warning('Invalid registry handle'), return, end

%-Condition nReg cell array & check validity of handles to be registered
%-----------------------------------------------------------------------
nReg = varargin(3:end);
if mod(length(nReg),2), error('Registry items must be in pairs'), end
if length(nReg)>2, nReg = reshape(nReg,length(nReg)/2,2)'; end
nReg = bspm_XYZreg('VReg',nReg,'Warn');

%-Set hReg registry link for registry candidates (Back registration)
%-----------------------------------------------------------------------
for i = 1:size(nReg,1)
    bspm_XYZreg('SetReg',nReg{i,1},hReg);
end

%-Append registry candidates to existing registry & write back to hReg
%-----------------------------------------------------------------------
RD     = get(hReg,'UserData');
Reg    = RD.Reg;
Reg    = cat(1,Reg,nReg);
RD.Reg = Reg;
set(hReg,'UserData',RD)

%-Synch co-ordinates of newly registered objects
%-----------------------------------------------------------------------
bspm_XYZreg('SetCoords',RD.xyz,hReg,hReg,nReg);

varargout = {Reg};



%=======================================================================
case 'add2reg'      % Add handle(s) & function(s) to registry
%=======================================================================
% nReg = bspm_XYZreg('Add2Reg',hReg,{h,Fcn}pairs)
if nargin<4, error('Insufficient arguments'), end
hReg = varargin{2};

%-Quick check of registry handle
%-----------------------------------------------------------------------
if isempty(hReg),   warning('Empty registry handle'), return, end
if ~ishandle(hReg), warning('Invalid registry handle'), return, end

%-Condition nReg cell array & check validity of handles to be registered
%-----------------------------------------------------------------------
nReg = varargin(3:end);
if mod(length(nReg),2), error('Registry items must be in pairs'), end
if length(nReg)>2, nReg = reshape(nReg,length(nReg)/2,2)'; end
nReg = bspm_XYZreg('VReg',nReg,'Warn');

%-Append to existing registry & put back in registry object
%-----------------------------------------------------------------------
RD     = get(hReg,'UserData');
Reg    = RD.Reg;
Reg    = cat(1,Reg,nReg);
RD.Reg = Reg;
set(hReg,'UserData',RD)

varargout = {Reg};



%=======================================================================
case 'setreg'           %-Set registry field of object's UserData
%=======================================================================
% bspm_XYZreg('SetReg',h,hReg)
if nargin<3, error('Insufficient arguments'), end
h    = varargin{2};
hReg = varargin{3};
if ( ~ishandle(h) | h==0 ), return, end
UD = get(h,'UserData');
if ~isstruct(UD) | ~any(strcmp(fieldnames(UD),'hReg'))
    error('No UserData structure with hReg field for this object')
end
UD.hReg = hReg;
set(h,'UserData',UD)



%=======================================================================
case 'unxreg'       % Un-cross register object handles & functions
%=======================================================================
% nReg = bspm_XYZreg('unXReg',hReg,hD1,hD2,hD3,...)
if nargin<3, error('Insufficient arguments'), end
hD   = [varargin{3:end}];
hReg = varargin{2};

%-Get Registry information
%-----------------------------------------------------------------------
RD         = get(hReg,'UserData');
Reg        = RD.Reg;

%-Find registry entires to delete
%-----------------------------------------------------------------------
[null,i,e] = intersect([Reg{:,1}],hD);
hD(e)      = [];
dReg       = bspm_XYZreg('VReg',Reg(i,:));
Reg(i,:)   = [];
if length(hD), warning('Not all handles were in registry'), end

%-Write back new registry
%-----------------------------------------------------------------------
RD.Reg = Reg;
set(hReg,'UserData',RD)

%-UnSet hReg registry link for hD's still existing (Back un-registration)
%-----------------------------------------------------------------------
for i = 1:size(dReg,1)
    bspm_XYZreg('SetReg',dReg{i,1},[]);
end

varargout = {Reg};



%=======================================================================
case 'del2reg'      % Delete handle(s) & function(s) from registry
%=======================================================================
% nReg = bspm_XYZreg('Del2Reg',hReg,hD)
if nargin<3, error('Insufficient arguments'), end
hD   = [varargin{3:end}];
hReg = varargin{2};

%-Get Registry information
%-----------------------------------------------------------------------
RD         = get(hReg,'UserData');
Reg        = RD.Reg;

%-Find registry entires to delete
%-----------------------------------------------------------------------
[null,i,e] = intersect([Reg{:,1}],hD);
Reg(i,:)   = [];
hD(e)      = [];
if length(hD), warning('Not all handles were in registry'), end

%-Write back new registry
%-----------------------------------------------------------------------
RD.Reg = Reg;
set(hReg,'UserData',RD)

varargout = {Reg};



%=======================================================================
case 'unsetreg'         %-Unset registry field of object's UserData
%=======================================================================
% bspm_XYZreg('UnSetReg',h)
if nargin<2, error('Insufficient arguments'), end
bspm_XYZreg('SetReg',varargin{2},[])



%=======================================================================
case 'cleanreg'     % Clean invalid handles from registry
%=======================================================================
% bspm_XYZreg('CleanReg',hReg)
%if ~strcmp(get(hReg,'Tag'),'hReg'), error('Not a registry'), end
hReg = varargin{2};
RD = get(hReg,'UserData');
RD.Reg = bspm_XYZreg('VReg',RD.Reg,'Warn');
set(hReg,'UserData',RD)


%=======================================================================
case 'vreg'     % Prune invalid handles from registry cell array
%=======================================================================
% Reg = bspm_XYZreg('VReg',Reg,Warn)
if nargin<3, Warn=0; else, Warn=1; end
Reg = varargin{2};
if isempty(Reg), varargout={Reg}; return, end
i = find(~ishandle([Reg{:,1}]));
%-***check existance of handling functions : exist('','file')?
if Warn & length(i), warning([...
    sprintf('%s: Disregarding invalid registry handles:\n\t',...
        mfilename),sprintf('%.4f',Reg{i,1})]), end
Reg(i,:)  = [];
varargout = {Reg};



%=======================================================================
case 'findreg'          % Find/check registry object
%=======================================================================
% hReg = bspm_XYZreg('FindReg',h)
if nargin<2, h=get(0,'CurrentFigure'); else, h=varargin{2}; end
if ischar(h), h=spm_figure('FindWin',h); end
if ~ishandle(h), error('invalid handle'), end
if ~strcmp(get(h,'Tag'),'hReg'), h=findobj(h,'Tag','hReg'); end
if isempty(h), error('Registry object not found'), end
if length(h)>1, error('Multiple registry objects found'), end
varargout = {h};



%=======================================================================
otherwise
%=======================================================================
warning('Unknown action string')

%=======================================================================
end

% | screencapture (file exchange)
% ========================================================================= 
function imageData = screencapture(varargin)
% screencapture - get a screen-capture of a figure frame, component handle, or screen area rectangle
%
% ScreenCapture gets a screen-capture of any Matlab GUI handle (including desktop, 
% figure, axes, image or uicontrol), or a specified area rectangle located relative to
% the specified handle. Screen area capture is possible by specifying the root (desktop)
% handle (=0). The output can be either to an image file or to a Matlab matrix (useful
% for displaying via imshow() or for further processing) or to the system clipboard.
% This utility also enables adding a toolbar button for easy interactive screen-capture.
%
% Syntax:
%    imageData = screencapture(handle, position, target, 'PropName',PropValue, ...)
%
% Input Parameters:
%    handle   - optional handle to be used for screen-capture origin.
%                 If empty/unsupplied then current figure (gcf) will be used.
%    position - optional position array in pixels: [x,y,width,height].
%                 If empty/unsupplied then the handle's position vector will be used.
%                 If both handle and position are empty/unsupplied then the position
%                   will be retrieved via interactive mouse-selection.
%                 If handle is an image, then position is in data (not pixel) units, so the
%                   captured region remains the same after figure/axes resize (like imcrop)
%    target   - optional filename for storing the screen-capture, or the 'clipboard' string.
%                 If empty/unsupplied then no output to file will be done.
%                 The file format will be determined from the extension (JPG/PNG/...).
%                 Supported formats are those supported by the imwrite function.
%    'PropName',PropValue - 
%               optional list of property pairs (e.g., 'target','myImage.png','pos',[10,20,30,40],'handle',gca)
%               PropNames may be abbreviated and are case-insensitive.
%               PropNames may also be given in whichever order.
%               Supported PropNames are:
%                 - 'handle'    (default: gcf handle)
%                 - 'position'  (default: gcf position array)
%                 - 'target'    (default: '')
%                 - 'toolbar'   (figure handle; default: gcf)
%                      this adds a screen-capture button to the figure's toolbar
%                      If this parameter is specified, then no screen-capture
%                        will take place and the returned imageData will be [].
%
% Output parameters:
%    imageData - image data in a format acceptable by the imshow function
%                  If neither target nor imageData were specified, the user will be
%                    asked to interactively specify the output file.
%
% Examples:
%    imageData = screencapture;  % interactively select screen-capture rectangle
%    imageData = screencapture(hListbox);  % capture image of a uicontrol
%    imageData = screencapture(0,  [20,30,40,50]);  % capture a small desktop region
%    imageData = screencapture(gcf,[20,30,40,50]);  % capture a small figure region
%    imageData = screencapture(gca,[10,20,30,40]);  % capture a small axes region
%      imshow(imageData);  % display the captured image in a matlab figure
%      imwrite(imageData,'myImage.png');  % save the captured image to file
%    img = imread('cameraman.tif');
%      hImg = imshow(img);
%      screencapture(hImg,[60,35,140,80]);  % capture a region of an image
%    screencapture(gcf,[],'myFigure.jpg');  % capture the entire figure into file
%    screencapture('handle',gcf,'target','myFigure.jpg');  % same as previous
%    screencapture('handle',gcf,'target','clipboard');     % copy to clipboard
%    screencapture('toolbar',gcf);  % adds a screen-capture button to gcf's toolbar
%    screencapture('toolbar',[],'target','sc.bmp'); % same with default output filename
%
% Technical description:
%    http://UndocumentedMatlab.com/blog/screencapture-utility/
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% See also:
%    imshow, imwrite, print
%
% Release history:
%    1.7 2014-04-28: Fixed bug when capturing interactive selection
%    1.6 2014-04-22: Only enable image formats when saving to an unspecified file via uiputfile
%    1.5 2013-04-18: Fixed bug in capture of non-square image; fixes for Win64
%    1.4 2013-01-27: Fixed capture of Desktop (root); enabled rbbox anywhere on desktop (not necesarily in a Matlab figure); enabled output to clipboard (based on Jiro Doke's imclipboard utility); edge-case fixes; added Java compatibility check
%    1.3 2012-07-23: Capture current object (uicontrol/axes/figure) if w=h=0 (e.g., by clicking a single point); extra input args sanity checks; fix for docked windows and image axes; include axes labels & ticks by default when capturing axes; use data-units position vector when capturing images; many edge-case fixes
%    1.2 2011-01-16: another performance boost (thanks to Jan Simon); some compatibility fixes for Matlab 6.5 (untested)
%    1.1 2009-06-03: Handle missing output format; performance boost (thanks to Urs); fix minor root-handle bug; added toolbar button option
%    1.0 2009-06-02: First version posted on <a href="http://www.mathworks.com/matlabcentral/fileexchange/authors/27420">MathWorks File Exchange</a>

% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.

% Programmed and Copyright by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.7 $  $Date: 2014/04/28 21:10:12 $

    % Ensure that java awt is enabled...
    if ~usejava('awt')
        error('YMA:screencapture:NeedAwt','ScreenCapture requires Java to run.');
    end

    % Ensure that our Java version supports the Robot class (requires JVM 1.3+)
    try
        robot = java.awt.Robot; %#ok<NASGU>
    catch
        uiwait(msgbox({['Your Matlab installation is so old that its Java engine (' version('-java') ...
                        ') does not have a java.awt.Robot class. '], ' ', ...
                        'Without this class, taking a screen-capture is impossible.', ' ', ...
                        'So, either install JVM 1.3 or higher, or use a newer Matlab release.'}, ...
                        'ScreenCapture', 'warn'));
        if nargout, imageData = [];  end
        return;
    end

    % Process optional arguments
    paramsStruct = processArgs(varargin{:});
    
    % If toolbar button requested, add it and exit
    if ~isempty(paramsStruct.toolbar)
        
        % Add the toolbar button
        addToolbarButton(paramsStruct);
        
        % Return the figure to its pre-undocked state (when relevant)
        redockFigureIfRelevant(paramsStruct);
        
        % Exit immediately (do NOT take a screen-capture)
        if nargout,  imageData = [];  end
        return;
    end
    
    % Convert position from handle-relative to desktop Java-based pixels
    [paramsStruct, msgStr] = convertPos(paramsStruct);
    
    % Capture the requested screen rectangle using java.awt.Robot
    imgData = getScreenCaptureImageData(paramsStruct.position);
    
    % Return the figure to its pre-undocked state (when relevant)
    redockFigureIfRelevant(paramsStruct);
    
    % Save image data in file or clipboard, if specified
    if ~isempty(paramsStruct.target)
        if strcmpi(paramsStruct.target,'clipboard')
            if ~isempty(imgData)
                imclipboard(imgData);
            else
                msgbox('No image area selected - not copying image to clipboard','ScreenCapture','warn');
            end
        else  % real filename
            if ~isempty(imgData)
                imwrite(imgData,paramsStruct.target);
            else
                msgbox(['No image area selected - not saving image file ' paramsStruct.target],'ScreenCapture','warn');
            end
        end
    end

    % Return image raster data to user, if requested
    if nargout
        imageData = imgData;
        
    % If neither output formats was specified (neither target nor output data)
    elseif isempty(paramsStruct.target) & ~isempty(imgData)  %#ok ML6
        % Ask the user to specify a file
        %error('YMA:screencapture:noOutput','No output specified for ScreenCapture: specify the output filename and/or output data');
        %format = '*.*';
        formats = imformats;
        for idx = 1 : numel(formats)
            ext = sprintf('*.%s;',formats(idx).ext{:});
            format(idx,1:2) = {ext(1:end-1), formats(idx).description}; %#ok<AGROW>
        end
        [filename,pathname] = uiputfile(format,'Save screen capture as');
        if ~isequal(filename,0) & ~isequal(pathname,0)  %#ok Matlab6 compatibility
            imwrite(imgData,fullfile(pathname,filename));
        else
            % TODO - copy to clipboard
        end
    end

    % Display msgStr, if relevant
    if ~isempty(msgStr)
        uiwait(msgbox(msgStr,'ScreenCapture'));
        drawnow; pause(0.05);  % time for the msgbox to disappear
    end

    return;  % debug breakpoint
function paramsStruct = processArgs(varargin)

    % Get the properties in either direct or P-V format
    [regParams, pvPairs] = parseparams(varargin);

    % Now process the optional P-V params
    try
        % Initialize
        paramName = [];
        paramsStruct = [];
        paramsStruct.handle = [];
        paramsStruct.position = [];
        paramsStruct.target = '';
        paramsStruct.toolbar = [];
        paramsStruct.wasDocked = 0;       % no false available in ML6
        paramsStruct.wasInteractive = 0;  % no false available in ML6

        % Parse the regular (non-named) params in recption order
        if ~isempty(regParams) & (isempty(regParams{1}) | ishandle(regParams{1}(1)))  %#ok ML6
            paramsStruct.handle = regParams{1};
            regParams(1) = [];
        end
        if ~isempty(regParams) & isnumeric(regParams{1}) & (length(regParams{1}) == 4)  %#ok ML6
            paramsStruct.position = regParams{1};
            regParams(1) = [];
        end
        if ~isempty(regParams) & ischar(regParams{1})  %#ok ML6
            paramsStruct.target = regParams{1};
        end

        % Parse the optional param PV pairs
        supportedArgs = {'handle','position','target','toolbar'};
        while ~isempty(pvPairs)

            % Disregard empty propNames (may be due to users mis-interpretting the syntax help)
            while ~isempty(pvPairs) & isempty(pvPairs{1})  %#ok ML6
                pvPairs(1) = [];
            end
            if isempty(pvPairs)
                break;
            end

            % Ensure basic format is valid
            paramName = '';
            if ~ischar(pvPairs{1})
                error('YMA:screencapture:invalidProperty','Invalid property passed to ScreenCapture');
            elseif length(pvPairs) == 1
                if isempty(paramsStruct.target)
                    paramsStruct.target = pvPairs{1};
                    break;
                else
                    error('YMA:screencapture:noPropertyValue',['No value specified for property ''' pvPairs{1} '''']);
                end
            end

            % Process parameter values
            paramName  = pvPairs{1};
            if strcmpi(paramName,'filename')  % backward compatibility
                paramName = 'target';
            end
            paramValue = pvPairs{2};
            pvPairs(1:2) = [];
            idx = find(strncmpi(paramName,supportedArgs,length(paramName)));
            if ~isempty(idx)
                %paramsStruct.(lower(supportedArgs{idx(1)})) = paramValue;  % incompatible with ML6
                paramsStruct = setfield(paramsStruct, lower(supportedArgs{idx(1)}), paramValue);  %#ok ML6

                % If 'toolbar' param specified, then it cannot be left empty - use gcf
                if strncmpi(paramName,'toolbar',length(paramName)) & isempty(paramsStruct.toolbar)  %#ok ML6
                    paramsStruct.toolbar = getCurrentFig;
                end

            elseif isempty(paramsStruct.target)
                paramsStruct.target = paramName;
                pvPairs = {paramValue, pvPairs{:}};  %#ok (more readable this way, although a bit less efficient...)

            else
                supportedArgsStr = sprintf('''%s'',',supportedArgs{:});
                error('YMA:screencapture:invalidProperty','%s \n%s', ...
                      'Invalid property passed to ScreenCapture', ...
                      ['Supported property names are: ' supportedArgsStr(1:end-1)]);
            end
        end  % loop pvPairs

    catch
        if ~isempty(paramName),  paramName = [' ''' paramName ''''];  end
        error('YMA:screencapture:invalidProperty','Error setting ScreenCapture property %s:\n%s',paramName,lasterr); %#ok<LERR>
    end
function [paramsStruct, msgStr] = convertPos(paramsStruct)
    msgStr = '';
    try
        % Get the screen-size for later use
        screenSize = get(0,'ScreenSize');

        % Get the containing figure's handle
        hParent = paramsStruct.handle;
        if isempty(paramsStruct.handle)
            paramsStruct.hFigure = getCurrentFig;
            hParent = paramsStruct.hFigure;
        else
            paramsStruct.hFigure = ancestor(paramsStruct.handle,'figure');
        end

        % To get the acurate pixel position, the figure window must be undocked
        try
            if strcmpi(get(paramsStruct.hFigure,'WindowStyle'),'docked')
                set(paramsStruct.hFigure,'WindowStyle','normal');
                drawnow; pause(0.25);
                paramsStruct.wasDocked = 1;  % no true available in ML6
            end
        catch
            % never mind - ignore...
        end

        % The figure (if specified) must be in focus
        if ~isempty(paramsStruct.hFigure) & ishandle(paramsStruct.hFigure)  %#ok ML6
            isFigureValid = 1;  % no true available in ML6
            figure(paramsStruct.hFigure);
        else
            isFigureValid = 0;  % no false available in ML6
        end

        % Flush all graphic events to ensure correct rendering
        drawnow; pause(0.01);

        % No handle specified
        wasPositionGiven = 1;  % no true available in ML6
        if isempty(paramsStruct.handle)
            
            % Set default handle, if not supplied
            paramsStruct.handle = paramsStruct.hFigure;
            
            % If position was not specified, get it interactively using RBBOX
            if isempty(paramsStruct.position)
                [paramsStruct.position, jFrameUsed, msgStr] = getInteractivePosition(paramsStruct.hFigure); %#ok<ASGLU> jFrameUsed is unused
                paramsStruct.wasInteractive = 1;  % no true available in ML6
                wasPositionGiven = 0;  % no false available in ML6
            end
            
        elseif ~ishandle(paramsStruct.handle)
            % Handle was supplied - ensure it is a valid handle
            error('YMA:screencapture:invalidHandle','Invalid handle passed to ScreenCapture');
            
        elseif isempty(paramsStruct.position)
            % Handle was supplied but position was not, so use the handle's position
            paramsStruct.position = getPixelPos(paramsStruct.handle);
            paramsStruct.position(1:2) = 0;
            wasPositionGiven = 0;  % no false available in ML6
            
        elseif ~isnumeric(paramsStruct.position) | (length(paramsStruct.position) ~= 4)  %#ok ML6
            % Both handle & position were supplied - ensure a valid pixel position vector
            error('YMA:screencapture:invalidPosition','Invalid position vector passed to ScreenCapture: \nMust be a [x,y,w,h] numeric pixel array');
        end
        
        % Capture current object (uicontrol/axes/figure) if w=h=0 (single-click in interactive mode)
        if paramsStruct.position(3)<=0 | paramsStruct.position(4)<=0  %#ok ML6
            %TODO - find a way to single-click another Matlab figure (the following does not work)
            %paramsStruct.position = getPixelPos(ancestor(hittest,'figure'));
            paramsStruct.position = getPixelPos(paramsStruct.handle);
            paramsStruct.position(1:2) = 0;
            paramsStruct.wasInteractive = 0;  % no false available in ML6
            wasPositionGiven = 0;  % no false available in ML6
        end

        % First get the parent handle's desktop-based Matlab pixel position
        parentPos = [0,0,0,0];
        dX = 0;
        dY = 0;
        dW = 0;
        dH = 0;
        if ~isa(handle(hParent),'figure')
            % Get the reguested component's pixel position
            parentPos = getPixelPos(hParent, 1);  % no true available in ML6

            % Axes position inaccuracy estimation
            deltaX = 3;
            deltaY = -1;
            
            % Fix for images
            %isAxes  = isa(handle(hParent),'axes');
            isImage = isa(handle(hParent),'image');
            if isImage  % | (isAxes & strcmpi(get(hParent,'YDir'),'reverse'))  %#ok ML6

                % Compensate for resized image axes
                hAxes = get(hParent,'Parent');
                if all(get(hAxes,'DataAspectRatio')==1)  % sanity check: this is the normal behavior
                    % Note 18/4/2013: the following fails for non-square images
                    %actualImgSize = min(parentPos(3:4));
                    %dX = (parentPos(3) - actualImgSize) / 2;
                    %dY = (parentPos(4) - actualImgSize) / 2;
                    %parentPos(3:4) = actualImgSize;

                    % The following should work for all types of images
                    actualImgSize = size(get(hParent,'CData'));
                    dX = (parentPos(3) - min(parentPos(3),actualImgSize(2))) / 2;
                    dY = (parentPos(4) - min(parentPos(4),actualImgSize(1))) / 2;
                    parentPos(3:4) = actualImgSize([2,1]);
                    %parentPos(3) = max(parentPos(3),actualImgSize(2));
                    %parentPos(4) = max(parentPos(4),actualImgSize(1));
                end

                % Fix user-specified img positions (but not auto-inferred ones)
                if wasPositionGiven

                    % In images, use data units rather than pixel units
                    % Reverse the YDir
                    ymax = max(get(hParent,'YData'));
                    paramsStruct.position(2) = ymax - paramsStruct.position(2) - paramsStruct.position(4);

                    % Note: it would be best to use hgconvertunits, but:
                    % ^^^^  (1) it fails on Matlab 6, and (2) it doesn't accept Data units
                    %paramsStruct.position = hgconvertunits(hFig, paramsStruct.position, 'Data', 'pixel', hParent);  % fails!
                    xLims = get(hParent,'XData');
                    yLims = get(hParent,'YData');
                    xPixelsPerData = parentPos(3) / (diff(xLims) + 1);
                    yPixelsPerData = parentPos(4) / (diff(yLims) + 1);
                    paramsStruct.position(1) = round((paramsStruct.position(1)-xLims(1)) * xPixelsPerData);
                    paramsStruct.position(2) = round((paramsStruct.position(2)-yLims(1)) * yPixelsPerData + 2*dY);
                    paramsStruct.position(3) = round( paramsStruct.position(3) * xPixelsPerData);
                    paramsStruct.position(4) = round( paramsStruct.position(4) * yPixelsPerData);

                    % Axes position inaccuracy estimation
                    if strcmpi(computer('arch'),'win64')
                        deltaX = 7;
                        deltaY = -7;
                    else
                        deltaX = 3;
                        deltaY = -3;
                    end
                    
                else  % axes/image position was auto-infered (entire image)
                    % Axes position inaccuracy estimation
                    if strcmpi(computer('arch'),'win64')
                        deltaX = 6;
                        deltaY = -6;
                    else
                        deltaX = 2;
                        deltaY = -2;
                    end
                    dW = -2*dX;
                    dH = -2*dY;
                end
            end

            %hFig = ancestor(hParent,'figure');
            hParent = paramsStruct.hFigure;

        elseif paramsStruct.wasInteractive  % interactive figure rectangle

            % Compensate for 1px rbbox inaccuracies
            deltaX = 2;
            deltaY = -2;

        else  % non-interactive figure

            % Compensate 4px figure boundaries = difference betweeen OuterPosition and Position
            deltaX = -1;
            deltaY = 1;
        end
        %disp(paramsStruct.position)  % for debugging
        
        % Now get the pixel position relative to the monitor
        figurePos = getPixelPos(hParent);
        desktopPos = figurePos + parentPos;

        % Now convert to Java-based pixels based on screen size
        % Note: multiple monitors are automatically handled correctly, since all
        % ^^^^  Java positions are relative to the main monitor's top-left corner
        javaX  = desktopPos(1) + paramsStruct.position(1) + deltaX + dX;
        javaY  = screenSize(4) - desktopPos(2) - paramsStruct.position(2) - paramsStruct.position(4) + deltaY + dY;
        width  = paramsStruct.position(3) + dW;
        height = paramsStruct.position(4) + dH;
        paramsStruct.position = round([javaX, javaY, width, height]);
        %paramsStruct.position

        % Ensure the figure is at the front so it can be screen-captured
        if isFigureValid
            figure(hParent);
            drawnow;
            pause(0.02);
        end
    catch
        % Maybe root/desktop handle (root does not have a 'Position' prop so getPixelPos croaks
        if isequal(double(hParent),0)  % =root/desktop handle;  handles case of hParent=[]
            javaX = paramsStruct.position(1) - 1;
            javaY = screenSize(4) - paramsStruct.position(2) - paramsStruct.position(4) - 1;
            paramsStruct.position = [javaX, javaY, paramsStruct.position(3:4)];
        end
    end
function [positionRect, jFrameUsed, msgStr] = getInteractivePosition(hFig)
    msgStr = '';
    try
        % First try the invisible-figure approach, in order to
        % enable rbbox outside any existing figure boundaries
        f = figure('units','pixel','pos',[-100,-100,10,10],'HitTest','off');
        drawnow; pause(0.01);
        jf = get(handle(f),'JavaFrame');
        try
            jWindow = jf.fFigureClient.getWindow;
        catch
            try
                jWindow = jf.fHG1Client.getWindow;
            catch
                jWindow = jf.getFigurePanelContainer.getParent.getTopLevelAncestor;
            end
        end
        com.sun.awt.AWTUtilities.setWindowOpacity(jWindow,0.05);  %=nearly transparent (not fully so that mouse clicks are captured)
        jWindow.setMaximized(1);  % no true available in ML6
        jFrameUsed = 1;  % no true available in ML6
        msg = {'Mouse-click and drag a bounding rectangle for screen-capture ' ...
               ... %'or single-click any Matlab figure to capture the entire figure.' ...
               };
    catch
        % Something failed, so revert to a simple rbbox on a visible figure
        try delete(f); drawnow; catch, end  %Cleanup...
        jFrameUsed = 0;  % no false available in ML6
        msg = {'Mouse-click within any Matlab figure and then', ...
               'drag a bounding rectangle for screen-capture,', ...
               'or single-click to capture the entire figure'};
    end
    uiwait(msgbox(msg,'ScreenCapture'));
    
    k = waitforbuttonpress;  %#ok k is unused
    %hFig = getCurrentFig;
    %p1 = get(hFig,'CurrentPoint');
    positionRect = rbbox;
    %p2 = get(hFig,'CurrentPoint');

    if jFrameUsed
        jFrameOrigin = getPixelPos(f);
        delete(f); drawnow;
        try
            figOrigin = getPixelPos(hFig);
        catch  % empty/invalid hFig handle
            figOrigin = [0,0,0,0];
        end
    else
        if isempty(hFig)
            jFrameOrigin = getPixelPos(gcf);
        else
            jFrameOrigin = [0,0,0,0];
        end
        figOrigin = [0,0,0,0];
    end
    positionRect(1:2) = positionRect(1:2) + jFrameOrigin(1:2) - figOrigin(1:2);

    if prod(positionRect(3:4)) > 0
        msgStr = sprintf('%dx%d area captured',positionRect(3),positionRect(4));
    end
function hFig = getCurrentFig
    oldState = get(0,'showHiddenHandles');
    set(0,'showHiddenHandles','on');
    hFig = get(0,'CurrentFigure');
    set(0,'showHiddenHandles',oldState);
function hObj = ancestor(hObj,type)
    if ~isempty(hObj) & ishandle(hObj)  %#ok for Matlab 6 compatibility
        try
            hObj = get(hObj,'Ancestor');
        catch
            % never mind...
        end
        try
            %if ~isa(handle(hObj),type)  % this is best but always returns 0 in Matlab 6!
            %if ~isprop(hObj,'type') | ~strcmpi(get(hObj,'type'),type)  % no isprop() in ML6!
            try
                objType = get(hObj,'type');
            catch
                objType = '';
            end
            if ~strcmpi(objType,type)
                try
                    parent = get(handle(hObj),'parent');
                catch
                    parent = hObj.getParent;  % some objs have no 'Parent' prop, just this method...
                end
                if ~isempty(parent)  % empty parent means root ancestor, so exit
                    hObj = ancestor(parent,type);
                end
            end
        catch
            % never mind...
        end
    end
function pos = getPos(hObj,field,units)
    % Matlab 6 did not have hgconvertunits so use the old way...
    oldUnits = get(hObj,'units');
    if strcmpi(oldUnits,units)  % don't modify units unless we must!
        pos = get(hObj,field);
    else
        set(hObj,'units',units);
        pos = get(hObj,field);
        set(hObj,'units',oldUnits);
    end
function pos = getPixelPos(hObj,varargin)
    persistent originalObj
    try
        stk = dbstack;
        if ~strcmp(stk(2).name,'getPixelPos')
            originalObj = hObj;
        end

        if isa(handle(hObj),'figure') %| isa(handle(hObj),'axes')
        %try
            pos = getPos(hObj,'OuterPosition','pixels');
        else  %catch
            % getpixelposition is unvectorized unfortunately!
            pos = getpixelposition(hObj,varargin{:});

            % add the axes labels/ticks if relevant (plus a tiny margin to fix 2px label/title inconsistencies)
            if isa(handle(hObj),'axes') & ~isa(handle(originalObj),'image')  %#ok ML6
                tightInsets = getPos(hObj,'TightInset','pixel');
                pos = pos + tightInsets.*[-1,-1,1,1] + [-1,1,1+tightInsets(1:2)];
            end
        end
    catch
        try
            % Matlab 6 did not have getpixelposition nor hgconvertunits so use the old way...
            pos = getPos(hObj,'Position','pixels');
        catch
            % Maybe the handle does not have a 'Position' prop (e.g., text/line/plot) - use its parent
            pos = getPixelPos(get(hObj,'parent'),varargin{:});
        end
    end

    % Handle the case of missing/invalid/empty HG handle
    if isempty(pos)
        pos = [0,0,0,0];
    end
function addToolbarButton(paramsStruct)
    % Ensure we have a valid toolbar handle
    hFig = ancestor(paramsStruct.toolbar,'figure');
    if isempty(hFig)
        error('YMA:screencapture:badToolbar','the ''Toolbar'' parameter must contain a valid GUI handle');
    end
    set(hFig,'ToolBar','figure');
    hToolbar = findall(hFig,'type','uitoolbar');
    if isempty(hToolbar)
        error('YMA:screencapture:noToolbar','the ''Toolbar'' parameter must contain a figure handle possessing a valid toolbar');
    end
    hToolbar = hToolbar(1);  % just in case there are several toolbars... - use only the first

    % Prepare the camera icon
    icon = ['3333333333333333'; ...
            '3333333333333333'; ...
            '3333300000333333'; ...
            '3333065556033333'; ...
            '3000000000000033'; ...
            '3022222222222033'; ...
            '3022220002222033'; ...
            '3022203110222033'; ...
            '3022201110222033'; ...
            '3022204440222033'; ...
            '3022220002222033'; ...
            '3022222222222033'; ...
            '3000000000000033'; ...
            '3333333333333333'; ...
            '3333333333333333'; ...
            '3333333333333333'];
    cm = [   0      0      0; ...  % black
             0   0.60      1; ...  % light blue
          0.53   0.53   0.53; ...  % light gray
           NaN    NaN    NaN; ...  % transparent
             0   0.73      0; ...  % light green
          0.27   0.27   0.27; ...  % gray
          0.13   0.13   0.13];     % dark gray
    cdata = ind2rgb(uint8(icon-'0'),cm);

    % If the button does not already exit
    hButton = findall(hToolbar,'Tag','ScreenCaptureButton');
    tooltip = 'Screen capture';
    if ~isempty(paramsStruct.target)
        tooltip = [tooltip ' to ' paramsStruct.target];
    end
    if isempty(hButton)
        % Add the button with the icon to the figure's toolbar
        hButton = uipushtool(hToolbar, 'CData',cdata, 'Tag','ScreenCaptureButton', 'TooltipString',tooltip, 'ClickedCallback',['screencapture(''' paramsStruct.target ''')']);  %#ok unused
    else
        % Otherwise, simply update the existing button
        set(hButton, 'CData',cdata, 'Tag','ScreenCaptureButton', 'TooltipString',tooltip, 'ClickedCallback',['screencapture(''' paramsStruct.target ''')']);
    end
function imgData = getScreenCaptureImageData(positionRect)
    if isempty(positionRect) | all(positionRect==0) | positionRect(3)<=0 | positionRect(4)<=0  %#ok ML6
        imgData = [];
    else
        % Use java.awt.Robot to take a screen-capture of the specified screen area
        rect = java.awt.Rectangle(positionRect(1), positionRect(2), positionRect(3), positionRect(4));
        robot = java.awt.Robot;
        jImage = robot.createScreenCapture(rect);

        % Convert the resulting Java image to a Matlab image
        % Adapted for a much-improved performance from:
        % http://www.mathworks.com/support/solutions/data/1-2WPAYR.html
        h = jImage.getHeight;
        w = jImage.getWidth;
        %imgData = zeros([h,w,3],'uint8');
        %pixelsData = uint8(jImage.getData.getPixels(0,0,w,h,[]));
        %for i = 1 : h
        %    base = (i-1)*w*3+1;
        %    imgData(i,1:w,:) = deal(reshape(pixelsData(base:(base+3*w-1)),3,w)');
        %end

        % Performance further improved based on feedback from Urs Schwartz:
        %pixelsData = reshape(typecast(jImage.getData.getDataStorage,'uint32'),w,h).';
        %imgData(:,:,3) = bitshift(bitand(pixelsData,256^1-1),-8*0);
        %imgData(:,:,2) = bitshift(bitand(pixelsData,256^2-1),-8*1);
        %imgData(:,:,1) = bitshift(bitand(pixelsData,256^3-1),-8*2);

        % Performance even further improved based on feedback from Jan Simon:
        pixelsData = reshape(typecast(jImage.getData.getDataStorage, 'uint8'), 4, w, h);
        imgData = cat(3, ...
            transpose(reshape(pixelsData(3, :, :), w, h)), ...
            transpose(reshape(pixelsData(2, :, :), w, h)), ...
            transpose(reshape(pixelsData(1, :, :), w, h)));
    end
function redockFigureIfRelevant(paramsStruct)
  if paramsStruct.wasDocked
      try
          set(paramsStruct.hFigure,'WindowStyle','docked');
          %drawnow;
      catch
          % never mind - ignore...
      end
  end
function imclipboard(imgData)
    % Import necessary Java classes
    import java.awt.Toolkit.*
    import java.awt.image.BufferedImage
    import java.awt.datatransfer.DataFlavor

    % Add the necessary Java class (ImageSelection) to the Java classpath
    if ~exist('ImageSelection', 'class')
        javaaddpath(fileparts(which(mfilename)), '-end');
    end
        
    % Get System Clipboard object (java.awt.Toolkit)
    cb = getDefaultToolkit.getSystemClipboard;  % can't use () in ML6!
    
    % Get image size
    ht = size(imgData, 1);
    wd = size(imgData, 2);
    
    % Convert to Blue-Green-Red format
    imgData = imgData(:, :, [3 2 1]);
    
    % Convert to 3xWxH format
    imgData = permute(imgData, [3, 2, 1]);
    
    % Append Alpha data (not used)
    imgData = cat(1, imgData, 255*ones(1, wd, ht, 'uint8'));
    
    % Create image buffer
    imBuffer = BufferedImage(wd, ht, BufferedImage.TYPE_INT_RGB);
    imBuffer.setRGB(0, 0, wd, ht, typecast(imgData(:), 'int32'), 0, wd);
    
    % Create ImageSelection object
    %    % custom java class
    imSelection = ImageSelection(imBuffer);
    
    % Set clipboard content to the image
    cb.setContents(imSelection, []);

