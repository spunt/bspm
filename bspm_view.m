function S = bspm_view(ol)
close all; 
% =========================================================================
% DEFAULTS
% =========================================================================
ul= fullfile(fileparts(which('spm.m')), 'canonical', 'single_subj_T1.nii');
if nargin==0
    fname = uigetvol('Select an Image File for Overlay', 0);
    ol = read_overlay(fname, .001, 20, 'both');
else
    if iscell(ol), ol = char(ol); end
    ol = read_overlay(ol);
end
pos     = default_positions; 
color   = default_colors; 
% =========================================================================
% CREATE THE FIGURE
% =========================================================================
try
    S.hFig    = figure(...
    'Units', 'pixels', ...
    'Position',pos.gui,...
    'Resize','on',...
    'Color',color.bg,...
    'ColorMap',gray(64),...
    'NumberTitle','off',...
    'DockControls','off',...
    'MenuBar','none',...
    'Name','bspmVIEW',...
    'CloseRequestFcn',@cb_closegui,...
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
    'DefaultUicontrolFontName','Arial',...
    'DefaultUicontrolFontSize',12,...
    'DefaultUicontrolInterruptible','on',...
    'Visible','off',...
    'Toolbar','none');
    set(S.hFig, 'ResizeFcn', @cb_resizegui); 
catch err
    rethrow(err)
end
% =========================================================================
% CREATE MENU ITEMS
% =========================================================================
try
    %% Main Menu
    S.menu1 = uimenu('Parent', S.hFig, 'Label', 'bspmVIEW');
    S.opencode  = uimenu(S.menu1, 'Label','Open GUI M-File', 'Callback', @cb_opencode); 
    S.exit      = uimenu(S.menu1, 'Label', 'Exit', 'Callback', {@cb_closegui, S});
    %% Load Menu
    S.load = uimenu(S.hFig,'Label','Load');
    S.loadol = uimenu(S.load,'Label','Overlay Image','CallBack', @cb_loadol);
    S.loadul = uimenu(S.load,'Label','Underlay Image', 'Separator', 'on', 'CallBack', @cb_loadul);
    %% Save Menu
    S.save              = uimenu(S.hFig,'Label','Save');
    S.saveintensity     = uimenu(S.save,'Label','Save Thresholded Map','CallBack', @cb_saveintensity);
    S.savemask          = uimenu(S.save,'Label','Save as Binary Mask', 'Separator', 'on', 'CallBack', @cb_savemask);
    %% Options Menu
    S.options   = uimenu(S.hFig,'Label','Options');
    S.crosshair = uimenu(S.options,'Label','Toggle Crosshairs','Checked', 'on', 'CallBack', @cb_crosshair);
catch err
    rethrow(err)
end
% =========================================================================
% INITIALIZE SPM REGISTRY & ORTHVIEWS
% =========================================================================
% try
%% REGISTRY OBJECT (HREG)
hReg = uipanel('Parent',S.hFig,'Units','Norm','Position',[0 0 1 1],...
        'BorderType', 'none', 'BackgroundColor',color.bg);
    
%% CREATE GLOBAL VARIABLE ST, PREVSECT
bspm_orthviews('Reset');
global st prevsect
prevsect    = ul;
st          = struct( ...
            'fig',          S.hFig,...
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
            'direct',       'both',...
            'snap',         []);
st.vols     = cell(24,1);
st.ol       = ol; 
spm_XYZreg('InitReg',hReg,st.ol.M,st.ol.DIM,[0;0;0]); % initialize registry object
st.ho = bspm_orthviews('Image', ul, [.025 .025 .95 .95]);
bspm_orthviews('MaxBB');
bspm_orthviews('AddBlobs', st.ho, st.ol.XYZ, st.ol.Z, st.ol.M);
bspm_orthviews('Register', hReg);
position_axes; 
position_colorbar;
addxyz; 
menu_axes;
% catch err
%     rethrow(err)
% end
% =========================================================================
% UPPER BUTTON PANEL
% =========================================================================
try
[h,axpos] = get_axes_handles;
uppos = axpos(2,:);
uppos(2) = sum(uppos([2 4])) + .01; 
uppos(3) = 1 - 2*uppos(1);
uppos(4) = 1 - uppos(2) - .01; 
S.uppan = uipanel('Parent',S.hFig,'Units','Norm',...
    'Position',uppos,...
    'BorderType', 'none', ...
    'BackgroundColor',color.bg);
S.contrasttitle = uicontrol('Parent', S.uppan, 'Units', 'Normal', 'FontUnits', 'Normal', ...
        'Style', 'Text', ...
        'Position', [0 .1 .6 .8], ...
        'BackgroundColor', color.bg, ...
        'ForegroundColor', color.fg, ...
        'Horizontal', 'Left', ...
        'Clip', 'off', ...
        'String', st.ol.descrip, ...
        'FontUnits', 'norm', ...
        'FontName', 'Arial', ...
        'FontSize', .45, ...
        'Enable', 'Inactive', ...
        'Tag', 'ContrastName', ...
        'FontWeight', 'Normal');
RADIO   = {'parent',S.uppan, 'style','radio', 'units','norm', 'visible','on', 'backg', color.bg, 'fontn', 'fixed-width', 'foreg', color.fg, 'fonts',20, 'pos'};
S.pos1 = uicontrol(RADIO{:}, [.620 .1 .10 .80], 'tag', 'direct', 'str', 'pos', 'callback', @cb_directmenu);
S.pos2 = uicontrol(RADIO{:}, [.730 .1 .10 .80], 'tag', 'direct', 'str', 'neg', 'callback', @cb_directmenu); 
S.pos3 = uicontrol(RADIO{:}, [.840 .1 .15 .80], 'tag', 'direct', 'str', 'pos/neg', 'value', 1, 'enable', 'inactive', 'callback', @cb_directmenu); 

catch err
    rethrow(err)
end
% =========================================================================
% LOWER BUTTON PANE
% =========================================================================
try
lowpos = axpos(1,:);
lowpos(1) = axpos(3, 1); 
lowpos(3) = 1 - lowpos(1) - axpos(1,2);
S.lowpan = uipanel('Parent',S.hFig, 'Units','Norm', 'Position',lowpos,...
        'BorderType','none', 'BackgroundColor',color.bg);
EDIT    = {'parent',S.lowpan, 'style','edit', 'units','norm', 'clip','off', 'visible','on', 'backg','w', 'fonts',16, 'pos'};
TEXT    = {'parent',S.lowpan, 'style','text', 'units','norm', 'clip','off', 'visible','on', 'backg', color.bg, 'foreg', color.fg, 'fonts',14, 'fontw','bold', 'pos'};
POPMENU = {'parent',S.lowpan, 'style','popupmenu', 'units','norm', 'fontunits', 'norm', 'visible','on', 'fonts',.30, 'pos'};
SLIDER  = {'parent',S.lowpan, 'style','slider', 'units', 'norm', 'visible','on', 'min', 1.0000e-20, 'max', 1, 'sliderstep', [.0001 .001], 'value', st.ol.P, 'pos'};
Tpos = {pos.k pos.tval pos.pval pos.df};
Tstr = {'Extent' 'Thresh' 'P-Value' 'DF'};
Tdefvalues = [st.ol.K st.ol.U st.ol.P st.ol.DF];
Tstrform = {'%d' '%2.3f' '%d' '%d'}; 
Ecallback   = {@cb_updateoverlay, @cb_updateoverlay, @cb_updateoverlay, @cb_updateoverlay};
for i = 1:length(Tstr)
    txpos = Tpos{i};
    txpos(2) = sum(txpos([2 4]));
    txpos(4) = txpos(4)*.50; 
    S.tx(i) = uicontrol(TEXT{:}, txpos, 'string', Tstr{i});
    S.ed(i) = uicontrol(EDIT{:}, Tpos{i}, 'string', sprintf(Tstrform{i}, Tdefvalues(i)), 'Tag', Tstr{i}, 'callback', Ecallback{i}); 
end
% S.pslider       = uicontrol(SLIDER{:}, pos.pslider, 'callback', @cb_pslider);

catch err
    rethrow(err)
end
% =========================================================================
% *
% * CALLBACK FUNCTIONS
% *
% =========================================================================
function cb_minmax(varargin)
global st
lab = get(varargin{1}, 'label');
if regexp(lab, 'global max')
    centre = st.ol.XYZmm(:,st.ol.Z==max(st.ol.Z));
elseif regexp(lab, 'global min')
    centre = st.ol.XYZmm(:,st.ol.Z==min(st.ol.Z)); 
end
bspm_orthviews('reposition', centre); 
function cb_directmenu(varargin)
    global st
    cbh     = varargin{1};
    str     = get(varargin{1}, 'string');
    allh = findobj(st.fig, 'Tag', 'direct'); 
    allhstr = get(allh, 'String');
    set(allh(strcmp(allhstr, str)), 'Value', 1, 'Enable', 'inactive'); 
    set(allh(~strcmp(allhstr, str)), 'Value', 0, 'Enable', 'on');
    T = getthreshinfo;
    if strcmp(str, 'pos/neg')
        changecolormap(jet(64)); 
        st.ol = read_overlay(st.ol.fname, T.pval, T.extent, 'both'); 
    else
        changecolormap(hot(64)); 
        st.ol = read_overlay(st.ol.fname, T.pval, T.extent, str);
    end
    bspm_orthviews('RemoveBlobs', st.ho); 
    bspm_orthviews('MaxBB');
    bspm_orthviews('AddBlobs', st.ho, st.ol.XYZ, st.ol.Z, st.ol.M);
    bspm_orthviews('Register', st.registry.hReg);
    position_axes;
    position_colorbar; 
function cb_loadol(varargin)
    global st
    fname = uigetvol('Select an Image File for Overlay', 0);
    T = getthreshinfo; 
    if strcmpi(T.direct, 'Show +'), direc = 'pos'; 
    elseif strcmpi(T.direct, 'Show -'), direc = 'neg';
    else direc = 'both'; end
    st.ol = read_overlay(fname, T.pval, T.extent, direc);
    bspm_orthviews('RemoveBlobs', st.ho); 
    bspm_orthviews('MaxBB');
    bspm_orthviews('AddBlobs', st.ho, st.ol.XYZ, st.ol.Z, st.ol.M);
    bspm_orthviews('Register', st.registry.hReg);
    position_axes;
    updatecontrastname; 
    position_colorbar;
function cb_loadul(varargin)
    global st
    ul = uigetvol('Select an Image File for Underlay', 0);
    bspm_orthviews('Delete', st.ho); 
    st.ho = bspm_orthviews('Image', ul, [.025 .025 .95 .95]);
    bspm_orthviews('MaxBB');
    bspm_orthviews('AddBlobs', st.ho, st.ol.XYZ, st.ol.Z, st.ol.M);
    bspm_orthviews('Register', st.registry.hReg);
    position_axes; 
    position_colorbar;
    addxyz;
    menu_axes;
function cb_opencode(varargin)
    open(mfilename('fullpath'));
function cb_closegui(varargin)
    if length(varargin)==3
        parent = varargin{3};
        h = parent.hFig;
    else
        h = varargin{1};
    end
    delete(h); % Bye-bye figure
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
    addxyz; 
function cb_updateoverlay(varargin)   
global st
T = getthreshinfo;
tag = get(varargin{1}, 'tag');
switch tag
    case {'Thresh'}
        T.pval = bob_t2p(T.thresh, T.df);
    case {'P-Value'}
        T.thresh = bob_p2t(T.pval, T.df); 
    case {'DF'}
        T.thresh = bob_p2t(T.pval, T.df); 
        T.pval = bob_t2p(T.thresh, T.df);
    case {'Extent'}
end
updatethreshinfo(T);
[X,Y,Z]     = ndgrid(1:st.ol.DIM(1),1:st.ol.DIM(2),1:st.ol.DIM(3));
st.ol.XYZ   = [X(:)';Y(:)';Z(:)'];
RCP         = st.ol.XYZ; 
RCP(4,:)    = 1;
st.ol.XYZmm       = st.ol.M(1:3,:)*RCP;
if strcmpi('both', T.direct)
    st.ol.XYZmm       = st.ol.XYZmm(:,abs(st.ol.Y(:))>=T.thresh);
    st.ol.XYZ         = st.ol.XYZ(:,abs(st.ol.Y(:))>=T.thresh);
    st.ol.Z           = st.ol.Y(abs(st.ol.Y(:))>=T.thresh);
else
    st.ol.XYZmm       = st.ol.XYZmm(:,st.ol.Y(:)>=T.thresh);
    st.ol.XYZ         = st.ol.XYZ(:,st.ol.Y(:)>=T.thresh);
    st.ol.Z           = st.ol.Y(st.ol.Y(:)>=T.thresh);
end
cl_index   = spm_clusters(st.ol.XYZ);
cl_bin     = repmat(1:max(cl_index), length(cl_index), 1)==repmat(cl_index', 1, max(cl_index));
nopassidx  = ismember(cl_index, find(sum(cl_bin) < T.extent));  
st.ol.Z(nopassidx,:) = []; 
st.ol.XYZ(:,nopassidx) = []; 
st.ol.XYZmm(:,nopassidx) = []; 
bspm_orthviews('RemoveBlobs', st.ho);
bspm_orthviews('MaxBB');
bspm_orthviews('AddBlobs', st.ho, st.ol.XYZ, st.ol.Z, st.ol.M);
position_axes;
updatecontrastname; 
position_colorbar;
% =========================================================================
% *
% * SUBFUNCTIONS
% *
% =========================================================================
function OL = read_overlay(fname, pval, k, direct)
    if nargin<4, direct = 'both'; end
    if nargin<3, k = 20; end
    if nargin<2, pval = .001; end
    [od, oh] = bspm_read_vol(fname);
    tmp = oh.descrip;
    idx1 = regexp(tmp,'[','ONCE');
    idx2 = regexp(tmp,']','ONCE');
    df = str2num(tmp(idx1+1:idx2-1));
    u = bob_p2t(pval, df);
    if strcmpi('neg', direct), od = od*-1; end
    
    M           = oh.mat;         %-voxels to mm matrix
    DIM         = oh.dim';
    VOX         = abs(diag(M(:,1:3))); 
    [X,Y,Z]     = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
    XYZ         = [X(:)';Y(:)';Z(:)'];
    RCP         = XYZ; 
    RCP(4,:)    = 1;
    XYZmm       = M(1:3,:)*RCP;
    if strcmpi('both', direct)
        XYZmm       = XYZmm(:,abs(od(:))>=u);
        XYZ         = XYZ(:,abs(od(:))>=u);
        Z           = od(abs(od(:))>=u);
    else
        XYZmm       = XYZmm(:,od(:)>=u);
        XYZ         = XYZ(:,od(:)>=u);
        Z           = od(od(:)>=u);
    end
    cl_index    = spm_clusters(XYZ);
    cl_bin      = repmat(1:max(cl_index), length(cl_index), 1)==repmat(cl_index', 1, max(cl_index));
    nopassidx  = ismember(cl_index, find(sum(cl_bin) < k));  
    Z(nopassidx,:) = []; 
    XYZ(:,nopassidx) = []; 
    XYZmm(:,nopassidx) = []; 
    OL          = struct( ...
                'fname',    fname,...
                'descrip',  oh.descrip, ...
                'DF',       df, ...
                'U',        u, ...
                'P',        pval, ...
                'K',        k, ...
                'Y',        od, ...
                'M',        M,...
                'Z',        Z,...
                'XYZmm',    XYZmm,...
                'XYZ',      XYZ,...
                'DIM',      DIM,...
                'VOX',      VOX);   
function T = getthreshinfo
global st
T.extent = str2num(get(findobj(st.fig, 'Tag', 'Extent'), 'String')); 
T.thresh = str2num(get(findobj(st.fig, 'Tag', 'Thresh'), 'String'));
T.pval = str2num(get(findobj(st.fig, 'Tag', 'P-Value'), 'String'));
T.df = str2num(get(findobj(st.fig, 'Tag', 'DF'), 'String'));
tmph = findobj(st.fig, 'Tag', 'direct'); 
opt = get(tmph, 'String');
T.direct = opt(find(cell2mat(get(tmph, 'Value'))));
if strcmp(T.direct, 'pos/neg'), T.direct = 'both'; end
function updatethreshinfo(T)
global st
Tval = [T.extent T.thresh T.pval T.df]; 
Tstr = {'Extent' 'Thresh' 'P-Value' 'DF'};
Tstrform = {'%d' '%2.3f' '%d' '%d'}; 
for i = 1:length(Tstr)
    set(findobj(st.fig, 'Tag', Tstr{i}), 'String', sprintf(Tstrform{i}, Tval(i)));
end
function updatecontrastname
global st
connamh = findobj(st.fig, 'Tag', 'ContrastName'); 
set(connamh, 'String', st.ol.descrip); 
function changecolormap(newmap)
global st
cmap = [gray(64); newmap]; 
set(st.fig,'Colormap', cmap);
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
function [h, axpos] = get_axes_handles(varargin)
    global st
    axpos = zeros(3,4);
    for a = 1:3
        tmp = st.vols{1}.ax{a};
        h.ax(a) = tmp.ax; 
        h.d(a)  = tmp.d;
        h.lx(a) = tmp.lx; 
        h.ly(a) = tmp.ly;
        axpos(a,:) = get(h.ax(a), 'position');
    end
function position_axes
    %% Handles for axes
    % 1 - transverse
    % 2 - coronal
    % 3 - sagittal 
    % st.vols{1}.ax{1}.ax   - axes
    % st.vols{1}.ax{1}.d    - image
    % st.vols{1}.ax{1}.lx   - crosshair (x)
    % st.vols{1}.ax{1}.ly   - crosshair (y)
    [h,axpos] = get_axes_handles;
    MARG    = .01; 
    RAT     = 1.1990; 
    SZ1     = .405; 
    SZ2     = SZ1*RAT; 
    axpos(1:2,1) = MARG; 
    axpos(1, 2)  = MARG; 
    axpos(:, 3) = [SZ1 SZ1 SZ2]; 
    axpos(:, 4) = [SZ2 SZ1 SZ1];
    axpos(3, 1) = SZ1 + MARG*2; 
    axpos(2:3,2) = SZ2 + MARG*2;   
    for a = 1:3, set(h.ax(a), 'position', axpos(a,:)); end
    bspm_orthviews('Redraw');
function menu_axes
    [h,axpos] = get_axes_handles;
    cmenu = uicontextmenu;
    ctmax = uimenu(cmenu, 'Label', 'Go to global max', 'callback', @cb_minmax);
    ctmin = uimenu(cmenu, 'Label', 'Go to global min', 'callback', @cb_minmax);
    ctsavemap = uimenu(cmenu, 'Label', 'Save map', 'separator', 'on');
    ctsavemask = uimenu(cmenu, 'Label', 'Save mask');
    for a = 1:3
        set(h.ax(a), 'uicontextmenu', cmenu); 
    end
function addxyz
    global st
    h = get_axes_handles;
    xyz = round(spm_XYZreg('GetCoords',st.registry.hReg));
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
function position_colorbar
global st
[h,axpos] = get_axes_handles;

%% Handle for colorbar
cbh = st.vols{1}.blobs{1}.cbar; 
cbpos = axpos(3,:); 
cbpos(4) = cbpos(4)*.9; 
cbpos(2) = cbpos(2) + (axpos(3,4)-cbpos(4))/2; 
cbpos(1) = sum(cbpos([1 3])); 
cbpos(3) = (1 - cbpos(1))/2; 
cbpos(1) = cbpos(1) + (cbpos(3)/4); 
yl = get(st.vols{1}.blobs{1}.cbar, 'ylim');

set(st.vols{1}.blobs{1}.cbar, 'ycolor', st.color.fg, 'fontsize', 12, ...
    'ytick', [ceil(min(yl)) round(median(yl)) floor(max(yl))], ...
    'position', cbpos, 'YAxisLocation', 'right'); 

%% Handle for crosshairs
set(h.lx, 'color', st.color.xhair); 
set(h.ly, 'color', st.color.xhair);
function color = default_colors 
    color.bg        = [0/255 0/255 0/255];
    color.fg        = [248/255 248/255 248/255];
    color.border    = [023/255 024/255 020/255]*2;
    color.xhair     = [0.7020    0.8039    0.8902];
function pos = default_positions 
    screensize      = get(0, 'ScreenSize');
    pos.ss          = screensize(3:4);
    pos.gui         = [pos.ss(1:2)*.5 pos.ss(2)*.55 pos.ss(2)*.5];
    pos.aspratio    = pos.gui(3)/pos.gui(4);
    %% BUTTONS
%     bpos = position_row(.2, .125, .75, 4); 
%     pos.k           = bpos(1,:);
%     pos.tval        = bpos(2,:);
%     pos.pval        = bpos(3,:);
%     pos.df          = bpos(4,:);
    pos.k              = [.025 .200 .200 .125];
    pos.tval           = [.250 .200 .225 .125];
    pos.pval           = [.500 .200 .275 .125];
    pos.df             = [.800 .200 .175 .125];
    pos.pslider        = [.050 .025 .900 .150];
function OL = thresh_overlay(in, u, k, direct)
    % THRESH_OVERLAY
    %
    % USAGE: out = thresh_overlay(in, u, k)
    %
    %   ARGUMENTS
    %       in:     3D matrix to threshold
    %       u:      height threshold
    %       k:      extent threshold
    %       direct: direction, 'pos', 'neg', or 'both' (default = 'both')
    %       
   
    if strcmpi('neg', direct), in = in*-1; end
    imdims = size(in);
    global st
    imdims = size(in);
    % if necessary, calculate critical t
    if ismember(u,[.10 .05 .01 .005 .001 .0005 .0001]);
        tmp = in_hdr.descrip;
        idx1 = regexp(tmp,'[','ONCE');
        idx2 = regexp(tmp,']','ONCE');
        df = str2num(tmp(idx1+1:idx2-1));
        u = bob_p2t(u, df);
    end
    in(in<u) = NaN;
    in(in==0)=NaN;

    % grab voxels
    % ------------------------------------------------------
    [X Y Z] = ind2sub(size(in), find(in > 0));
    voxels = sortrows([X Y Z])';

    % get cluster indices of voxels
    % ------------------------------------------------------
    cl_index = spm_clusters(voxels);

    % find index of clusters of sufficient size
    % ------------------------------------------------------
    for i = 1:max(cl_index)
        a(cl_index == i) = sum(cl_index == i);
    end
    which_vox = (a >= k);
    cl_vox = voxels(:,which_vox);
    cl_vox = cl_vox';
    roi_mask = zeros(imdims);
    for i = 1:size(cl_vox,1)
        roi_mask(cl_vox(i,1),cl_vox(i,2),cl_vox(i,3)) = in(cl_vox(i,1),cl_vox(i,2),cl_vox(i,3));
    end
    out = double(roi_mask);
    out(out==0) = NaN;
    st.ol.tY    = out;
    st.ol.Z     = out(~isnan(out(:)));
    st.ol.XYZ   = cl_vox; 
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
function matlab_fslview(over, under)
% matlab_fslview Call fslview from MATLAB
% fslview [-m 3d|ortho|lightbox] <baseimage> [-l lutname] [-b low,hi]
% 	[ <overlay> [-l lutname] [-b low,hi] ] ...
% fslview -m ortho,lightbox filtered_func_data thresh_zstat1 -t 0.5 thresh_zstat2 -l "Cool" -t 0.5
% 
% Optional arguments (You may optionally specify one or more of):
% 	-V,--verbose	switch on diagnostic messages
% 	-h,--help	display this message
% 	-m,--mode	Initial viewer mode. Comma separated list of: 3d; single, ortho; lightbox
% 
% 
% Per-image options
% 
% Usage:
% image [-l GreyScale] [-t 0.1] [-b 2.3,6]
% 	-l,--lut	Lookup table name. As per GUI, one of: Greyscale;
% 			"Red-Yellow"; "Blue-Lightblue"; Red; Green;
% 			Blue; Yellow; Pink; Hot; Cool; Copper, etc.
% 	-b,--bricon	Initial bricon range, e.g., 2.3,6
% 	-t,--trans	Initial transparency, e.g., 0.2
        % -----------------------------------------------------
    if iscell(over), over = char(over); end
    if iscell(under), under = char(under); end
    htmp = spm_vol(under); dtmp = spm_read_vols(htmp);
    dtmp = dtmp(:); 
    dtmp(dtmp < nanmean(dtmp)/10) = [];
    command = sprintf('fslview -m ortho %s -l "Greyscale" -b %2.3f,%2.3f %s -l "Blue" -t .15 &', under, min(dtmp), max(dtmp), over);
    system(command);
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
function bpos = position_row(buttonw, buttonh, lowerh, nbutton)
    bm = .02;
    for i = 1:nbutton
        bpos(i,:) = [(i*bm)+(i-1)*buttonw lowerh buttonw buttonh]; 
    end
    lm = bpos(1,1);
    rm = 1 - sum(bpos(nbutton,[1 3])); 
    bpos(:,1) = bpos(:,1) + (rm+lm)*.75; 
% =========================================================================
% *
% * SUBFUNCTIONS (MIP)
% *
% =========================================================================
function addmip(varargin)
%% MIP

% Add MIP to Panel
hMIPax  = axes('Parent',hReg,'Position',[.05 .05 .9 .9],'Visible','off');

%% MIP
hMIPax  = bspm_mip_ui(Z,XYZmm,M,DIM,hMIPax,units);

%% Cross Register the MIP with Registry, Push Coords
spm_XYZreg('XReg',hReg, hMIPax,'bspm_mip_ui');
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

    xyz = spm_XYZreg('RoundCoords',[0;0;0],M,DIM);

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
            [xyz,d] = spm_XYZreg('RoundCoords',xyz,MD.M,MD.DIM);
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
        if ~isempty(MD.hReg) && MD.hReg~=hC, spm_XYZreg('SetCoords',xyz,MD.hReg,h); end

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
                [xyz,i,d] = spm_XYZreg('NearestXYZ',oxyz,MD.XYZ);
            case 'nrmax'
                str       = 'nearest local maximum';
                iM        = inv(MD.M);
                XYZvox    = iM(1:3,:)*[MD.XYZ; ones(1,size(MD.XYZ,2))];
                [null,null,XYZvox] = spm_max(MD.Z,XYZvox);
                XYZ       = MD.M(1:3,:)*[XYZvox; ones(1,size(XYZvox,2))];
                [xyz,i,d] = spm_XYZreg('NearestXYZ',oxyz,XYZ);
            case 'glmax'
                str       = 'global maximum';
                [null, i] = max(MD.Z); i = i(1);
                xyz       = MD.XYZ(:,i);
                d         = sqrt(sum((oxyz-xyz).^2));
            case 'nrchan'
                str       = 'nearest suprathreshold channel';
                if ~isfield(MD, 'hChanPlot'), bspm_mip_ui('Channels',h); end
                MD        = get(h,'UserData');
                [xyz,i,d] = spm_XYZreg('NearestXYZ',[oxyz(1); oxyz(2); 0],MD.Channels.pos);
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
            xyz    = spm_XYZreg('RoundCoords',xyz,MS.M,MS.DIM);
        elseif DragType==1
            xyz    = spm_XYZreg('RoundCoords',xyz,MS.M,MS.DIM);
            hMIPax = get(cO,'Parent');
            MD     = get(hMIPax,'UserData');
            i      = spm_XYZreg('FindXYZ',xyz,MD.XYZ);
        elseif DragType==2
            hMIPax = get(cO,'Parent');
            MD     = get(hMIPax,'UserData');
            [xyz,i,d] = spm_XYZreg('NearestXYZ',xyz,MD.XYZ);
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
            spm_XYZreg('SetCoords',get(MS.hMIPxyz,'UserData'),MS.hReg,hMIPax);
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
        [XYZmm,i] = spm_XYZreg('NearestXYZ',XYZmm,xSPM.XYZmm);
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
% =========================================================================
% *
% * SUBFUNCTIONS (SPM_ORTHVIEWS)
% *
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
            spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
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
            xyz = spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
        end
        cm_pos;
        xyzstr = num2str([-99; round(xyz)]); 
        xyzstr(1,:) = [];
        
        axidx = [3 2 1];
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
function H = specify_image(img)
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
    ax = axes('Visible','off', 'Parent', st.fig, ...
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
cbpos(1) = cbpos(1) + (cbpos(3)/4); 
yl = [st.vols{vh}.blobs{bh}.min st.vols{vh}.blobs{bh}.max]; 
st.vols{vh}.blobs{bh}.cbar = axes('Parent',st.fig, 'ycolor', [1 1 1], ...
    'position', cbpos, 'YAxisLocation', 'right', ...
    'ytick', [ceil(min(yl)) round(median(yl)) floor(max(yl))], ...
    'Box','off', 'YDir','normal', 'XTickLabel',[], 'XTick',[]); 
if isfield(st.vols{vh}.blobs{bh},'name')
    ylabel(st.vols{vh}.blobs{bh}.name,'parent',st.vols{vh}.blobs{bh}.cbar);
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
yl = [st.vols{vh}.blobs{bh}.min st.vols{vh}.blobs{bh}.max]; 
image([0 1],interval,cdata,'Parent',st.vols{vh}.blobs{bh}.cbar);
set(st.vols{vh}.blobs{bh}.cbar, 'ycolor', [1 1 1], ...
    'position', cbpos, 'YAxisLocation', 'right', ...
    'ytick', [ceil(min(yl)) round(median(yl)) floor(max(yl))], ...
    'Box','off', 'YDir','normal', 'XTickLabel',[], 'XTick',[]); 
if isfield(st.vols{vh}.blobs{bh},'name')
    ylabel(st.vols{vh}.blobs{bh}.name,'parent',st.vols{vh}.blobs{bh}.cbar);
end
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
function register(hreg)
global st
%tmp = uicontrol('Position',[0 0 1 1],'Visible','off','Parent',st.fig);
h   = valid_handles;
if ~isempty(h)
    tmp = st.vols{h(1)}.ax{1}.ax;
    st.registry = struct('hReg',hreg,'hMe', tmp);
    spm_XYZreg('Add2Reg',st.registry.hReg,st.registry.hMe, 'bspm_orthviews');
else
    warning('Nothing to register with');
end
st.centre = spm_XYZreg('GetCoords',st.registry.hReg);
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
function H = pos(handle)
global st
H = [];
for i=valid_handles(handle)
    is = inv(st.vols{i}.premul*st.vols{i}.mat);
    H = is(1:3,1:3)*st.centre(:) + is(1:3,4);
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
                    changecolormap(jet(64));
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
function img = scaletocmap(inpimg,mn,mx,cmap,miscol)
if nargin < 5, miscol=1; end
cml = size(cmap,1);
scf = (cml-1)/(mx-mn);
img = round((inpimg-mn)*scf)+1;
img(img<1)   = 1;
img(img>cml) = cml;
img(~isfinite(img)) = miscol;
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
function addcontexts(handles)
for ii = valid_handles(handles)
    addcontext(ii);
end
bspm_orthviews('reposition',bspm_orthviews('pos'));
function rmcontexts(handles)
global st
for ii = valid_handles(handles)
    for i=1:3
        set(st.vols{ii}.ax{i}.ax,'UIcontextmenu',[]);
        try, st.vols{ii}.ax{i} = rmfield(st.vols{ii}.ax{i},'cm'); end
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
function cm_handles = get_cm_handles
global st
cm_handles = [];
for i = valid_handles
    cm_handles = [cm_handles st.vols{i}.ax{1}.cm];
end
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
function m = max_img
m = 24;









    
   


