function a = getset(h)
% GETSET Current and possible values of object properties
%
%  USAGE: getset(H, varargin)
% ________________________________________________________________________________________
%  INPUTS
%	 H: object   
%

% ----------------------------- Copyright (C) 2015 Bob Spunt -----------------------------
%	Created:  2015-12-26
%	Email:     spunt@caltech.edu
% ________________________________________________________________________________________
if nargin < 1, mfile_showhelp; return; end
if ~ishandle(h), disp('Input is not a valid MATLAB object! Doing nothing...'); a = []; return; end
colorb      = [23/255 24/255 20/255];
colorf      = [1 1 1];
[sn, sv]    = divvystruct(set(h), 1);
[gn, gv]    = divvystruct(get(h), 1);
tmpidx      = ~cellfun('isempty', regexp(sn, '.*Fcn'));
sn(tmpidx)  = [];
sv(tmpidx)  = [];
exclude     = props2exclude(h)
tmpidx      = find(ismember(sn, exclude));
sn(tmpidx)  = []; 
sv(tmpidx)  = [];
ln          = cellfun('length', sv);
nidx        = find(ismember(gn, sn));
gn          = gn(nidx); 
gv          = gv(nidx);
for i = 1:length(sv)
    cv = sv{i};
    if isempty(cv)
        sv{i} = gv{i};
    else
        cv(strcmp(cv, gv{i})) = [];
        sv{i} = vertcat(gv(i), cv);   
    end
end
svclass     = cellfun(@class, gv, 'unif', false);
valididx    = find(ismember(svclass, {'double' 'char' 'logical'}));
svclass     = svclass(valididx);
gv          = gv(valididx);
gn          = gn(valididx);
sv          = sv(valididx);
sn          = sn(valididx);
ln          = ln(valididx);
nrow        = length(sv);
svstyle     = repmat({'edit'}, nrow, 1);
svstyle(ln > 1) = {'popup'};
a           = cell2struct(sv, sn);
a           = orderfields(a);
% | FIGURE
hfig        = putfigure(sn, h, colorb);
% | PANEL GRID
hgrid = pgrid(nrow, 2       , ...
'parent',       hfig        , ...
'foreg',        colorf      , ...
'backg',        colorb      , ...
'marginsep',    .025        , ...
'panelsep',     .010        , ...
'relwidth',     [3 4]         ...
);
pn          = hgrid(1:2:end);
pv          = hgrid(2:2:end);
hnprop      = {'style', 'text', 'units', 'norm', 'pos', [.025 .025 .95 .95], 'fontsize', 12, 'fontweight', 'bold', 'foreg', colorf, 'backg', colorb, 'horiz', 'right'}; 
hvprop      = {'units', 'norm', 'pos', [.01 .025 .98 .95], 'fontsize', 10, 'foreg', [0 0 0], 'backg', [1 1 1]};
for i = 1:nrow 
    lab   = sn{i};
    val   = sv{i};
    if ismember(svclass{i}, {'double' 'logical'}), val = propnum2str(val); end
    hn(i) = uibutton(pn(i), hnprop{:}, 'string', lab);
    hv(i) = uicontrol(pv(i), 'style', svstyle{i}, hvprop{:}, 'value', 1, 'string', val, 'tag', sn{i}, 'Callback', {@cb_set, h, svclass{i}});
    drawnow; 
end
set(findall(hfig, '-property', 'Units'), 'Units', 'Norm'); drawnow; 
set(findall(hfig, '-property', 'FontUnits'), 'FontUnits', 'Norm'); drawnow; 
set(hfig, 'visible', 'on');
end
% ========================================================================================
%
% ------------------------------------- SUBFUNCTIONS -------------------------------------
%
% ========================================================================================
function cb_set(varargin)
opt     = get(varargin{1}, 'string');
prop    = get(varargin{1}, 'tag'); 
try
    if iscell(opt)
        set(varargin{3}, prop, opt{get(varargin{1}, 'value')});
    else
        if strcmpi(varargin{4}, 'double')
            opt = propstr2num(opt);
        else
            opt = propstr2cellstr(opt);
        end
        set(varargin{3}, prop, opt);
    end
catch
    opt = get(varargin{3}, prop);
    if isnumeric(opt), opt = propnum2str(opt); end
    set(varargin{1}, 'String', opt);
end
if strcmpi(prop, 'units')
    p   = get(get(varargin{1}, 'parent'), 'parent');
    ph  = findall(p, 'tag', 'Position');
    set(ph, 'String', propnum2str(get(varargin{3}, 'Position'))); 
elseif strcmpi(prop, 'fontunits')
    p   = get(get(varargin{1}, 'parent'), 'parent');
    ph  = findall(p, 'tag', 'FontSize');
    set(ph, 'String', propnum2str(get(varargin{3}, 'FontSize'))); 
end

end
function mfile_showhelp(varargin)
% MFILE_SHOWHELP
ST = dbstack('-completenames');
if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
eval(sprintf('help %s', ST(2).file));
end
function exclude    = props2exclude(h)
exclude = {
    'ALim'
    'ALimMode'
    'CLim'
    'CLimMode'
    'FontSmoothing'
    'ColorOrder'
    'ColorOrderIndex'
    'SortMethod'
    'Selected'
    'YMinorGrid'
    'XMinorGrid'
    'YMinorTick'
    'XMinorTick'
    'SelectionHighlight'
    'PickableParts'
    'XScale'
    'YScale'
    'ActivePositionProperty'
    'AmbientLightColor'
    'Alphamap'
    'Callback'
    'DataAspectRatio'
    'DataAspectRatioMode'
    'BusyAction'
    'Clipping'
    'GridAlpha'
    'CameraPosition'
    'CameraPositionMode'
    'CameraTarget'
    'CameraTargetMode'
    'CameraUpVector'
    'CameraUpVectorMode'
    'CameraViewAngle'
    'CameraViewAngleMode'
    'PlotBoxAspectRatio'
    'PlotBoxAspectRatioMode'
    'Projection'
    'ZAxis'
    'View'
    'TickLabelInterpreter'
    'ZColor'
    'ZColorMode'
    'Layer'
    'LineStyleOrder'
    'LineStyleOrderIndex'
    'ZDir'
    'ZGrid'
    'ZLabel'
    'ZLim'
    'ZLimMode'
    'ZMinorGrid'
    'ZMinorTick'
    'ZScale'
    'ZTick'
    'ZTickLabel'
    'ZTickLabelMode'
    'XTickLabelMode'
    'YTickLabelMode'
    'ZTickLabelRotation'
    'TickDirMode'
    'ZTickMode'
    'XColorMode'
    'YColorMode'
    'YLimMode'
    'XLimMode' 
    'GridAlphaMode'
    'GridColor'
    'GridColorMode'
    'GridLineStyle'
    'MinorGridAlpha'
    'MinorGridAlphaMode'
    'MinorGridColor'
    'MinorGridColorMode'
    'MinorGridLineStyle'
    'NextPlot'
    'HandleVisibility'
    'HitTest'
    'Interruptible'
    'ClippingStyle'
    'CData'
    'CellEditCallback'
    'CellSelectionCallback'
    'Children'
    'Colormap'
    'ColumnFormat'
    'ColumnName'
    'CurrentAxes'
    'CurrentObject'
    'Data'
    'Parent'
    'PointerShapeCData'
    'UIContextMenu'
    'UserData'
};
end
function hfig       = putfigure(sn, h, colorb)
padh    = 10;
padw    = 30; 
strln   = cellfun('length', sn); 
strtmp  = sn(strln==max(strln));
[sw,sh] = strsize(strtmp{1});
ss      = get(0, 'ScreenSize'); 
fw      = (sw+padw)*3; 
fh      = (sh+padh)*(length(sn));
hfig  = figure(...
        'Name', sprintf('getset: %s', get(h, 'Type')), ...
        'Units', 'pixels', ...
        'Position', [1 1 fw fh],...
        'Resize','on',...
        'Color', colorb, ...
        'NumberTitle','off',...
        'DockControls','off',...
        'MenuBar','none',...
        'Visible','off',...
        'Toolbar','none');
end
function val        = propnum2str(val)
    if size(val, 1)==1
        val = num2str(val); 
    else 
        tmp = []; 
        for r = 1:size(val, 1)
            tmp = [tmp num2str(val(r,:)) '; ']; 
        end
        tmp         = strtrim(tmp);
        tmp(end)    = []; 
        val         = tmp;
    end
    val = regexprep(val, '\s+', ' ');
end
function val        = propstr2num(val)
    if regexp(val, ';')
        tmp = strtrim(regexp(val, ';', 'split'));
        tmp = cellfun(@str2num, tmp, 'Unif', false);
        val = vertcat(tmp{:});
    else
        val = str2num(val);
    end
end
function val        = propstr2cellstr(val)
    if regexp(val, ';')
        val = strtrim(regexp(val, ';', 'split'));
    end
end
function [fn, fv]   = divvystruct(s, sortflag)
fn = fieldnames(s);
fv = struct2cell(s);
if sortflag
    tmp = [fn fv];
    tmp = sortrows(tmp, 1);
    fn = tmp(:,1);
    fv = tmp(:,2);
end
end
function [w, h]     = strsize(string, varargin)
% STRSIZE Calculate size of string
%
%  USAGE: strsize(string, varargin) 
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-07-14
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
def = { ... 
    'FontSize',         get(0, 'DefaultTextFontSize'),      ...
    'FontName',         get(0, 'DefaultTextFontname'),      ...
    'FontWeight',       get(0, 'DefaultTextFontWeight'),    ...
    'FontAngle',        get(0, 'DefaultTextFontAngle'),     ...
    'FontUnits',        get(0, 'DefaultTextFontUnits')     ...
	};
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end

% | Get text size in data units
tmpfig  = figure('visible', 'off');
hTest   = text(1,1, string, 'Units','Pixels', 'FontUnits',FontUnits,...
    'FontAngle',FontAngle,'FontName',FontName,'FontSize',FontSize,...
    'FontWeight',FontWeight,'Parent',gca, 'Visible', 'off');
textExt = get(hTest,'Extent');
delete(tmpfig);
h    = textExt(4);
w    = textExt(3);

% | If using a proportional font, shrink text width by a fudge factor to account for kerning.
if ~strcmpi(FontName, 'Fixed-Width'), w = w*0.9; end 
end