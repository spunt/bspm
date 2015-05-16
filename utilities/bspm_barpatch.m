function h = bspm_barpatch(data, varargin)
% BSPM_BARPATCH Create bar graph with error bars using patch and line objects
% 
% This function will create a grouped bar graph with error bars without
% using the standard plotting functions BAR and ERRORBAR. It uses PATCH to
% create the bars and LINE to construct the error bars.
% 
%  USAGE: h = bspm_barpatch(data, varargin)
% __________________________________________________________________________
%  OUTPUT
% 	h: handles to all graphics objects  
% __________________________________________________________________________
%  INPUTS
%   data:           data matrix to plot; rows are cases, cols are variables  
%   varargin:       optional arguments entered as "name,value" pairs: 
%       figh        - handle for figure to plot in
%       groupidx    - rows index columns of "data" to plot as a group
%       groupname   - labels for different groups of bars
%       grouptick   - flag to place tickmark between groups on x-axis
%       barname     - labels for different bars within groups (in legend)
%       barcmap     - colormap for distinguishing bars within a group
%       t           - figure title
%       xl          - x-axis label
%       yl          - y-axis label
%       fontsize    - base font size
%       fontname    - name of font to use
% 
% __________________________________________________________________________
% USAGE EXAMPLE
%   data        = randn(10, 8); 
%   groupidx    = [1 2; 3 4; 5 6; 7 8]; 
%   groupn      = {'Group A' 'Group B' 'Group C' 'Group D'};
%   xl          = 'X-Axis Label'; 
%   yl          = 'Y-Axis Label';
%   t           = 'The Figure Title';
%   h  = barpatch(data, 'groupidx', groupidx, 'groupname', groupn, 'xl', xl, 'yl', yl, 't', t); 
% 

% ---------------------- Copyright (C) 2014 Bob Spunt ----------------------
% 	Created:  2015-03-09
% 	Email:    spunt@caltech.edu
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or (at
%   your option) any later version.
%       This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%       You should have received a copy of the GNU General Public License
%   along with this program.  If not, see: http://www.gnu.org/licenses/.
% __________________________________________________________________________
if nargin < 1, mfile_showhelp; return; end
nvar = size(data, 2); 
def = { 'figh',         [], ...
        'groupidx',     1:nvar, ...
        'groupname',    [], ...
        'grouptick',    0, ...
        'barname',      [], ...
        'barcmap',      'gray', ...
        'xl',           [], ...
        'yl',           [], ...
        't',            [], ...
        'fontsize',     12, ...
        'fontname',     'Arial' };
    
% | Check Varargin
% | ========================================================================
setdefaults(def, varargin);
if any(size(groupidx)==1), ngroup = 1; else ngroup = size(groupidx, 1); end

% | Compute means, ses, etc.
% | ========================================================================
nbar        = nvar/ngroup;
allm        = nanmean(data);
allse       = nansem(data);

% | Check Validity of Arguments
% | ========================================================================
if ~isempty(barname) && length(barname)~=nbar, error('Length of "barname" must match group size'); end
if ~isempty(groupname) && length(groupname)~=ngroup, error('Length of "groupname" must match number of groups'); end
if numel(groupidx)~=nvar, error('Number of indices in "groupidx" does not match number of columns in input "data"'); end

% | Check Color Map
% | ========================================================================
if ischar(barcmap)
    map     = quantile(colormap(barcmap), .1:.801/nbar:.90);
elseif isnumeric(barcmap)
    mapdim  = size(barcmap);
    if any([mapdim(1)<nbar mapdim(2)~=3])
        error('You have specified an invalid colormap in input "barcmap"');
    end
    if mapdim(1)>nbar
        fprintf('\n- WARNING: Using just the first %d rows of the input colormap - \n\n', nbar);  
    end
    map = barcmap(1:nbar,:);
end

% | Make figure
% | ========================================================================
if isempty(figh) 
    h.fig = figure('color','white', 'visible', 'off'); 
else
    h.fig = figh; 
end
h.ax = gca; 
if length(findall(h.fig, 'type', 'axes'))==1
    set(h.ax, 'position', [.075 .125 .875 .825]);
end
set(h.ax, 'FontSize', fontsize); % base font size
propfs      = round(get(h.ax, 'FontSize')*1.20); % font size proportioned wrt to axis fontsize
set(gca, 'xticklabel', []); 
h.patch     = zeros(ngroup, nbar); 
h.error     = zeros(ngroup, nbar);
h.cap       = zeros(ngroup, nbar);
xcenter     = .25+.5:1:nbar;
if ngroup > 1
    xcenter = repmat(xcenter, ngroup, 1);
    gwidth  = repmat(nbar+.5, ngroup, 1); 
    xcenter(2:end,:) = xcenter(2:end,:) + repmat(cumsum(gwidth(2:end)), 1, nbar);
    tickcenter = xcenter(:,end) + .50 + .25; 
end
xc = [xcenter-.45 xcenter+.45 xcenter+.45 xcenter-.45];  
for g = 1:ngroup 
    m = allm(groupidx(g,:)); 
    se = allse(groupidx(g,:));
    for i = 1:nbar
        xcoord  = xc(g,i:nbar:end); 
        ycoord  = [0 0 m(i) m(i)];
        h.patch(g,i) = patch(xcoord, ycoord, map(i,:), 'linestyle', 'none', 'edgecolor', map(end,:)); 
        hold on
        % - ERROR BAR
        if m(i) < 0
            ylim = [m(i)-se(i) m(i)]; 
            ycap = ylim(1); 
        else
            ylim = [m(i) m(i) + se(i)];
            ycap = ylim(2); 
        end
        h.error(g,i) = line([xcenter(g,i) xcenter(g,i)], ylim); % line
        h.cap(g,i) = line([xcenter(g,i)-.05 xcenter(g,i)+.05], [ycap ycap]); % horizontal cap
        set(h.error(g,i), 'color', get(h.ax, 'ycolor'), 'linewidth', get(h.patch(g,i), 'linewidth'));
        set(h.cap(g,i), 'color', get(h.ax, 'ycolor'), 'linewidth', get(h.patch(g,i), 'linewidth'));  
        hold on
    end
end
set(h.ax, 'units', 'normal');

% | Make legend
% | ========================================================================
if ~isempty(barname)
    h.leg = legend(h.patch(1,:), barname);
    set(h.leg, ...
        'Location', 'Best', ...
        'FontWeight', 'normal', ...
        'FontSize', ceil(propfs*1.20), ...
        'EdgeColor', get(gca, 'color'));
    hold on
end

% | LINE FOR X-AXIS (BASELINE)
% | ========================================================================
h.xline  = line(get(h.ax, 'xlim'), [0 0]);
set(h.xline, 'Color', get(h.ax, 'ycolor'), 'LineWidth', get(h.ax, 'linewidth'));
set(h.ax, 'xcolor', get(gca, 'color')); 

% | GROUP DIVING TICK MARK
% | ========================================================================
if ngroup > 1 & grouptick
   tln = .025*range(get(h.ax, 'ylim'));
   for g = 1:ngroup-1
      h.xtick(g) = line([tickcenter(g) tickcenter(g)], [-tln tln]);  
      set(h.xtick(g), 'Color', get(h.ax, 'ycolor'), 'LineWidth', get(h.ax, 'linewidth'));
   end
end

% | GROUP LABEL
% | ========================================================================
if ~isempty(groupname)
    h.grouplabel       = zeros(ngroup);
    ylim = get(h.ax, 'ylim'); 
    for g = 1:ngroup
        h.grouplabel(g) = text(0, ylim(1), groupname{g}, 'margin', 1, 'horizontalalign', 'left', 'FontSize', ceil(propfs*1.2));
        strext = get(h.grouplabel(g), 'extent');
        stradj = (nbar - strext(3))/2; 
        set(h.grouplabel(g), 'position', [xcenter(g,1)+stradj-.5 ylim(1)], 'verticalalign', 'top');  
    end
end

% | X-LABEL
% | ========================================================================
if ~isempty(xl)
    ylim = get(h.ax, 'ylim'); 
    h.xlabel = xlabel(xl, 'FontSize', ceil(propfs*1.2), 'FontWeight', 'normal');
    xpos = get(h.xlabel, 'pos'); 
    xpos(2) = ylim(1)-(.10*range(ylim)); 
    set(h.xlabel, 'pos', xpos); 
end
   
% | Y-LABEL
% | ========================================================================
if ~isempty(yl)
    h.ylabel = ylabel(yl, 'FontSize', ceil(propfs*1.2), 'FontWeight', 'normal'); 
end

% | TITLE
% | ========================================================================
if ~isempty(t)
    h.title = title(t, 'FontSize', ceil(propfs*1.5), 'FontWeight', 'normal');
end

% | FINAL CLEANUP, PROPERTY SETTING
% | ========================================================================
set(findall(h.fig, '-property', 'FontName'), 'FontName', fontname); 
set(findall(h.fig, '-property', 'units'), 'units', 'pixels');
ae = get(h.ax, 'position');
oe = get(h.ax, 'outerposition');
ae(1:2) = abs(ae(1:2)-oe(1:2))*1.5; 
if isfield(h, 'xlabel')
    xe      = get(h.xlabel, 'extent');
    ae(2)   = ceil(abs(xe(2))*1.25);
elseif isfield(h, 'grouplabel')
    xe      = get(h.grouplabel(1), 'extent');
    ae(2)   = ceil(abs(xe(2))*1.25);
end
if isfield(h, 'ylabel');
    ye      = get(h.ylabel, 'extent');
    ae(1)   = ceil(abs(ye(1))*1.25); 
end
set(h.ax, 'position', ae);
oe = get(h.ax, 'outerposition'); 
fe = get(h.fig, 'position');
fe(3) = sum(oe([1 3]))*1.01;
fe(4) = sum(oe([2 4]))*1.01; 
set(h.fig, 'pos', fe); 
set(findall(h.fig, '-property', 'units'), 'units', 'norm'); 
set(h.fig, 'visible', 'on');

end
% =========================================================================
% * SUBFUNCTIONS
% =========================================================================
function defstruct = setdefaults(def, args)
% SETDEFAULTS Name/value parsing of varargin with default values
%
%  USAGE: defstruct = setdefaults(def, args)  
% __________________________________________________________________________
%  OUTPUT
% 	defstruct: structure containing the final settingss
% __________________________________________________________________________
%  INPUTS
% 	def:       cell array containing "name, value" pairs for all varargins
% 	args:      the argument array to parse (i.e., varargin)
% __________________________________________________________________________
%  USAGE EXAMPLE (WITHIN FUNCTION)
%     defaults    = {   'arg1',     0,      ...
%                       'arg2',     'zero', ...    
%                       'arg3',     rand    ...
%                   };
%     defstruct   = setdefaults(defaults, varargin);
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-03-11
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 1, disp('USAGE: defstruct = setdefaults(def, args)'); return; end
if nargin < 2, args = []; end
def = reshape(def, 2, length(def)/2)'; 
if ~isempty(args)
    if mod(length(args), 2)
        error('Optional inputs must be entered as Name, Value pairs, e.g., myfunction(''name'', value)'); 
    end
    arg = reshape(args, 2, length(args)/2)';
    for i = 1:size(arg,1)
       idx = strncmpi(def(:,1), arg{i,1}, length(arg{i,1}));
       if sum(idx) > 1
           error(['Input "%s" matches multiple valid inputs:' repmat('  %s', 1, sum(idx))], arg{i,1}, def{idx, 1});
       elseif ~any(idx)
           error('Input "%s" does not match a valid input.', arg{i,1});
       else
           def{idx,2} = arg{i,2};
       end  
    end
end
for i = 1:size(def,1), assignin('caller', def{i,1}, def{i,2}); end
if nargout>0, defstruct = cell2struct(def(:,2), def(:,1)); end
end
function mfile_showhelp(varargin)
% MFILE_SHOWHELP
ST = dbstack('-completenames');
if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
eval(sprintf('help %s', ST(2).file));  
end