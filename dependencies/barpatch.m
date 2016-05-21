function h = barpatch(data, varargin)
% BARPATCH Create bar graph with error bars using patch and line objects
% 
% This function will create a grouped bar graph with error bars without
% using the standard plotting functions BAR and ERRORBAR. It uses PATCH to
% create the bars and LINE to construct the error bars.
% 
%  USAGE: h = barpatch(data, varargin)
% __________________________________________________________________________
%  OUTPUT
% 	h: handles to all graphics objects  
% __________________________________________________________________________
%  INPUT
%   data: data matrix to plot; rows are cases, cols are variables 
% __________________________________________________________________________
%  OPTIONAL INPUTS (VARARGIN)
%   These are entered as `'name', value` argument pairs. Matching is not
%   case-sensitive and partial name matches are OK. To see default values,
%   run barpatch without any arguments.
%
%   | Name            | Description                                         |
%   | -----------     | --------------------------------------------------- |
%   | figh            | handle for figure to plot in                        |
%   | newfig          | flag to create and plot in a new figure             |
%   | groupidx        | rows index columns of "data" to plot as a group     |
%   | groupname       | labels for different groups of bars                 |
%   | grouptick       | flag to place tickmark between groups on x-axis     |
%   | groupspace      | controls spacing between groups of bars             |
%   | barname         | labels for different bars within groups (in legend) |
%   | barcmap         | colormap for distinguishing bars within a group     |
%   | barwidth,       | width of bars (>1 produces overlapping bars)        |
%   | errlinewidth    | width of error bar lines                            |
%   | detachlegend    | flag to plot legend in new figure (if applicable)   |
%   | t               | figure title                                        |
%   | xl              | x-axis label                                        |
%   | yl              | y-axis label                                        |
%   | fontsize        | base font size                                      |
%   | fontname        | name of font to use                                 |
%   | fontunits       | font units of finished product                      |
%   | ytickformat     | display formatting for yticklabels (e.g., '%.2f')   |
%   | yticklength     | # of yticks (if empty, determined automatically)    |
%   | fontmultiplier1 | multiplier for next largest size from basefontsize  |
%   | fontmultiplier2 | multiplier for next largest size from basefontsize  |
%   | fontmultiplier3 | multiplier for next largest size from basefontsize  |  
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
def = { ...
    'newfig',           1,                  ...
    'groupidx',         [],                 ...
    'groupname',        [],                 ...
    'grouptick',        0,                  ...
    'groupspace'        0.5,                ...
    'barname',          [],                 ...
    'barcmap',          'gray',             ...
    'barwidth',         .95,                ...
    'detachlegend',     0,                  ...
    'errlinewidth',     1.5,                ...
    'xl',               [],                 ...
    'yl',               [],                 ...
    't',                [],                 ...
    'fontunits',        'norm',             ...
    'fontsize',         12,                 ...
    'fontmultiplier1',  1.2,                ...
    'fontmultiplier2',  1.4,                ...
    'fontmultiplier3',  1.6,                ...
    'fontname',         'Arial',            ...
    'ytickformat',      '%.2f',             ...
    'yticklength',      []                  ...
     };
    
% | Check Varargin
% | ========================================================================
vals = setargs(def, varargin);
if nargin<1, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
nvar = size(data, 2); 
if isempty(groupidx), groupidx = 1:nvar; end
if any(size(groupidx)==1), ngroup = 1; else ngroup = size(groupidx, 1); end

% | Compute means, ses, etc.
% | ========================================================================
nbar        = nvar/ngroup;
allm        = nanmean(data);
allse       = nansem(data);

% | Check Validity of Arguments
% | ========================================================================
if ~isempty(barname) && length(barname)~=nbar, error('Length of "barname" must match group size'); end
if and(~isempty(groupname), ischar(groupname)), groupname = cellstr(groupname); end
if and(~isempty(groupname), length(groupname)~=ngroup), error('Length of "groupname" must match number of groups'); end
if numel(groupidx)~=nvar, error('Number of indices in "groupidx" does not match number of columns in input "data"'); end

% | Grab or Make Figure
% | ========================================================================
if newfig
    h.fig   = figure('color','white');
    h.ax    = gca;
    set(h.ax, 'position', [.075 .125 .875 .825]);
else
    h.fig   = gcf;
    h.ax    = gca; 
end
set(findall(h.fig, '-property', 'units'), 'units', 'norm');
set(h.ax, 'FontSize', fontsize); % base font size
propfs      = round(get(h.ax, 'FontSize')*fontmultiplier1); % font size proportioned wrt to axis fontsize
set(gca, 'xticklabel', []);

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

% | Create the Pretty Patch of Bars with Lines, Labels, Tigers, Bears, etc.
% | ========================================================================
h.patch     = zeros(ngroup, nbar); 
h.error     = zeros(ngroup, nbar);
h.cap       = zeros(ngroup, nbar);
xcenter     = .25+.5:1:nbar;
if ngroup > 1
    xcenter = repmat(xcenter, ngroup, 1);
    gwidth  = repmat(nbar + groupspace, ngroup, 1); 
    xcenter(2:end,:) = xcenter(2:end,:) + repmat(cumsum(gwidth(2:end)), 1, nbar);
    tickcenter = xcenter(:,end) + .50 + .25;
    gcenter     = sum(xcenter(:,[1 end]), 2)/2;
    fcenter     = sum(gcenter)/ngroup; 
else
    gcenter = sum(xcenter)/nbar;
    fcenter = sum(xcenter)/nbar;
end
halfbar     = barwidth/2; 
xc          = [xcenter-halfbar xcenter+halfbar xcenter+halfbar xcenter-halfbar];
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
        set(h.error(g,i), 'color', get(h.ax, 'ycolor'), 'linewidth', errlinewidth);
        set(h.cap(g,i), 'color', get(h.ax, 'ycolor'), 'linewidth', errlinewidth);  
        
%         set(h.error(g,i), 'color', get(h.ax, 'ycolor'), 'linewidth', get(h.patch(g,i), 'linewidth'));
%         set(h.cap(g,i), 'color', get(h.ax, 'ycolor'), 'linewidth', get(h.patch(g,i), 'linewidth'));  
        hold on
    end
end
set(h.ax, 'units', 'normal');

% | Make legend
% | ========================================================================
if ~isempty(barname)
    
    if detachlegend
    
        h.leg   = figure('color','white', 'units', 'pix');
        h.legax = gca;
        w = .9*(1/nbar);
        ybottom = ((1/nbar)-w)/2:(1/nbar):1;
        ytop = ybottom + w;
        xleft = ((1/nbar)-w)/2;
        xright = xleft + w; 
        ycoord = [ybottom' ybottom' ytop' ytop'];
        xcoord = [xleft' xright' xright' xleft'];
        for i = 1:length(barname)
            p(i) = patch(xcoord, ycoord(i,:), map(i,:), 'linestyle', 'none', 'edgecolor', map(end,:)); 
        end
        set(h.legax, 'ylim', [0 1], 'xlim', [0 1]);
        axis off; 
        labx = (1/nbar) + xleft;
        laby = ((1/nbar)/2):(1/nbar):1; 
        for g = 1:nbar
            h.leglabel(g) = text(labx, laby(g), barname{g}, 'margin', 1, 'horizontalalign', 'left',  'verticalalign', 'middle', 'FontSize', ceil(propfs*fontmultiplier1)); 
        end
        axes(h.ax);
        
    else
        
        h.leg = legend(h.patch(1,:), barname);
        set(h.leg, ...
            'Location', 'Best', ...
            'FontWeight', 'normal', ...
            'FontSize', ceil(propfs*fontmultiplier2), ...
            'EdgeColor', get(gca, 'color'));
        hold on
        
    end
end

% | LINE FOR X-AXIS (BASELINE)
% | ========================================================================
h.xline  = line(get(h.ax, 'xlim'), [0 0]);
set(h.xline, 'Color', get(h.ax, 'ycolor'), 'LineWidth', get(h.ax, 'linewidth'));
set(h.ax, 'xcolor', get(gca, 'color')); 

% | GROUP DIVIDING TICK MARK
% | ========================================================================
if and(ngroup>1, grouptick)
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
    barlim = cell2mat(get(h.patch, 'ydata'));
    errlim = cell2mat(get(h.error, 'ydata'));
    alllim = [barlim(:); errlim(:)];
    alllim = [min(alllim(:)) max(alllim(:))];
    for g = 1:ngroup
        h.grouplabel(g) = text(0, alllim(1), groupname{g}, 'margin', 1, 'horizontalalign', 'center', 'FontSize', ceil(propfs*fontmultiplier1));
        set(h.grouplabel(g), 'position', [gcenter(g,1) alllim(1)], 'verticalalign', 'top', 'tag', 'groupname', 'horizontalalign', 'center');
    end
    ext = get(h.grouplabel(1), 'extent');
    ylim = get(h.ax, 'ylim');
    ylim(1) = ext(2);
    set(h.ax, 'ylim', ylim);
end

% | X-LABEL
% | ========================================================================
if ~isempty(xl)
    ylim        = get(h.ax, 'ylim');
    if ~isempty(groupname), ylim(1) = ylim(1) - (range(ylim)*.05); end
    h.xlabel    = text(fcenter, ylim(1), xl, ...
        'margin', 1, ...
        'FontSize', ceil(propfs*fontmultiplier2), ...
        'horizontalalign', 'center', ...
        'verticalalign', 'top', ...
        'tag', 'xlabel' ...
    );
    ext = get(h.xlabel, 'extent');
    ylim = get(h.ax, 'ylim');
    ylim(1) = ext(2);
    set(h.ax, 'ylim', ylim);
end
   
% | Y-LABEL
% | ========================================================================
if ~isempty(yl)
    h.ylabel = ylabel(yl, 'FontSize', ceil(propfs*fontmultiplier1), 'FontWeight', 'normal'); 
end

% | TITLE
% | ========================================================================
if ~isempty(t) 
    barlim = cell2mat(get(h.patch, 'ydata'));
    errlim = cell2mat(get(h.error, 'ydata'));
    alllim = [barlim(:); errlim(:)];
    alllim = [min(alllim(:)) max(alllim(:))];
    h.title    = text(fcenter, alllim(2)+(.025*range(alllim)), t, ...
        'margin', 1, ...
        'FontSize', ceil(propfs*fontmultiplier3), ...
        'horizontalalign', 'center', ...
        'verticalalign', 'bottom', ...
        'tag', 'title' ...
    );
    ext     = get(h.title(1), 'extent');
    ylim    = get(h.ax, 'ylim');
    ylim(2) = sum(ext([2 4])); 
    set(h.ax, 'ylim', ylim);
end

% | FORMAT TICKLABELS
% | ========================================================================
if yticklength
    ylim    = get(h.ax, 'ylim');
    set(h.ax, 'ytick', linspace(ylim(1), ylim(2), yticklength));
end
if and(~isempty(ytickformat), ~isempty(get(h.ax, 'yticklabel')))
    yt  = get(h.ax, 'ytick');
    yts = cell(size(yt)); 
    for i = 1:length(yt)
       yts{i} = sprintf(ytickformat, yt(i));
    end
    set(h.ax, 'yticklabel', yts);
end


% | FINAL CLEANUP, PROPERTY SETTING
% | ========================================================================
set(findall(h.fig, '-property', 'FontName'), 'FontName', fontname, 'FontUnits', 'norm'); 
set(findall(h.fig, '-property', 'units'), 'units', fontunits);
if newfig, set(h.ax, 'OuterPosition', [0 0 1 1]); end

end
% ==========================================================================
%
% ------------------------------ SUBFUNCTIONS ------------------------------
%
% ==========================================================================
function argstruct = setargs(defaultargs, varargs)
% SETARGS Name/value parsing and assignment of varargin with default values
% 
% This is a utility for setting the value of optional arguments to a
% function. The first argument is required and should be a cell array of
% "name, default value" pairs for all optional arguments. The second
% argument is optional and should be a cell array of "name, custom value"
% pairs for at least one of the optional arguments.
% 
%   USAGE: argstruct = setargs(defaultargs, varargs)
% __________________________________________________________________________
% OUTPUT
% 
% 	ARGSTRUCT
%    structure containing the final argument values
% __________________________________________________________________________
% INPUTS
% 
% 	DEFAULTARGS  
%     cell array of "'Name', value" pairs for all variables with default
%     values
% 
% 	VARARGS [optional]     
%     cell array of user-specified "'Name', value" pairs for one or more of
%     the variables with default values. this will typically be the
%     "varargin" cell array. for each pair, SETARGS determines if the
%     specified variable name can be uniquely matched to one of the default
%     variable names specified in DEFAULTARGS. matching uses STRNCMPI and
%     thus is case-insensitive and open to partial name matches (e.g.,
%     default variable name 'FontWeight' would be matched by 'fontweight',
%     'Fontw', etc.). if a match is found, the user-specified value is then
%     used in place of the default value. if no match is found or if
%     multiple matches are found, SETARGS returns an error and displays in
%     the command window information about the argument that caused the
%     problem.
% __________________________________________________________________________
% USAGE EXAMPLE (TO BE USED AT TOP OF FUNCTION WITH VARARGIN)
% 
%     defaultargs = {'arg1', 0, 'arg2', 'words', 'arg3', rand}; 
%     argstruct   = setargs(defaultargs, varargin)
%


% ---------------------- Copyright (C) 2015 Bob Spunt -----------------------
%	Created:  2015-03-11
%	Email:    spunt@caltech.edu
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
% 
if nargin < 1, mfile_showhelp; return; end
if nargin < 2, varargs = []; end
defaultargs = reshape(defaultargs, 2, length(defaultargs)/2)'; 
if ~isempty(varargs)
    if mod(length(varargs), 2)
        error('Optional inputs must be entered as "''Name'', Value" pairs, e.g., myfunction(''arg1'', val1, ''arg2'', val2)'); 
    end
    arg = reshape(varargs, 2, length(varargs)/2)';
    for i = 1:size(arg,1)
       idx = strncmpi(defaultargs(:,1), arg{i,1}, length(arg{i,1}));
       if sum(idx) > 1
           error(['Input "%s" matches multiple valid inputs:' repmat('  %s', 1, sum(idx))], arg{i,1}, defaultargs{idx, 1});
       elseif ~any(idx)
           error('Input "%s" does not match a valid input.', arg{i,1});
       else
           defaultargs{idx,2} = arg{i,2};
       end
    end
end
for i = 1:size(defaultargs,1), assignin('caller', defaultargs{i,1}, defaultargs{i,2}); end
if nargout>0, argstruct = cell2struct(defaultargs(:,2), defaultargs(:,1)); end
end
function mfile_showhelp(varargin)
% MFILE_SHOWHELP
ST = dbstack('-completenames');
if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
eval(sprintf('help %s', ST(2).file));  
end