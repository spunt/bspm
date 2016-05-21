function h = cell2uitable(data, varargin)
% CELL2UITABLE Display cell array in uitable with menu option to print to CSV
% 
%  USAGE: h = cell2uitable(data, varargin)
% __________________________________________________________________________
%  INPUTS
%     data:       cell array of data (if otherwise, will convert to cell)
% __________________________________________________________________________
%  VARARGIN - run cell2uitable with no arguments to see defaults
%     parent               handle to parent (creates new fig if empty)
%     colnames             cell array of column names
%     rownames             cell array of row names
%     rowstriping          'on' | 'off'
%     fontsize             font size for table contents 
%     fontname             font name for table contents 
%     rearrangeablecols    'on' | 'off'
%     oversizecolfactor    factor to mulitply auto-computed column width
%                          (useful for making all cell contents visible)
%     addsaveuimenu        logical flag to include/exclude uimenu for saving
%     editable             used to set 'ColumnEditable' property
%     backgroundcolor      uitable background color
%     foregroundcolor      uitable foreground color
%     emptypadsize         if >0, pads with that # empty editable rows/cols 
% __________________________________________________________________________
%  EXAMPLE USAGE
%     mydata      = num2cell(randn(20, 3));
%     mycolnames  = {'Variable 1' 'Variable 2' 'Variable 3'}; 
%     h = cell2uitable(mydata, 'colnames', mycolnames); 
% 

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
% 	Created:  2015-02-03
% 	Email:    spunt@caltech.edu
% __________________________________________________________________________
def = { ...
    'parent',               [],                 ...
    'colnames',             'numbered',         ...
    'rownames',             'numbered',         ...
    'rowstriping',          'on',               ...
    'fontsize',             12,                 ...
    'fontname',             'Fixed-Width',      ...
    'rearrangeablecols',    'on',               ...
    'oversizecolfactor',    1.25,               ...
    'backgroundcolor',      [.9 .9 .9],         ...
    'foregroundcolor',      [0 0 0],            ...
    'addsaveuimenu',        true,               ...
    'editable',             true                ...
    'emptypadsize',         0                   ...
    };

% | Check Varargin
vals = setargs(def, varargin);
if nargin < 1, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
if ~iscell(data), if ischar(data), data = cellstr(data); else data = num2cell(data); end; end
data(cellfun('isempty', data)) = {''};
if and(iscell(colnames), length(colnames)~=size(data, 2)) 
    error('Length of colnames (%d) does not equal number of columns in data (%d)', length(colnames), size(data, 2));
end
if and(iscell(rownames), length(rownames)~=size(data, 1)) 
    error('Length of colnames (%d) does not equal number of columns in data (%d)', length(rownames), size(data, 1));
end

% | Create Figure
if isempty(parent)
    parent  = figure('Units', 'Norm', ...
        'pos', [0 0 1 1], ...
        'DockControls','off', ...
        'MenuBar', 'none', ...
        'Name', 'cell2uitable', ...
        'Resize', 'on', ...
        'NumberTitle', 'off', ...
        'Visible', 'off');
end

% | Menu for Save to CSV
if addsaveuimenu
    tfig = parent; 
    if ~strcmpi(get(parent, 'type'), 'figure')
        badfig = 1;
        while badfig
            tfig = get(tfig, 'parent');
            if strcmpi(get(tfig, 'type'), 'figure'), badfig = 0; end
        end
    end
    tfigmenu = uimenu(tfig,'Label','Option');
    savemenu = uimenu(tfigmenu,'Label','Write to CSV', 'CallBack', {@cb_savetable, parent});
end

% | Create Table
th = uitable('Parent', parent, ...
    'Units', 'norm', ...
    'Data', data, ...
    'ColumnName', colnames, ...
    'RowName', rownames, ...
    'Pos', [0 0 1 1], ...
    'RearrangeableColumns', rearrangeablecols, ...
    'ColumnEditable', editable, ...
    'BackgroundColor', backgroundcolor, ...
    'ForegroundColor', foregroundcolor, ...
    'RowStriping', rowstriping, ...
    'FontUnits', 'Points', ...
    'FontSize', fontsize, ...
    'FontName', fontname);

% | COLUMN WIDTHS
strdata     = get(th, 'data');
[nrow,ncol] = size(strdata); 
colwidth    = cell(1,ncol);
for i = 1:ncol
   ln = cellfun('length', strdata(:,i)); 
   idx = find(ln==max(ln)); 
   colwidth{i} = strsize(strdata(idx(1), i), 'FontSize', fontsize', 'FontName', fontname)*oversizecolfactor; 
end
% | PADDING 
if emptypadsize
   pdata = cell(size(strdata)+emptypadsize);
   pdata(:) = {''};
   pdata(1:nrow, 1:ncol) = strdata;
   colwidth = [colwidth repmat({'auto'}, 1, emptypadsize)]; 
   set(th, 'data', pdata);  
end
set(th, 'columnwidth', colwidth); 
% | Resize Table to Fit Figure
set(th, 'units', 'pix');
set(parent, 'units', 'pix'); 
tpos    = get(th, 'extent');
fpos    = get(parent, 'pos'); 
figpos  = align_figure(tpos(3), tpos(4), 'middle', 'center'); 
set(parent, 'pos', figpos);
set(th, 'units', 'norm');
set(findall(parent, '-property', 'units'), 'units', 'norm');
set(th, 'pos', [0 0 1 1]); 
set(parent, 'Visible', 'on');

% | Save Handles
if nargout
    h.fig   = parent;
    h.tab   = th;
    if addsaveuimenu, h.menu  = savemenu; end
end
drawnow; 
end
% ==========================================================================
%
% ------------------------------ SUBFUNCTIONS ------------------------------
%
% ==========================================================================
function cb_savetable(varargin)
    t = findobj(varargin{3}, 'type', 'uitable');
    data = get(t, 'data'); 
    colnames = get(t, 'columnname'); 
    if size(colnames,2)~=size(data,2), colnames = colnames'; end
    outcell = [colnames; data]; 
    outname = sprintf('%s.csv', regexprep(get(varargin{3}, 'name'), ' ', '_')); 
    [fname, pname] = uiputfile({'*.csv', 'CSV File'; '*.*', 'All Files (*.*)'}, 'Save Table As', outname);
    writereport(outcell, fullfile(pname, fname)); 
end
function writereport(incell, outname)
% WRITEREPORT Write cell array to CSV file
%
%  USAGE: outname = writereport(incell, outname)	*optional input
% __________________________________________________________________________
%  INPUTS
%	incell:     cell array of character arrays
%	outname:   base name for output csv file 
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-02-02
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 2, disp('USAGE: outname = writereport(incell, outname)'); return; end

% | Convert all cell contents to character arrays
% | ========================================================================
[nrow, ncol] = size(incell);
for i = 1:numel(incell)
    if isnumeric(incell{i}), incell{i} = num2str(incell{i}); end
    if strcmp(incell{i},'NaN'), incell{i} = ''; end
end
incell = regexprep(incell, ',', '');

% | Write to file
% | ========================================================================
fid = fopen(outname,'w');
for r = 1:nrow
    fprintf(fid,['%s' repmat(',%s',1,ncol-1) '\n'],incell{r,:});
end
fclose(fid);
end
function figpos = align_figure(uiW, uiH, valign, halign)
    screenPos   = get(0, 'ScreenSize');
    screenW     = screenPos(3);
    screenH     = screenPos(4);
    figpos      = [0 0 uiW uiH];
    switch lower(valign)
        case 'middle'
          figpos(2) = (screenH/2)-(uiH/2);
        case {'top', 'upper'}
          figpos(2) = screenH-uiH; 
        case {'bottom', 'lower'}
          figpos(2) = 1; 
        otherwise
          error('VALIGN options are: middle, top, upper, bottom, lower')
    end
    switch lower(halign)
        case 'center'
          figpos(1) = (screenW/2) - (uiW/2);
        case 'right'
          figpos(1) = screenW - uiW; 
        case 'left'
          figpos(1) = 1; 
        otherwise
          error('HALIGN options are: center, left, right')
    end
end
function argstruct = setargs(defaults, optargs)
% SETARGS Name/value parsing and assignment of varargin with default values
% 
% This is a utility for setting the value of optional arguments to a
% function. The first argument is required and should be a cell array of
% "name, default value" pairs for all optional arguments. The second
% argument is optional and should be a cell array of "name, custom value"
% pairs for at least one of the optional arguments.
% 
%  USAGE: argstruct = setargs(defaults, args)  
% __________________________________________________________________________
%  OUTPUT
% 
% 	argstruct: structure containing the final argument values
% __________________________________________________________________________
%  INPUTS
% 
% 	defaults:  
%       cell array of "name, default value" pairs for all optional arguments
% 
% 	optargs [optional]     
%       cell array of "name, custom value" pairs for at least one of the
%       optional arguments. this will typically be the "varargin" array. 
% __________________________________________________________________________
%  USAGE EXAMPLE (WITHIN FUNCTION)
% 
%     defaults    = {'arg1', 0, 'arg2', 'words', 'arg3', rand}; 
%     argstruct   = setargs(defaults, varargin)
%


% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
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
% __________________________________________________________________________
if nargin < 1, mfile_showhelp; return; end
if nargin < 2, optargs = []; end
defaults = reshape(defaults, 2, length(defaults)/2)'; 
if ~isempty(optargs)
    if mod(length(optargs), 2)
        error('Optional inputs must be entered as Name, Value pairs, e.g., myfunction(''name'', value)'); 
    end
    arg = reshape(optargs, 2, length(optargs)/2)';
    for i = 1:size(arg,1)
       idx = strncmpi(defaults(:,1), arg{i,1}, length(arg{i,1}));
       if sum(idx) > 1
           error(['Input "%s" matches multiple valid inputs:' repmat('  %s', 1, sum(idx))], arg{i,1}, defaults{idx, 1});
       elseif ~any(idx)
           error('Input "%s" does not match a valid input.', arg{i,1});
       else
           defaults{idx,2} = arg{i,2};
       end  
    end
end
for i = 1:size(defaults,1), assignin('caller', defaults{i,1}, defaults{i,2}); end
if nargout>0, argstruct = cell2struct(defaults(:,2), defaults(:,1)); end
end
function mfile_showhelp(varargin)
% MFILE_SHOWHELP
ST = dbstack('-completenames');
if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
eval(sprintf('help %s', ST(2).file));  
end
function out = cellnum2str(in, ndec)
% CELLNUM2STR 
%
%  USAGE: out = cellnum2str(in, ndec)
% __________________________________________________________________________
%  INPUTS
%	in:     numeric cell array
%   ndec:   number of decimal points to display
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-01-13
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 2, ndec = 3; end
if nargin < 1, mfile_showhelp; return; end
if ~iscell(in), error('Input array must be cell!'); end
out = cellfun(@num2str, in, repmat({['%2.' num2str(ndec) 'f']}, size(in)), 'Unif', false); 
out = regexprep(out, '0\.', '\.');
end
function [strw, strh] = strsize(string, varargin)
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
strh    = textExt(4);
strw    = textExt(3);

% | If using a proportional font, shrink text width by a fudge factor to account for kerning.
if ~strcmpi(FontName, 'Fixed-Width'), strw = strw*0.9; end 
end
