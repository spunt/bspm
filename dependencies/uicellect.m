function [theChosen, theChosenIDX] = uicellect(theCell, varargin)
% UICELLECT Present dialogue for selecting cells from a cell array
%
%     USAGE: [theChosen, theChosenIDX] = uicellect(theCell, varargin)
%
% ________________________________________________________________________________________
% OUTPUTS
%
%	  theChosen:     cell array of chosen items (empty if none chosen)
%	  theChosenIDX:  idx to input cell array of choices
%
% ________________________________________________________________________________________
% INPUTS
%
%	  theCell:  cell array of items to choose from
%
% ________________________________________________________________________________________
% VARARGIN (partial matches OK; run without arguments to see default values)
%
% | NAME            | DEFAULT       
% |-----------------|---------------------------------------------------------------------
% | DefaultIDX      | optional indices to inCell for items selected by default
% | Prompt          | message to present to user at top of gui
% | StripPath       | if true (1) and inputs are paths to files, strips path from items 
% | MultiSelect     | if true (1), user allowed to select multiple items           
% | MaxPerColumn    | max items per column (if > # of items, then one column layout)
% | RowPixelHeight  | height of gui rows (one item per row/column), in pixels     
% | ColPixelWidth   | width of gui columns, in pixels        
% | BaseFontSize    | base font size (used for item labels)  
% | hAlign          | gui horizontal alignment, can be: middle,top,upper,bottom,lower
% | vAlign          | gui vertical alignment, can be: center, left, right
% | BackgroundColor | gui background color             
% | ForegroundColor | gui foreground color
% | ItemBackgroundColor | item background color             
% | ItemForegroundColor | item foreground color
%
% ________________________________________________________________________________________
% EXAMPLES
%
%   % - Create a length 25 cell array of Items
%   theCell = cellfun(@sprintf,repmat({'Item %d'},25,1), num2cell((1:25)'),'Unif',false);
%
%   % - Present in GUI using Default Settings
%   [theChosen, theChosenIDX] = uicellect(theCell); 
%
%   % - Present in GUI and disable Multi-Selection
%   [theChosen, theChosenIDX] = uicellect(theCell, 'Multi', 0); 
%
%   % - Present in GUI but Change How it Looks
%   [theChosen, theChosenIDX] = uicellect(theCell,'MaxPer',15,'RowPix',35,'ColPix',150); 
%

% ----------------------------- Copyright (C) 2015 Bob Spunt -----------------------------
%	Created:  2015-08-23
%	Email:     spunt@caltech.edu
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
% ________________________________________________________________________________________

% | DEFAULTS FOR VARARGIN
def = { ... 
      'DefaultIDX',       [], ...
      'Prompt',           'Select from:',	...
      'StripPath',        true, ...
      'MultiSelect',      true, ...
      'MaxPerColumn',     10, ...
      'RowPixelHeight',   40, ...
      'ColPixelWidth',    250, ...
      'BaseFontSize',     14, ...
      'FontName',         'fixed-width', ...
      'hAlign',           'center', ...
      'vAlign',           'middle', ...
      'BackgroundColor',  [20/255 23/255 24/255], ...
      'ForegroundColor',  [1 1 1], ...
      'ItemBackgroundColor',  [40/255 46/255 48/255], ...
      'ItemForegroundColor',  [1 1 1] ...
	};

% | UPDATE VALUES FOR VARARGIN WHERE NECESSARY
vals = setargs(def, varargin);

% | CHECK ARGUMENTS
if nargin < 1, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
if ischar(theCell), fprintf('\nInput is not a Cell Array! Doing nothing...\n'); return; end
nopt = length(theCell);

% | DETERMINE IF INPUTS ARE PATHS TO FILES AND STRIP PATH IF NECESSARY
if size(theCell, 1)==1, theCell = theCell'; end
if all(cellfun(@exist, theCell, repmat({'file'}, nopt, 1))) && StripPath
    [theCellPath, n, e] = cellfun(@fileparts, theCell, 'Unif', false);
    theCell             = strcat(n, e);
else
    StripPath = 0; 
end

% | CALCULATE GRID
if nopt > MaxPerColumn
    ncol = ceil(nopt/MaxPerColumn);
    nrow = ceil(nopt/ncol); 
else
    ncol = 1; 
    nrow = nopt; 
end

% | CALCULATE FIGURE POSITION
uiW       = ColPixelWidth*ncol;
uiH       = RowPixelHeight*nrow;
figpos    = align_figure(uiW, uiH, vAlign, hAlign);

% | MAKE FIGURE
fige = figure( ...
            'Name'         , 'UICELLECT'     ,...
            'Units'        , 'pix'           ,...
            'Position'     , figpos          ,...
            'Resize'       , 'on'            ,...
            'Color'        , BackgroundColor ,...
            'NumberTitle'  , 'off'           ,...
            'DockControls' , 'off'           ,...
            'MenuBar'      , 'none'          ,...
            'Toolbar'      , 'none'          ,...
            'Tag'          , 'uicellect      dialogue',...
            'WindowStyle'  , 'normal'        ,...
            'UserData'     , 1);
                    
% | USE PGRID TO CREATE UIPANEL LAYOUT
pbase = pgrid(3, 1, ...
            'margin',    .01              ,...
            'panelsep',  .025             ,...
            'parent',    fige             ,...
            'relheight', [1 nrow 1]       ,...
            'backg',     BackgroundColor  ,...
            'foreg',     ForegroundColor);
set(pbase(2), 'BackgroundColor', ItemBackgroundColor, 'ForegroundColor', ItemForegroundColor);
        
% | UPPER UIPANEL:  PROMPT
ht  =  uicontrol( ...
             'Style'               , 'text'            ,...
             'Parent'              , pbase(1)          ,...
             'Tag'                 , 'Prompt'          ,...
             'String'              , Prompt            ,...
             'Units'               , 'normalized'      ,...
             'Position'            , [0 0 1 1]         ,...
             'BackgroundColor'     , BackgroundColor   ,...
             'ForegroundColor'     , ForegroundColor   ,...
             'Value'               , 0                 ,...
             'HorizontalAlignment' , 'left'            ,...
             'Enable'              , 'on'              ,...
             'FontAngle'           , 'normal'          ,...
             'FontSize'            , BaseFontSize*1.25 ,...
             'FontUnits'           , 'points'          ,...
             'FontWeight'          , 'normal');

% | MIDDLE UIPANELS:  CELL ITEMS
[popt, pidx] = pgrid(nrow, ncol, 'parent', pbase(2), 'backg', ItemBackgroundColor, 'foreg', ItemForegroundColor);

% | RE-ARRANGE SO THAT COLS ARE TRAVERSED FIRST 
pidx(:,3) = 1:length(popt);
pidx      = sortrows(pidx,2); 
popt      = popt(pidx(:,3)); 
set(pbase(2), 'bordertype', 'line');

% | CREATE UICONTROL FOR EACH CELL
if MultiSelect
    theItemCallback = '';
else
    theItemCallback = {@cb_changeselection, pbase(2)};  
end
isSelected = zeros(nopt, 1);
if DefaultIDX, isSelected(DefaultIDX) = 1; end
for i = 1:nopt
    opth(i)  =  uicontrol( ...
                   'Style'   ,     'check'                                ,...
                  'Parent'   ,     popt(i)                                ,...
                     'Tag'   ,     'opt'                                  ,...
                  'String'   ,      theCell{i}                            ,...
                   'Units'   ,     'normalized'                           ,...
                'Position'   ,     [0 0 1 1]                              ,...
         'BackgroundColor'   ,     ItemBackgroundColor                        ,...
         'ForegroundColor'   ,     ItemForegroundColor                        ,...
                   'Value'   ,     isSelected(i)                          ,...
     'HorizontalAlignment'   ,     'center'                               ,...
                  'Enable'   ,     'on'                                   ,...
           'TooltipString'   ,     ''                                     ,...
                'Callback'   ,     theItemCallback                        ,...
               'FontAngle'   ,     'normal'                               ,...
                'FontName'   ,     FontName                          ,...
                'FontSize'   ,     BaseFontSize                           ,...
               'FontUnits'   ,     'points'                               ,...
              'FontWeight'   ,     'normal'                               ,...
                'UserData'   ,     MultiSelect                            ,...
                 'Visible'   ,     'on');
     drawnow; 
end


     
% | LOWER UIPANELS:  PUSH BUTTONS
pui = pgrid(1, ncol+1, 'parent', pbase(3), 'backg', BackgroundColor, 'foreg', ForegroundColor); 

% | Select All
if MultiSelect
    htog  =  uicontrol( ...
                   'Style'   ,     'push'                               ,...
                  'Parent'   ,     pui(end-1)                             ,...
                     'Tag'   ,     'toggle'                               ,...
                  'String'   ,     'Select All'                           ,...
                   'Units'   ,     'normalized'                           ,...
                'Position'   ,     [.05 0 .90 1]                          ,...
         'BackgroundColor'   ,     ForegroundColor                        ,...
         'ForegroundColor'   ,     BackgroundColor                        ,...
                   'Value'   ,     0                                      ,...
     'HorizontalAlignment'   ,     'center'                               ,...
                  'Enable'   ,     'on'                                   ,...
           'TooltipString'   ,     ''                                     ,...
                'Callback'   ,     ''                                     ,...
               'FontAngle'   ,     'normal'                               ,...
                'FontSize'   ,     BaseFontSize*1.1                       ,...
               'FontUnits'   ,     'points'                               ,...
              'FontWeight'   ,     'bold'                                 ,...
                'UserData'   ,     []                                     ,...
                 'Visible'   ,     'on'                                    ...
                        );
end

% | FINISH
hok  =  uicontrol( ...
               'Style'   ,     'push'                                 ,...
              'Parent'   ,     pui(end)                               ,...
                 'Tag'   ,     'okbutton'                             ,...
              'String'   ,     'Done'                               ,...
               'Units'   ,     'normalized'                           ,...
            'Position'   ,     [.05 0 .90 1]                          ,...
     'BackgroundColor'   ,     ItemBackgroundColor                        ,...
     'ForegroundColor'   ,     ItemForegroundColor                        ,...
               'Value'   ,     0                                      ,...
 'HorizontalAlignment'   ,     'center'                               ,...
              'Enable'   ,     'on'                                   ,...
       'TooltipString'   ,     ''                                     ,...
            'Callback'   ,     ''                                     ,...
           'FontAngle'   ,     'normal'                               ,...
            'FontSize'   ,     BaseFontSize*1.1                       ,...
           'FontUnits'   ,     'points'                               ,...
          'FontWeight'   ,     'bold'                                 ,...
            'UserData'   ,     []                                     ,...
             'Visible'   ,     'on'                                    ...
                    );

% | FINISH UP
set(findall(fige, '-property', 'FontName'), 'FontName', FontName); 
set(fige, 'CloseRequestFcn', {@cb_closefig, fige, 0})
if MultiSelect, set(htog, 'Callback', {@cb_selectall, opth}); end
set(hok, 'Callback', {@cb_closefig, fige, 1})
drawnow
uiwait(fige)
theChosen       = [];
theChosenIDX    = [];
if StripPath, theCell = fullfile(theCellPath, theCell); end
if ishandle(fige)
    if get(fige, 'UserData')
        idx = cell2mat(get(opth, 'Value'));
        if any(idx)
            theChosen       = theCell(find(idx));
            theChosenIDX    = find(idx);
        end
    end
    delete(fige)
end
end
% ========================================================================================
%
% ------------------------------------- SUBFUNCTIONS -------------------------------------
%
% ========================================================================================
function [ph, pidx] = pgrid(nrow, ncol, varargin)
% PGRID Create a grid of of UIPANELs
%
%  USAGE: [phandle, pidx] = pgrid(nrow, ncol, varargin)
%
%  OUTPUT
%   hpanel: array of handles to uipanels comprising the grid
%   hidx:   [row,col] indices for the returned uipanel handles
% ________________________________________________________________________________________
%  INPUTS
%   nrow:   number of rows in grid
%   ncol:   number of cols in grid
% ________________________________________________________________________________________
%  VARARGIN
% | NAME            | DEFAULT       | DESCRIPTION 
% |-----------------|---------------|-----------------------------------------------------
% | parent          | gcf           | parent object for grid 
% | relwidth        | ones(1, ncol) | relative width of columns (arbitrary units)            
% | relheight       | ones(1, nrow) | relative height of rows (arbitrary units)            
% | marginsep       | 0.0100        | size of margin surrounding grid (normalized units)           
% | panelsep        | 0.0100        | size of space between panels (normalized units)            
% | backgroundcolor | [.08 .09 .09] | uipanel background color             
% | foregroundColor | [.97 .97 .97] | uipanel foreground color
% | bordertype      | 'none'        | etchedin, etchedout, beveledin, beveledout, line 
% | borderwidth     | 1             | uipanel border width in pixels
% ________________________________________________________________________________________
%

% ----------------------------- Copyright (C) 2015 Bob Spunt -----------------------------
%	Created:  2015-08-23
%	Email:     spunt@caltech.edu
% ________________________________________________________________________________________

% | Defaults for VARARGIN
def = { ...
'parent',              []                                   ,...
'relwidth',            []                                   ,...
'relheight',		   []                                   ,...
'marginsep',          .01                                   ,...
'panelsep',           .01                                   ,...
'backgroundcolor',    [20/255 23/255 24/255]                ,...
'foregroundcolor',    [248/255 248/255 248/255]             ,...
'bordertype',         'none'                                ,...
'borderwidth',         1                                     ...
};

% | Update values for VARARGIN where necessary
vals = setargs(def, varargin);

% | Check arguments
if nargin < 2, mfile_showhelp; return; end
if isempty(parent), parent          = gcf; end
if isempty(relwidth), relwidth      = ones(1, ncol); end
if isempty(relheight), relheight    = ones(1, nrow); end
if length(relwidth)~=ncol, printmsg('Length of RELWIDTH must equal NCOL. Try again!'); ph = []; pidx = []; return; end
if length(relheight)~=nrow, printmsg('Length of RELHEIGHT must equal NROW. Try again!'); ph = []; pidx = []; return; end

% | Get normalized positions for each panel 
pos         = getpositions(relwidth, relheight, marginsep, panelsep);
pidx        = pos(:,1:2);
hpos        = pos(:,3:end);

% | pgrid loop
npanel      = size(hpos, 1);
ph     = gobjects(npanel, 1);
for i = 1:npanel
    ph(i)  =  uipanel( ...
                  'Parent'   ,     parent                                 ,...
                   'Units'   ,     'normalized'                           ,...
                     'Tag'   ,     sprintf('[%d] x [%d]', pidx(i,:))        ,...
                   'Title'   ,     ''                                     ,...
           'TitlePosition'   ,     'centertop'                            ,...
                'Position'   ,     hpos(i,:)                              ,...
         'BackgroundColor'   ,     backgroundcolor                        ,...
         'ForegroundColor'   ,     foregroundcolor                        ,...
              'BorderType'   ,     bordertype                             ,...
             'BorderWidth'   ,     borderwidth                            ,...
                'UserData'   ,     pidx(i,:)                              ,...
                 'Visible'   ,     'off'                                   ...
                            );
end

% | Make Visible
for i = 1:npanel, set(ph(i), 'visible', 'on'); drawnow; end

end
function pos        = getpositions(relwidth, relheight, marginsep, uicontrolsep, top2bottomidx)
if nargin<2, relheight = [6 7]; end
if nargin<3, marginsep = .025; end
if nargin<4, uicontrolsep = .01; end
if nargin<5, top2bottomidx = 1; end
if size(relheight,1) > 1, relheight = relheight'; end
if size(relwidth, 1) > 1, relwidth = relwidth'; end
ncol        = length(relwidth);
nrow        = length(relheight); 
if top2bottomidx, relheight = relheight(end:-1:1); end

% width
rowwidth    = 1-(marginsep*2)-(uicontrolsep*(ncol-1));  
uiwidths    = (relwidth/sum(relwidth))*rowwidth;
allsep      = [marginsep repmat(uicontrolsep, 1, ncol-1)];
uilefts     = ([0 cumsum(uiwidths(1:end-1))]) + cumsum(allsep); 

% height
colheight   = 1-(marginsep*2)-(uicontrolsep*(nrow-1));
uiheights   = (relheight/sum(relheight))*colheight;
allsep      = [marginsep repmat(uicontrolsep, 1, nrow-1)];
uibottoms   = ([0 cumsum(uiheights(1:end-1))]) + cumsum(allsep);
if top2bottomidx, uiheights = uiheights(end:-1:1); end
if top2bottomidx, uibottoms = uibottoms(end:-1:1); end

% combine
pos = zeros(ncol*nrow, 6);
pos(:,1) = reshape(repmat(nrow:-1:1, ncol, 1), size(pos,1), 1);
pos(:,2) = reshape(repmat(1:ncol, 1, nrow), size(pos,1), 1);
pos(:,3) = uilefts(pos(:,2)); 
pos(:,4) = uibottoms(pos(:,1)); 
pos(:,5) = uiwidths(pos(:,2)); 
pos(:,6) = uiheights(pos(:,1));
pos      = sortrows(pos, 1);
end
function argstruct  = setargs(defaults, optargs)
% SETARGS Name/value parsing and assignment of varargin with default values
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
function figpos     = align_figure(uiW, uiH, vAlign, hAlign)
    screenPos   = get(0, 'ScreenSize');
    screenW     = screenPos(3);
    screenH     = screenPos(4);
    figpos      = [0 0 uiW uiH];
    switch lower(vAlign)
        case 'middle'
          figpos(2) = (screenH/2)-(uiH/2);
        case {'top', 'upper'}
          figpos(2) = screenH-uiH; 
        case {'bottom', 'lower'}
          figpos(2) = 1; 
        otherwise
          error('vAlign options are: middle, top, upper, bottom, lower')
    end
    switch lower(hAlign)
        case 'center'
          figpos(1) = (screenW/2) - (uiW/2);
        case 'right'
          figpos(1) = screenW - uiW; 
        case 'left'
          figpos(1) = 1; 
        otherwise
          error('hAlign options are: center, left, right')
    end
end
function cb_closefig(varargin)
  set(varargin{3}, 'UserData', varargin{4});
  uiresume(varargin{3})
end
function cb_selectall(varargin)
    h = varargin{3};
    str = {'De-Select All' 'Select All'};
    idx = strcmp(str, get(varargin{1}, 'String'));
    for i = 1:length(h)
        set(h, 'Value', find(idx)-1);
        drawnow; 
    end
    set(varargin{1}, 'String', str(~idx));
    drawnow;   
end
function cb_changeselection(varargin)
    cval    = get(varargin{1}, 'value'); 
    h       = findall(varargin{3}, 'tag', 'opt');
    set(h, 'value', 0);
    set(varargin{1}, 'value', cval); 
    drawnow;
end
function mfile_showhelp(varargin)
% MFILE_SHOWHELP
ST = dbstack('-completenames');
if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
eval(sprintf('help %s', ST(2).file));
end
