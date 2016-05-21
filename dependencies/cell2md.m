function varargout = cell2md(X, varargin)
% CELL2MD Save a cell array as a markdown table
% 
%  USAGE: cell2md(X, varargin)
% __________________________________________________________________________
%  INPUTS
%   X:              cell array to write (can be numeric or char array, in 
%                   which case it will be converted to a cell array)
%   VARARGIN:       optional arguments entered as "name,value" pairs: 
%       outfile     - name for output file (default: cell2md_YYYY_MM_dd)
%       hdrnames    - name for column headers (default: first row of X)
%       alignment   - horz alignment for columns (can be 1x1 char to use 
%                     same alignment for all columns, or 1 x size(X,2) cell 
%                     array to use column-specific alignments)
% __________________________________________________________________________
%  EXAMPLE USAGE
%    X = num2cell(rand(10, 3));
%    Save with default filename and alignment with column headers A, B, C
%    cell2md(X, 'hdr', {'A' 'B' 'C'});
% 
% 
% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
% 	Created:  2015-06-11
% 	Email:    spunt@caltech.edu
% __________________________________________________________________________

% | Default VARARGIN
defaults = {
            'outfile',         fullfile(pwd, ['cell2md_' datestr(now, 'YYYY-mm-DD') '.md']) ,            ...
            'hdrnames',        [],            ...
            'alignment',       'center',              ...
            'dispoutput',         true            ...
            };
vals = setargs(defaults, varargin);
if nargin==0, mfile_showhelp; fprintf('\t= DEFAULT VARARGIN =\n'); disp(vals); return; end

% | Check X
if ~iscell(X)
    fprintf('| WARNING: Input is not a cell array and will be converted to cell. Check output carefully.\n'); 
    if ischar(X) 
        X = cellstr(X); 
    elseif isnumeric(X) 
        X = cellfun(@num2str, num2cell(X), 'Unif', false); 
    end
else
    cellclass = cellfun(@class, X, 'Unif', false);
    notchar = ~strcmpi(cellclass, 'char'); 
    if any(notchar(:))
        fprintf('| WARNING: Some cells do not contain CHAR data and will be converted to CHAR. Check output carefully.\n');
        try
            X(notchar) = cellfun(@num2str, X(notchar), 'Unif', false); 
        catch
           fprintf('| ERROR: Could not converted some cells to CHAR. Exiting.\n'); return;  
        end
    end
end

% | Check HDRNAMES
if isempty(hdrnames)
  hdrnames = X(1,:);
  X(1,:) = [];
elseif length(hdrnames)~=size(X, 2)
  error('Length of colnames (%d) does not equal number of columns in X (%d)', length(hdrnames), size(X, 2)); 
end

% | Check ALIGN
if ischar(alignment), alignment = cellstr(alignment); end
if length(alignment)==1
  alignment = repmat(alignment, 1, size(X, 2));
elseif length(alignment)~=size(X, 2)
  error('Length of ALIGNMENT (%d) does not equal number of columns in X (%d)', length(alignment), size(X, 2)); 
end

% | Check OUTFILE
[p,n,e] = fileparts(outfile); 
if isempty(p), p = pwd; end
if isempty(e), e = '.md'; end
outfile = fullfile(p,[n e]); 

% | X to MarkDown Table
cellln = cellfun('length', [hdrnames; X]);
colsizes = max(cellln); 
wspace = repmat(colsizes, size(cellln, 1), 1) - cellln; 
mdcell = cell(size(X,1) + 2, size(X, 2)*2 + 1);
mdcell(:,1:2:end) = repmat({'|'}, size(X, 1) + 2, size(X, 2)+1);
colpos = 2:2:size(mdcell, 2);
alignlab = {'left' 'center' 'right'}; 
alignopt = {':-' '::' '-:'};  
for i = 1:length(colpos)
  div = repmat('-', 1, colsizes(i)); 
  div([1 end]) = alignopt{strcmpi(alignlab, alignment{i})}; 
  mdcell{2, colpos(i)} = div; 
  mdcell{1, colpos(i)} = [hdrnames{i} repmat(' ', 1, wspace(1,i))];
  for ii = 1:size(X, 1) 
    mdcell{ii + 2, colpos(i)} = [X{ii, i} repmat(' ', 1, wspace(ii+1,i))];
  end
end
if nargout>0, varargout{1} = mdcell; end

% | Write to file
if ~dispoutput
    fprintf('%s', horzcat(mdcell{1, :})); 
    for r = 2:size(mdcell, 1), fprintf('\n%s', horzcat(mdcell{r, :})); end
end

fid = fopen(outfile, 'w');
fprintf(fid, '%s', horzcat(mdcell{1, :})); 
for r = 2:size(mdcell, 1), fprintf(fid, '\n%s', horzcat(mdcell{r, :})); end
fclose(fid);
fprintf('| OUTPUT: %s\n', outfile); 

end
% ==========================================================================
%
% ------------------------------ SUBFUNCTIONS ------------------------------
%
% ==========================================================================
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

