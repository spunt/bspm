 function newfunc(name, varargin)
% NEWFUNC Customize and create a formatted file for a new MATLAB function
%
%  USAGE:   newfunc(name, varargin)
%           *newfunc with no args prints help and shows default VARARGINs
% 
% _________________________________________________________________________
%  NECESSARY ARGUMENT
%     name          = function name e.g. ['myfunc.m','myfunc',{'myfunc'}]
%                     if path is omitted, arg "outdir" is used (see below)
% 
% _________________________________________________________________________
%  VARARGIN (run newfunc w/no arguments to see default values)
%     descrip       = brief description to include at top of doc
%     numargin      = number of necessary arguments
%     numvarargin   = number of optional arguments (name-value pairs)
%     numargout     = number of outputs
%     usage_example = flag to incl USAGE EXAMPLE section in doc
%     credits       = flag to incl CREDITS section in doc
%     authorname    = author name for Copyright seciton of doc
%     authoremail   = author email for Copyright seciton of doc
%     outdir        = path to save the new function file
%     linewidth     = width (in chars) for doc and section dividers
%     editafter     = flag to open created file in default editor
% 
% _________________________________________________________________________
%  EXAMPLES
%     newfunc
%     newfunc('mynewfunc'); 
%     newfunc('mynewfunc', ...
%             'descrip',          'This is my new function', ...
%             'numargin',         2, ...
%             'numvarargin',      6, ...
%             'numargout',        1, ...
%             'usage_example',    true, ...
%             'credits',          true  ...
%             );         
% 

% ---------------------- Copyright (C) 2014 Bob Spunt ----------------------
%	Created:  2014-09-27
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
% _________________________________________________________________________

% | DEFAULT VALUES FOR VARARGIN NAME-VALUE PAIRS
% | Customize these  to save time (e.g., author name, email, outdir)
def = { ...
        'descrip',          '',         ...
        'numargin',         1,          ...
        'numvarargin',      3,          ...
        'numargout',        0,          ...
        'usage_example',    false,      ...
        'credits',          false,      ...
        'authorname',       'Bob Spunt',   ...
        'authoremail',      'spunt@caltech.edu', ...
        'outdir',           '~/Github/matlab-general/library',        ...
        'linewidth',        90,         ...
        'editafter',        true,       ...
       };
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
option = [numargout usage_example credits];
subfunctions = true; 
if iscell(name), name = char(name);end
if iscell(descrip), descrip = char(descrip); end

% | CREATE OPEN-MFILE
% | =======================================================================
[tmp, name]     = fileparts(name);
outname         = fullfile(outdir, [name '.m']);
if exist(outname, 'file')
    resp = inputoption('File by that name exists. Overwrite it?', {'y' 'n'}); 
    if strcmp(resp, 'n'), return; end
end
fid                 = fopen(outname,'w+');

arginstr            = strtrim(sprintf(repmat('arg%d, ', 1, numargin), 1:numargin));
if numvarargin
    arginstr        = [arginstr ' varargin']; 
else
    arginstr(end)   = []; 
end
if numargout
    argoutstr       = strtrim(sprintf(repmat('argout%d, ', 1, numargout), 1:numargout));
    argoutstr(end)  = []; 
    if numargout > 1, argoutstr = ['[' argoutstr ']']; end
    funcstring  = sprintf('%s = %s(%s)', argoutstr, name, arginstr); 
else
    funcstring = sprintf('%s(%s)', name, arginstr); 
end
fprintf(fid, 'function %s', funcstring);

% | CREATE FORMATTED SECTION DIVIDERS & TITLES
% | =======================================================================
crln = ['Copyright (C) ' datestr(now,'YYYY') ' ' authorname];
sfln = 'SUBFUNCTIONS'; 
sidelength = floor(((linewidth) - length(crln) - 1)/2) - 1;
crln = ['% ' repmat('-',1,sidelength) ' ' crln ' ' repmat('-',1,sidelength)];
sidelength = floor(((linewidth) - length(sfln) - 1)/2) - 1;
sfln = ['% ' repmat('-',1,sidelength) ' ' sfln ' ' repmat('-',1,sidelength)];
fmtline1 = ['% ' sprintf(repmat('_',1,linewidth-2))];
fmtline2 = sprintf('%% | %s', repmat('=', 1, linewidth-4));
fmtline3 = sprintf('%% %s', repmat('=', 1, linewidth-2)); 
all = {crln sfln fmtline1 fmtline2 fmtline3};
allln = cellfun('length', all); 
if length(unique(allln)) > 1
    linewidth = max(allln);
    crln = ['Copyright (C) ' datestr(now,'YYYY') ' ' authorname];
    sfln = 'SUBFUNCTIONS'; 
    sidelength = floor(((linewidth) - length(crln) - 1)/2) - 1;
    crln = ['% ' repmat('-',1,sidelength) ' ' crln ' ' repmat('-',1,sidelength)];
    sidelength = floor(((linewidth) - length(sfln) - 1)/2) - 1;
    sfln = ['% ' repmat('-',1,sidelength) ' ' sfln ' ' repmat('-',1,sidelength)];
    fmtline1 = ['% ' sprintf(repmat('_',1,linewidth-2))];
    fmtline2 = sprintf('%% | %s', repmat('=', 1, linewidth-4));
    fmtline3 = sprintf('%% %s', repmat('=', 1, linewidth-2));    
end

% | DOCUMENTATION
% | =======================================================================
fprintf(fid, '\n%% %s %s\n%%\n%%  USAGE: %s', upper(name), descrip, funcstring);
if option(1), fprintf(fid, ['\n%%\n%%  OUTPUT\n' repmat('%%\targout%d:  \n', 1, option(1)) '%%'], 1:numargout); end
fprintf(fid, ['\n%s\n%%  INPUTS\n' repmat('%%\targ%d:  \n', 1, numargin) '%%'], fmtline1, 1:numargin);
if numvarargin, fprintf(fid, ['\n%s\n%%  VARARGIN\n' repmat('%%\tvararg%d:  \n', 1, numvarargin) '%%'], fmtline1, 1:numvarargin); end
if option(2), fprintf(fid, '\n%s\n%%  EXAMPLES\n%%\t>> %s\n%%', fmtline1, name); end
if option(3), fprintf(fid, '\n%s\n%%  CREDITS\n%%\tProps go to...\n%%', fmtline1); end

% | COPYRIGHT
% | =======================================================================
fprintf(fid, '\n\n%s', crln);
fprintf(fid, '\n%%\tCreated:  %s', datestr(now, 'YYYY-mm-DD')); 
fprintf(fid, '\n%%\tEmail:     %s', authoremail);
fprintf(fid, '\n%s', fmtline1);

% | BODY
% | =======================================================================
if numvarargin
    fprintf(fid, '\ndef = { ... '); 
    for i = 1:numvarargin-1
        fprintf(fid, '\n\t''vararg%d'',\t\t''def%d'',\t...', i, i); 
    end
    fprintf(fid, '\n\t''vararg%d'',\t\t''def%d''\t...\n\t};', numvarargin, numvarargin);
    fprintf(fid, '\nvals = setargs(def, varargin);'); 
    fprintf(fid, '\nif nargin < %d, mfile_showhelp; fprintf(''\\t| - VARARGIN DEFAULTS - |\\n''); disp(vals); return; end\n', numargin);
else
    fprintf(fid, '\nif nargin < %d, mfile_showhelp; return; end\n', numargin);
end
fprintf(fid, '\n%% | SECTION\n%s', fmtline2); 
fprintf(fid, '\n\n%% | Subsection');
fprintf(fid, '\n\nend'); 

% | SUBFUNCTIONS
% | =======================================================================
if subfunctions
    fprintf(fid, '\n\n\n%s\n%%\n%s\n%%\n%s\n', fmtline3, sfln, fmtline3);
    m = mfile_showhelp_contents; 
    for i = 1:length(m)
       fprintf(fid, '%s\n', m{i});  
    end   
    if numvarargin
        m = setargs_contents; 
        for i = 1:length(m)
           fprintf(fid, '%s\n', m{i});  
        end
    end
end

% | Done
fclose(fid);
fprintf('\nFunction saved to: %s\n', name); 
if editafter, edit(outname); end
end
% =========================================================================
% *
% * SUBFUNCTIONS
% *
% =========================================================================
function resp       = inputoption(prompt, respset)
% INPUTOPTION Like INPUT, but with default and response option set
%
%  USAGE: resp = inputoption(prompt, respset)
% __________________________________________________________________________
%  INPUTS
%	prompt:     message to user 
%	respset:    valid responses (case is ignored when validating)
%
% __________________________________________________________________________
%  EXAMPLES
%	>> resp = inputoption('Would you like to quit?', {'Y' 'N'}, 'Y')
%	>> resp = inputoption('1, 2, or 3 apples?', [1 2 3])
%

% ---------------------- Copyright (C) 2014 Bob Spunt ----------------------
%	Created:  2014-09-29
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 2, disp('USAGE: resp = inputoption(prompt, respset)'); return; end
if iscell(prompt), prompt = char(prompt); end
numtag = 0; 
if any(isnumeric(respset))
    numtag = 1; 
    if size(respset, 1)==1, respset = respset'; end
    respset = cellstr(num2str(respset)); 
end
if ischar(respset), respset = cellstr(respset); end
nopt = length(respset);
respset = strtrim(respset); 
if nopt > 2
    optstr = sprintf(['(' repmat('%s, ', 1, nopt-1) 'or %s)'], respset{:});
else
    optstr = sprintf('(%s or %s)', respset{:}); 
end
resp = input(sprintf('%s %s ', prompt, optstr), 's'); 
while ~any(strcmpi(respset, resp))
    resp = input(sprintf('Invalid response. %s ', optstr), 's'); 
end
if numtag, resp = str2double(resp); end
end
function argstruct  = setargs(defaults, optargs)
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
function m = mfile_showhelp_contents
m = {                                                                                      
     'function mfile_showhelp(varargin)'                                                   
     '% MFILE_SHOWHELP'                                                                    
     'ST = dbstack(''-completenames'');'                                                   
     'if isempty(ST), fprintf(''\nYou must call this within a function\n\n''); return; end'
     'eval(sprintf(''help %s'', ST(2).file));'                                             
     'end'                                                                                 
     };    
end
function m = setargs_contents
m = {                                                                                                                              
     'function argstruct = setargs(defaults, optargs)'                                                                             
     '% SETARGS Name/value parsing and assignment of varargin with default values'                                                 
     'if nargin < 1, mfile_showhelp; return; end'                                                                                  
     'if nargin < 2, optargs = []; end'                                                                                            
     'defaults = reshape(defaults, 2, length(defaults)/2)''; '                                                                     
     'if ~isempty(optargs)'                                                                                                        
     '    if mod(length(optargs), 2)'                                                                                              
     '        error(''Optional inputs must be entered as Name, Value pairs, e.g., myfunction(''''name'''', value)''); '            
     '    end'                                                                                                                     
     '    arg = reshape(optargs, 2, length(optargs)/2)'';'                                                                         
     '    for i = 1:size(arg,1)'                                                                                                   
     '       idx = strncmpi(defaults(:,1), arg{i,1}, length(arg{i,1}));'                                                           
     '       if sum(idx) > 1'                                                                                                      
     '           error([''Input "%s" matches multiple valid inputs:'' repmat(''  %s'', 1, sum(idx))], arg{i,1}, defaults{idx, 1});'
     '       elseif ~any(idx)'                                                                                                     
     '           error(''Input "%s" does not match a valid input.'', arg{i,1});'                                                   
     '       else'                                                                                                                 
     '           defaults{idx,2} = arg{i,2};'                                                                                      
     '       end  '                                                                                                                
     '    end'                                                                                                                     
     'end'                                                                                                                         
     'for i = 1:size(defaults,1), assignin(''caller'', defaults{i,1}, defaults{i,2}); end'                                         
     'if nargout>0, argstruct = cell2struct(defaults(:,2), defaults(:,1)); end'                                                    
     'end'                                                                                                                         
     };  
end

 
 
 
 
 
 
 
 
