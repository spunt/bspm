function set_defaults(defaults, var_names)
%  set_defaults -- Create variables with default values if non-existant in workspace.
%
%****f* wmtsa.utils/set_defaults
%
% NAME
%    set_defaults -- Create variables with default values if non-existant in workspace.
%
% SYNOPSIS
%   set_defaults(defaults, [var_names])
%
% INPUTS
%   * defaults   -- strut of name-value pairs of defaults (struct).
%   * var_names  -- (optional) names of specific variables to set defaults
%                   (string or cell array of strings).
%
% OUTPUTS
%
% SIDE EFFECTS
%   Function call requires a minimum of 1 input arguments; otherwise error.
%
% DESCRIPTION
%   set_defaults checks whether variable(s) exist in the caller's workspace,
%   and, if not, creates them with the supplied default values.  The default
%   values are passed as a set of name-value pairs via the 'defaults' struct
%   argument.  If the 'var_names' argument is specified, then only those 
%   variables in var_names are created.  Other all variables with the
%   fieldnames in 'defaults' are created.
%
% USAGE
%
%
% WARNINGS
%
%
% ERRORS
%
%
% EXAMPLE
%   defaults.a = 1;
%   defaults.b = 'abc';
%   defaults.c = {'xyz', 1, 2, 3};
%   
%   % Example 1 - Create all variables with defaults;
%   %             Variables do not exist in caller's workspace.              
%   clear a b c
%   set_defaults(defaults);
%   whos a b c
%
%   % Example 2 - Create those variables that do not exist in caller's workspace.
%   clear a b c
%   b = 2
%   set_defaults(defaults);
%   whos a b c
%   b
%
%   % Example 3 - Create specific variables with defaults.
%   %             Some variables do exist in caller's workspace.              
%   clear a b c
%   a = 9
%   var_names = {'a', 'b'}
%   set_defaults(defaults);
%   whos a b c
% 
% NOTES
%
%
% BUGS
%
%
% TODO
%
%
% ALGORITHM
%
%
% REFERENCES
%
%
% SEE ALSO
%
%
% TOOLBOX
%
%
% CATEGORY
%
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%
%
% COPYRIGHT
%   (c) 2005 Charles R. Cornish
%
% MATLAB VERSION
%   7.0
%
% CREDITS
%
%
% REVISION
%   $Revision: 612 $
%
%***

%   $Id: set_defaults.m 612 2005-10-28 21:42:24Z ccornish $

%% Set Defaults

usage_str = [ mfilename, '(defaults, [var_names])'];
  
%% Check arguments.
error(nargerr(mfilename, nargin, [1:2], nargout, [0:0], 1, usage_str, 'struct'));

if (~isstruct(defaults))
  error('WMTSA:set_defaults:invalidArgumentDataType', ...
        'Argument (defaults) must be a struct.');
end


% If var_names is emtpy, set defaults for all variable in defaults struct.
if (~exist('var_names', 'var') || isempty(var_names))
  var_names = fieldnames(defaults);
elseif (ischar(var_names))
  var_names = {var_names};
else
%  if (~iscellstr(var_names) & ~wmtsa_isvector(var_names))
  if (~iscellstr(var_names) & ~wmtsa_isvector(var_names))
    error('WMTSA:set_defaults:invalidArgumentDataType', ...
          'Argument (var_names) must be a string or cell array of strings.');
  end
end

%% Set defaults for each var in var_names
for (i = 1:length(var_names))
  var_name = var_names{i};
  evalstr = ['~exist(''', var_name, ''', ''var'') || isempty(', var_name, ')'];
  do_set_default = evalin('caller', evalstr);
  
  if (do_set_default)
    if (isfield(defaults, var_name))
      var_value = defaults.(var_name);
    else
      error('WMTSA:set_defaults:noDefaultValue', ...
            ['No default value in struct (', inputname(1), ') for variable (', var_name, ').']);
    end
    assignin('caller', var_name, var_value);
  end
end


return

