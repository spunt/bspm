function set_default(var_name, default)
% set_default -- Set default value for a variable in workspace.
%
%****f* wmtsa.utils/set_default
%
% NAME
%   set_default -- Set default value for a variable in workspace.
%
% SYNOPSIS
%   set_default(var_name, default)
%
% INPUTS
%   * var_name   -- Name of variable to check and set value (string).
%   * default    -- Value to set (various datatypes).
%
% OUTPUTS
%   (none)
%
% SIDE EFFECTS
%
%
% DESCRIPTION
%   set_defaults checks whether a variable with the name 'var_name' exists
%   in the caller's workspace.  If the variable  does not exist, the function
%   creates a variable the 'var_name' name and sets its value to 'default.  
%   If the variable does exist, the function leaves the variable value unchanged.
%
% ERRORS
%
%
% EXAMPLE
%   set_default('x', 1);
%
% NOTES
%   set_default creates and sets the default value for a single variable.
%   The set_defaults function creates and sets the default values for multiple
%   variables.
%
% SEE ALSO
%   set_defaults
%
% TOOLBOX
%   wmtsa/utils
%
% CATEGORY
%
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2005-08-01
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

%   $Id: set_default.m 612 2005-10-28 21:42:24Z ccornish $

%% Set Defaults
  
usage_str = [mfilename, '(var_name, default)'];
  
%% Check arguments.
error(nargerr(mfilename, nargin, [2:2], nargout, [0:0], 1, usage_str, 'struct'));

evalstr = ['~exist(''', var_name, ''', ''var'') || isempty(', var_name, ')'];
do_set_default = evalin('caller', evalstr);

if (do_set_default)
  var_value = default;
  assignin('caller', var_name, var_value);
end
  
return

