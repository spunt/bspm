function [err] = validate_opts(opts, valid_opts, return_type)
% validate_opts -- Validate opts fieldnames.
%
%****f* wmtsa.utils/validate_opts
%
% NAME
%   validate_opts -- Validate opts fieldnames.
%
% SYNOPSIS
%   [errmsg] = validate_opts(opts, valid_opts)
%   [errmsg] = validate_opts(opts, valid_opts, 'string')
%   [errstruct] = validate_opts(opts, valid_opts, 'struct')
%
% INPUTS
%   * opts       -- name-value pairs to validate (struct).
%   * valid_opts -- valid opts names (struct or cell array of strings).
%
% OUTPUTS
%   * errmsg        -- error message (string).
%   * errstruct     -- error struct with fields:  message, identifier.
%
% SIDE EFFECTS
%   Function call requires a minimum of 2 input arguments; otherwise error.
%
% DESCRIPTION
%   validate_opts validates the fieldnames in the struct 'opts' against a
%   list of possible option names contained in valid_opts.  valid_opts may
%   be another struct or a cell array of strings.
%
%   If all fieldnames are valid, the function returns with empty output 
%   arguments.  If the function encounters an opt fieldname not found in
%   the set of valid fieldnames, it returns an error message encoding
%   the name of invalid opt fieldname.
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
%
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

%   $Id: validate_opts.m 612 2005-10-28 21:42:24Z ccornish $

%% Set Defaults
valid_return_type = {'string', 'struct'};
  
usage_str = ['[err] = ', mfilename, '(opts, valid_opts, [return_type])'];
  
%% Check arguments.
error(nargerr(mfilename, nargin, [2:3], nargout, [0:1], 1, usage_str, 'struct'));

if (~exist('return_type', 'var') || isempty(return_type))
  return_type = 'string';
end
if (isempty(strmatch(return_type, valid_return_type)))
  err_s.message = ['Invalid value for input argument ', 'return_type'];
  err_s.identifier = ['WMTSA:nargerr:invalidInputArgumentValue'];
  error(err_s);
end

%% Nothing to parse if opts is empty.
if (isempty(opts))
  return
end

if (isstruct(valid_opts))
  valid_names = fieldnames(valid_opts);
elseif (iscellstr(valid_opts))
  valid_names = valid_opts;
else
  error('WMTSA:invalidArgumentDataType', ...
        'Argument (valid_opts) must be a struct or cell array of strings');
end

names = fieldnames(opts);

err = [];

for (i = 1:length(names))
  name = names{i};
  if (isempty(strmatch(name, valid_names, 'exact')))
    err_msg = [name, ' is not a valid option.'];
    if (strcmp(return_type, 'string'))
      err = err_msg;
    elseif (strcmp(return_type, 'struct'))
      err.message = err_msg;
      err.identifier = ['WMTSA:validate_opts:invalidOption'];
    end
  end
end



return
