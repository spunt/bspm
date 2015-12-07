function [out_opts] = set_opts_defaults(in_opts, opts_defaults, default_fields_only)
% set_opts_defaults -- Set default values for opts.
%
%****f* wmtsa.utils/set_opts_defaults
%
% NAME
%   set_opts_defaults -- Set default values for opts.
%
% SYNOPSIS
%   [opts] = set_opts_defaults(opts, opts_defaults)
%
% INPUTS
%   * in_opts          -- set of name-value pairs (struct).
%   * opts_defaults -- default values for opts (struct).
%   * default_fields_only -- return opts struct with default fields (logical)
%
% OUTPUTS
%   * out_opts          -- name-value pairs with default values (struct).
%
% SIDE EFFECTS
%
%
% DESCRIPTION
%
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
%   $Revision: 614 $
%
%***

%   $Id: set_opts_defaults.m 614 2006-04-27 21:43:40Z ccornish $

%% Set Defaults
  
usage_str = ['[out_opts] = ', mfilename, '(in_opts, opts_defaults, [default_fields_only])'];
  
%% Check arguments.
error(nargerr(mfilename, nargin, [2:3], nargout, [0:1], 1, usage_str, 'struct'));

set_default('default_fields_only', 0);

fields = fieldnames(opts_defaults);

if (default_fields_only)
  out_opts = [];
else
  out_opts = in_opts;
end

for (i = 1:length(fields))
  field = fields{i};
  if (isfield(in_opts, field))
    out_opts.(field) = in_opts.(field);
  else
    out_opts.(field) = opts_defaults.(field);
  end
end


return
