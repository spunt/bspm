function [opts] = parse_opts(varargin)
% parse_opts -- Parse name-value pair options into a struct.
%
%****f* wmtsa.utils/parse_opts
%
% NAME
%   parse_opts -- Parse name-value pair options into a struct.
%
% SYNOPSIS
%   [opts] = parse_opts(varargin)
%
% INPUTS
%   * varargin   -- variable input argument list.
%
% OUTPUTS
%   * opts       -- struct of name-value pairs.
%
% SIDE EFFECTS
%   Function call requires a minimum of 2 input arguments; otherwise error.
%
% DESCRIPTION
%   parse_opts parses the input arguments for name-value pairs and returns the
%   struct 'opts' filled with name-value pairs.
%
%  Input arguments must be one of following:
%  * (1) varargin list of name/value pairs, or
%  * (2) a single argument of type struct with name/value pairs, or
%  * (3) a single argument even-length vector of type cell with name/value pairs.
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
%   % Example 1: Parse varargin 
%   opts = parse_opts('a', 1, 'b', 'xyz', 'c', {'abc', 'def', 'jkl'});
%
%   % Example 2: Parse a cell array
%   opts_list = {'a', 1, 'b', 'xyz', 'c', {'abc', 'def', 'jkl'}};
%    opts2 = parse_opts(opts_list);
%  
%   % Example 2: Parse a struct
%    opt3 = parse_opts(opts);
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
%   2005-07-19
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

%   $Id: parse_opts.m 612 2005-10-28 21:42:24Z ccornish $

%% Set Defaults
  
usage_str = ['[opts] = ', mfilename, '(varargin)'];
  
%% Check arguments.
error(nargerr(mfilename, nargin, ':', nargout, [0:1], 1, usage_str, 'struct'));

%% If argument is empty, create an empty opts struct
if (nargin == 1 & isempty(varargin{1}))
  opts = struct;
  return
end

%% Check if first and only argument is a struct.
if (nargin == 1 & isstruct(varargin{1}))
  opts = varargin{1};
  return
end

%% Check if varargin is passed as a single cell array; if so convert to cell array.
%if (nargin == 1 & iscell(varargin{1}))
if (nargin == 1 & iscell(varargin))
  list = varargin{:};
else
  list = varargin;
end

%% The list cell array is empty.
if (isempty(list))
  opts = struct;
  return
end

%% Parse if input is a vector cell array.
if (wmtsa_isvector(list) & iscell(list))
  try
    opts = list2struct(list);
  catch
    rethrow(lasterror);
  end
else
%% Otherwise throw error:  can only process simple cell array vectors.
  error('WMTSA:parse_opts:invalidInputFormat', ...
        ['Input must be one of following: ', ...
         '(1) varargin list of name/value pairs, or ', ...
         '(2) a single argument of type struct with name/value pairs, or', ...
         '(3) a single argument even-length vector of type cell with name/value pairs.']);
end

return
