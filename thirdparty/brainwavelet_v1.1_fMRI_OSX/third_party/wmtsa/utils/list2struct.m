function [s] = list2struct(list)
% liststruct -- Convert a list (cell array) of name-value pairs to a struct.
%
%****f* wmtsa.utils/list2struct
%
% NAME
%   cell2struct -- Convert a list (cell array) of name-value pairs to a struct.
%
% USAGE
%   [s] = list2struct(list)
%
% INPUTS
%   * list         -- list of name-value pairs (cell array)
%
% OUTPUTS
%   * s            -- structure fieldnames with values (struct).
%
%
% DESCRIPTION
%   list2struct converts a list (cell array) of name-value pair entries
%   into a structure with field names (name) and field values (value).
%
% EXAMPLE
%
%
% WARNINGS
%
%
% ERRORS
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
%   wmtsa/utils
%
% CATEGORY
%
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2005-07-18
%
% COPYRIGHT
%   (c) 2005 Charles R. Cornish
%
% CREDITS
%
% MATLAB VERSION
%   7.0
%
% REVISION
%   $Revision: 612 $
%
%***

%   $Id: list2struct.m 612 2005-10-28 21:42:24Z ccornish $

usage_str = ['Usage:  [s] = ', mfilename, ...
             '(list)'];
  
error(nargerr(mfilename, nargin, [1:1], nargout, [0:1], 1, usage_str, 'struct'));

if (~iscell(list) && wmtsa_isvector(list))
  error('list must be a cell array vector.');
end

fields = {};
values = {};


% Check list is a cell array of alternating name value pairs:

%  Ascertain that odd members of a cell array vector are valid fieldnames
%    (i.e. strings).
if (~iscellstr(list(1:2:length(list))))
  error('WMTSA:list2struct:invalidFieldName', ...
        'Names of name-value pairs in list must be strings.');
end

% Ascertain that length of list is even.
if (mod(length(list),2) ~= 0)
  error('WMTSA:list2strut:invalidListLength', ...
        'For a list of name-value pairs, length(list) must be even.');
end

% Convert into a lists of names and values;
fields = list(1:2:length(list));
values = list(2:2:length(list));

s = [];

for (i = 1:length(fields))
  field = fields{i};
  value = values{i};
  if (~ischar(field))
    error('WMTSA:list2struct:invalidFieldname', ...
          'Fieldname must be of type character.');
  end
  s.(field) = value;
end

return
