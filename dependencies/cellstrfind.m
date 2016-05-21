function idx = cellstrfind(cellarray, string)
% CELLSTRFIND Find cells containing a string
%
% USAGE: idx = cellstrfind(cellarray, string)
%
% ARGUMENTS:
%   cellarray:   cell array of strings
%   string:      string to locate
%
% ==============================================
if nargin<2, disp('USAGE: idx = cellstrfind(cellarray, string)'); return; end
if ~iscell(cellarray),
    disp('You should probably be using STRFIND.');
    return
end
idx = find(~cellfun('isempty',strfind(cellarray,string)));
