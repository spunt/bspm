function idx = cellismember(cellarray1, cellarray2)
% CELLISMEMBER
%
%   idx shows cells in cellarray1 that contain 
%   members of cellarray2
%
%   USAGE: idx = cellismember(cellarray1, cellarray2)
%
%   ARGUMENTS:
%       cellarray1:   array to search in
%       cellarray2:   array with members to search for
%
% ==============================================================
if nargin<2, error('USAGE: idx = cellismember(cellarray1, cellarray2)'); end
if ~iscell(cellarray1) || ~iscell(cellarray2), 
    error('Both inputs must be cell arrays!'); 
end

idx = [];
for i = 1:length(cellarray2)
    tmp = find(~cellfun('isempty',strfind(cellarray1,cellarray2{i})));
    if ~isempty(tmp), idx = [idx; tmp]; end
end
idx = sort(idx);
