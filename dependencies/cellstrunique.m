function out = cellstrunique(in, delimiter)
% CELLSTRUNIQUE Find unique parts of each string in a cell array
%
% USAGE: out = cellstrunique(in, delimiter)
%
%       in: cell array of strings
%       delimiter: for splitting strings (default = '/')
%
% ===================================================================
if nargin<1, error('USAGE: out = cellstrunique(in, delimiter)'); end
if nargin<2, delimiter = filesep; end
if ischar(in), in = cellstr(in); end
alls = [];
for i = 1:length(in)
    s = regexp(in{i}, delimiter);
    alls = [alls; s];
end
for i = 1:size(alls,2)
    u(i) = length(unique(alls(:,i)));
end
out = alls(:,u>1);

