function out = clean_text(in, remove)
% CLEAN_GRIDIN
%
% USAGE: out = clean_text(in, remove)
%
%      in = cell array text
%
% -------------------------------------------------------------------
if nargin < 1, error('USAGE: out = clean_text(in, remove)'); end
if nargin < 2, remove = []; end
for i = 1:length(remove), in = regexprep(in, remove{i}, '', 'ignorecase'); end
in = strtrim(in);
in = regexprep(in, '\.', '');
in = regexprep(in, ' ', '_');
in = regexprep(in, '-','');
for i = 1:10
    in = regexprep(in, repmat('_',1,12-i), '_');
end
out = in;
