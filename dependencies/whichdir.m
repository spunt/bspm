function directory = whichdir(file)
% WHICHDIR Return directory of file in Matlab path
%
% USAGE: directory = whichdir(file)
%
% ARGUMENTS
%
%   file = file to find directory for
%
% ==================================================================
if nargin<1
    disp('USAGE: directory = whichdir(file)');
    return
end
if iscell(file)
    file = char(file);
end
d = which(file);
if isempty(d)
    disp(sprintf('%s not found in path!', file));
    directory = '';
else
    [directory n e] = fileparts(d);
end


