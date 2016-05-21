function [parpath, branchpar] = parentpath(subpaths)
% PARENTPATH Find parent path from multiple subpaths
%
%   USAGE:      parpath = parentpath(subpaths)
% __________________________________________________________________________
%   SUBPATHS:   CHAR or CELL array containing strings for multiple paths
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-03-27
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if      nargin < 1, disp('USAGE: parpath = parentpath(subpaths)'); return; end
if      iscell(subpaths), subpaths = char(subpaths); end
if      size(subpaths, 1)==1
    parpath = fileparts(subpaths); 
    disp('Only one subpath!'); 
    return; 
end

% | Get indices of noncommon characters
if      size(subpaths, 1)==2, diffidx = find(diff(subpaths));
else    diffidx = find(sum(diff(subpaths))); end

% | Assign parent path
parpath = fileparts(subpaths(1, 1:diffidx(1)));

if nargout==2
   tmp = regexp(regexprep(cellstr(subpaths), parpath, ''), filesep, 'split');
   branchpar = cellfun(@(x) x(2), tmp); 
end

end

