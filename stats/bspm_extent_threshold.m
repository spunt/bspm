function out = bspm_extent_threshold(in, extent)
% BSPM_EXTENT_THRESHOLD
%
% USAGE: out = bspm_extent_threshold(in, extent)
%   
%   NOTE: "in" must be image data, not file
%

% --------------- Copyright (C) 2014 ---------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<2, mfile_showhelp; return; end
    
imdims = size(in);
[X Y Z] = ind2sub(size(in), find(in > 0));
voxels = sortrows([X Y Z])';
cl_index = spm_clusters(voxels);
cidx = unique(cl_index);
for i = cidx, cl_size(i) = sum(cl_index==cidx(i)); end
badidx = ismember(cl_index, find(cl_size<extent));
bad = voxels(:,badidx)';
out = in;
for i = 1:size(bad,1)
    out(bad(i,1),bad(i,2),bad(i,3)) = 0;
end





 
 
 
 
