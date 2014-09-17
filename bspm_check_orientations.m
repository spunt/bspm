function [sts, info] = bspm_check_orientations(V, verbose)
% Check the dimensions and orientations of the images
% FORMAT [sts, str] = spm_check_orientations(V [,verbose])
% V       - a struct array as returned by spm_vol
% verbose - [Default: true]
%
% sts     - status (true means OK)
% str     - string describing status, empty if OK
%
% When used without LHS, this function throws an error accordingly.
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_check_orientations.m 5097 2012-12-06 16:08:16Z guillaume $
if nargin < 1, error('USAGE: [sts, str] = bspm_check_orientations(images, verbose)'); end
if nargin < 2, verbose = true; end
if ~isstruct(V)
    if iscell(V), V = char(V); end
    try V = spm_vol(V); catch err, rethrow err; end
end
sts = true;
if nargout==2
    info = [{V.fname}' {V.dim}' {V.mat}']; 
end
dims = cat(1,V.dim);
if any(any(diff(dims,1,1),1))
    sts = false;
    if verbose, disp('The images do not all have the same dimensions.'); end
    info.im 
end
matx = reshape(cat(3,V.mat),[16,numel(V)]);
if any(any(abs(diff(matx,1,2))>1e-4))
    sts = false;
    if verbose, disp('The images do not all have same orientation and/or voxel sizes.'); end
end
end
