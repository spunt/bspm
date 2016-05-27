function [flag, volinfo] = bspm_check_orientations(images, verbose)
% BSPM_CHECK_ORIENTATIONS Compare dimensions and orientations across images
% 
%   USAGE [flag, volinfo] = bspm_check_orientations(images, [verbose])
%
%       FLAG    0: images have same orientations, dimensions, and voxel sizes
%               1: images have different dimensions
%               2: images have different orientations and/or voxel sizes
%               3: images have different dimensions and different
%               orientations and/or voxel sizes
%
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_check_orientations.m 5097 2012-12-06 16:08:16Z guillaume $
if nargin < 1, mfile_showhelp; return; end
if nargin < 2, verbose = true; end
if ~isstruct(images)
    if iscell(images), images = char(images); end
    try images = spm_vol(images); catch err, rethrow err; end
end
flag = 0;
volinfo.fname   = {images.fname}'; 
volinfo.dims    = cat(1,images.dim);
volinfo.mats    = reshape(cat(3,images.mat),[16,numel(images)]); 
if any(any(diff(volinfo.dims,1,1),1))
    flag = flag + 1;
    volinfo.uniquedims = unique(volinfo.dims, 'rows'); 
    if verbose, disp('The images do not all have the same dimensions.'); end 
end
if any(any(abs(diff(volinfo.mats,1,2))>1e-4))
    flag = flag + 2; 
    volinfo.uniquemats = unique(volinfo.mats', 'rows');
    if verbose, disp('The images do not all have same orientation and/or voxel sizes.'); end
end
end
function mfile_showhelp(varargin)
    % MFILE_SHOWHELP
    ST = dbstack('-completenames');
    if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
    eval(sprintf('help %s', ST(2).file));  
end