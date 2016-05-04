function info = bspm_check_vols(in, minnvalid)
% BSPM_CHECK_VOLS Check volumes for common problems
%
%  USAGE: info = bspm_check_vols(in, [minnvalid])
%

% ---------------------- Copyright (C) 2014 Bob Spunt ----------------------
%	Created:  2014-11-26
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 1, mfile_showhelp; return; end
if nargin < 2, minnvalid = 0; end
if iscell(in), in = char(in); end
info.fname  = in; 
V           = spm_vol(in);
info.nvol   = length(V);
if any(any(diff(cat(1, V.dim),1,1),1)), disp('The images do not all have the same dimensions.'); return; end
if any(any(abs(diff(reshape(cat(3,V.mat),[16,numel(V)]),1,2))>1e-4)), disp('The images do not all have same orientation.'); return; end
info.dim    = V(1).dim; 
info.mat    = V(1).mat; 
data        = spm_read_vols(V);
data        = reshape(data, prod(info.dim(1:3)), info.nvol);
data(data==0) = NaN;
info.nvox   = size(data,1); 
info.nvalidbysub = sum(~isnan(data)); 
info.nvalidbyvox = sum(~isnan(data), 2);
info.nvalidall   = sum(info.nvalidbyvox==info.nvol);
if minnvalid
    hdr = V(1); 
    hdr.fname = sprintf('mask_min%dsubs.nii', minnvalid); min
    hdr.descrip = sprintf('Mask - Voxels present in at least %d subjects', minnvalid); 
    mask = double(info.nvalidbyvox >= minnvalid); 
    spm_write_vol(hdr, reshape(mask, info.dim)); 
end

end
