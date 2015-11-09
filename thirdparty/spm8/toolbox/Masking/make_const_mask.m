function make_const_mask(ims, out, nans)
%make_const_mask(ims, out) -- generate mask excluding all constant voxels
% i.e. only voxels which changed over the scans will be 1; the rest 0.
%
% NB, automatically done in SPM's estimation machinery. This script is
% intended to let you see the effects in advance, and/or to exactly
% reproduce SPM's behaviour when using other software.
%
% Voxels with some NaN values and some changing values will be included,
% unless nans is false.
%
% See also: make_majority_mask, opt_thresh

[spm5 select fparts] = check_spm_setup;

if ~exist('ims', 'var') || isempty(ims)
    ims = spm_select(inf, 'image', 'Select unsmth unmodul warped segs');
end
if ~exist('out', 'var') || isempty(out)
    out = 'fslvbm_mask.img'; % (create in CWD)
end
if ~exist('nans', 'var') || isempty(nans)
    nans = true;
end

Vi = spm_vol(ims);

Vo = struct(...
    'fname',   out,         ...
    'mat',     Vi(1).mat,   ...
    'descrip', 'mask image' ...
    );
if spm5
    Vo.dim = Vi(1).dim(1:3);
    Vo.dt = [spm_type('uint8') Vi(1).dt(2)];
else
    Vo.dim = [Vi(1).dim(1:3) spm_type('uint8')];
end

interp = 0; % nearest neighbour -- should have identical dimensions anyway
flags = {true, false, interp}; % (dmtx, mask, hold)

if nans
    expr = 'any(diff(X))';
else
    expr = 'any(diff(X)) & ~any(isnan(X))';
end

spm_imcalc(Vi, Vo, expr, flags);
