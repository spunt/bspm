function make_fslvbm_mask(ims, out, maxthr, minthr)
%make_fslvbm_mask(ims, out) -- generate mask as FSL-4's fslvbm scripts do.
% Note that fslvbm used *unsmoothed* segmentations, and *unmodulated*
% normalised ones even when (smoothed) modulated normalised images
% are entered into the statistical analysis.
%
% fslvbm's mask includes voxels which meet *both* of these criteria:
%  - minimum GM prob over segs > minthr
%  - maximum GM prob over segs >= maxthr
% where maxthr is 100 for fsl's 1000-scaled ints, and hence 0.1 for SPM
%
%  make_fslvbm_mask(ims, out, maxthr, minthr) % allows different thresholds
%
% See also: make_majority_mask, opt_thresh

[spm5 select fparts] = check_spm_setup;

if ~exist('ims', 'var') || isempty(ims)
    ims = select('Select unsmth unmodul warped segs');
end
if ~exist('out', 'var') || isempty(out)
    out = 'fslvbm_mask.img'; % (create in CWD)
end
if ~exist('maxthr', 'var') || isempty(maxthr)
    maxthr = 0.1; % can override to e.g. 100 for use on fslvbm images
end
if ~exist('minthr', 'var') || isempty(minthr)
    minthr = 0;
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
expr = sprintf('min(X) > %g & max(X) >= %g', minthr, maxthr);

spm_imcalc(Vi, Vo, expr, flags);
