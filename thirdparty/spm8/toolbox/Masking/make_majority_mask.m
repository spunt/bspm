function make_majority_mask(thr, cons, bfn, files, userglobs)
%  make_majority_mask(thr, cons, bfn, files, userglobs)
%
% Create mask based on consensus of cons images being above threshold thr.
% cons may be specified as a number of images or as a fraction (<=1) of
% the total number later selected (cons=1 will be treated as 100%, while
% cons=1.1 will be rounded down to 1 image). Either/both thr and cons may
% be vectors, in which case all their combinations will be evaluated.
%
% If any thresholds are negative, their absolute values will be used to
% derive thresholds *relative* to each image's total GM volume in litres.
% Or relative to the elements of the userglobs vector if this is specified.
% Note that this differs from relative thresholding in SPM, which is based
% on spm_global's mean(mean/8) heuristic. To get that use imaginary thr(s).
%
% Mask images are saved with filename bfn (suffixed with cons and thr),
% bfn can include a path, which defaults to the current dir.
%
% If 'files' is a single SPM.mat file, its list of scans will be reused.
% If files is empty or unspecified, interactive GUI selection will be used.
%
% Examples:
%
%   make_majority_mask; % no arguments
%   % default to an absolute threshold of 0.1, a consensus of 70%,
%   % a basename of 'mjr_mask', with GUI selection of input images.
%
%   make_majority_mask([0.05 0.8i -0.4], [0.5 0.7 1], bfn, files)
%   % This will use an absolute threshold of 0.05, and thresholds of .8
%   % and .4 respectively, relative to the "global" mean value and total in
%   % litres. Consensus fractions of 50, 70, and 100 percent will be used,
%   % giving  a total of 9 mask options. thr0.05_cons1 and thr0.8i_cons1
%   % should match the equivalents from the usual SPM estimation.
%
% Regarding "implicit masking" (with NaN or zero), note that (nan > thr) is
% always false, so any non-zero threshold will include implicit masking.
%
% See also: opt_thresh, make_fslvbm_mask, make_const_mask
%
% Please reference: Ridgway et al. (in press) Neuroimage
%  http://dx.doi.org/10.1016/j.neuroimage.2008.08.045
%  pre-print available from http://eprints.ucl.ac.uk/13060/
%
% Email ged.ridgway@gmail.com with any questions
% http://www.cs.ucl.ac.uk/staff/g.ridgway/masking

% Note that SPM's estimation process also automatically masks out voxels
% where the data are constant. This means it is possible for the resulting
% mask.img (in the SPM results directory) to exclude some voxels present in
% the explicit mask. However, this is likely to be rare if the explicit
% mask is built using this script with any non-zero threshold, since voxels
% are only typically constant if they are zero in all scans. The script
% make_const_mask may be useful to see this effect in advance, though note
% that if you are building a mask from different scans to those you will
% test (e.g. a mask from a balanced subset, or one mask used for several
% techniques) you probably don't want to use this, as the test data might
% not be constant where the mask data was.

[spm5 select fparts] = check_spm_setup;

% Default arguments
if ~exist('thr','var') || isempty(thr)
    thr = 0.1; disp(['using default threshold ' num2str(thr)]);
end
if ~exist('cons','var') || isempty(cons)
    cons = 0.7; disp(['using default consensus fraction ' num2str(cons)]);
end
if ~exist('bfn', 'var') || isempty(bfn)
    bfn = 'mjr_mask';
end
[pth fnm ext] = fparts(bfn);
if isempty(ext);
    ext = '.img';
end
bfn = fullfile(pth, fnm); % (add ext back in create_masks)

thr = thr(:); cons = cons(:); % force column vectors (or scalars)

% Get filenames or SPM.mat
if ~exist('files', 'var') || isempty(files)
    files = select('Select images/masks or SPM.mat'); drawnow
end
if size(files, 1) == 1 && ~isempty(regexp(files, 'SPM.mat$', 'once' ))
    load(files); % load SPM.mat
    files = char(SPM.xY.P);
end
if size(files, 1) < 1, error('No files specified'); end
vols = spm_vol(files);
if spm5
    spm_check_orientations(vols); % (Could resample vols into space of 1st)
end

% Setup output volume
MV = vols(1);
MV.pinfo = [1 0 0]';
if spm5, MV.dt(1)  = spm_type('uint8');
else     MV.dim(4) = spm_type('uint8');
end

% Check userglobs if specified
if exist('userglobs', 'var')
    userglobs = userglobs(:);
    if length(userglobs) ~= size(vols, 1);
        error('size(userglobs) should match size(files)');
    elseif ~any(thr < 0)
        warning('make_majority_mask:unusedglobs', ...
            'userglobs vector specified but ignored, since no rel. thr(s)')
    end
else userglobs = [];
end

% Process images, recording information needed to build masks
num_above = compute_num_above(thr, vols, userglobs);

% Build masks and write to disk
create_masks(thr, cons, bfn, ext, vols, num_above, MV);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function num_above = compute_num_above(thr, vols, userglobs)
% Count supra-threshold voxels, i.e. each voxel in num_above is the number
% of images which exceeded the threshold at that voxel.

T = length(thr);
num_above = cell(T, 1);
for t = 1:T
    num_above{t} = zeros(vols(1).dim);
end
% Loop over images
for n = 1:size(vols, 1)
    img = spm_read_vols(vols(n));
    img(isnan(img)) = 0; % treat NaN as zero, leave +/- Inf (good idea?)
    Thr = thr; % copy so can override with per-image relative thresholds...
    if any(Thr < 0)
        if ~isempty(userglobs)
            glob = userglobs(n);
        else
            vxvol = abs(det(vols(n).mat(1:3,1:3))) / 1e6; % voxel volume
            glob = sum(img(:)) * vxvol; % total GM volume (in litres)
        end
        Thr(Thr < 0) = abs(Thr(Thr < 0)) * glob;
    end
    if any(imag(Thr)~=0)
        glob = spm_global(vols(n));
        Thr(imag(Thr)~=0) = abs(Thr(imag(Thr)~=0)) * glob;
    end
    for t = 1:T
        num_above{t} = num_above{t} + ( img > Thr(t) );
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_masks(thr, cons, bfn, ext, vols, num_above, maskvol)
N = size(vols, 1);
% Ensure consensus is expressed as number of images, assuming
% fraction rather than num of images if cons <= 1
consN = zeros(size(cons));
consN(cons <= 1) = max(round(cons(cons <= 1)*N), 1);
consN(cons > 1)  = min(round(cons(cons > 1)), N);

% Create masks
for t = 1:length(thr)
    if ~isreal(thr(t))
        thrstr = sprintf('%gi', abs(thr(t)));
    else
        thrstr = sprintf('%g', thr(t)); % (handles negative okay)
    end
    for c = 1:length(consN)
        mask = num_above{t} >= consN(c);
        MV = maskvol;
        MV.fname = sprintf('%s_thr%s_cons%g%s', bfn, thrstr, cons(c), ext);
        spm_write_vol(MV, mask);
    end
end
