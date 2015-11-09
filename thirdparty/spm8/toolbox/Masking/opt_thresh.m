function [thr Vthr] = opt_thresh(in, optfunc, out)
%opt_thresh - determine an optimal threshold to binarise a continuous image
% [thr Vthr] = opt_thresh(input, optfunc, outname)
%
% If outname is included (even if empty), then the optimally thresholded
% volume is written to disk, and a volume handle is returned in Vthr
%
% optfunc is a function handle (or string containing one) that finds the
% optimal threshold. Two alternatives are included:
%
% (i)  Maximising correlation between original and binarised result
%      '@opt_thr_corr' (also default if optfunc is empty or unspecified)
% (ii) An estimate of the statistical anti-mode '@opt_thr_antimode'
%
% The first option is equivalent to maximising the t-value for a two-sample
% t-test between the below- and above-threshold parts of the image, which
% in turn is similar to the histogram-based criterion in Otsu's method.
%
% See also: make_average, make_majority_mask
%
% Please reference:
%  Ridgway et al. (in press) Neuroimage
%  http://dx.doi.org/10.1016/j.neuroimage.2008.08.045
%  pre-print available from http://eprints.ucl.ac.uk/13060/
%
% If using opt_thr_antimode, please also reference:
%  Luo and Nichols (2003) Neuroimage 19, 1014-1032
%  http://dx.doi.org/10.1016/S1053-8119(03)00149-6
%  pre-print available from http://www.sph.umich.edu/ni-stat/SPMd/SPMd.pdf
%
% Email ged.ridgway@gmail.com with any questions
% http://www.cs.ucl.ac.uk/staff/g.ridgway/masking

[spm5 select fparts] = check_spm_setup;

% Check if called from batch system
if nargin == 1 && isstruct(in) && isfield(in, 'inname')
    job     = in;
    in      = job.inname;
    optfunc = job.optfunc;
    [pth nam ext num] = spm_fileparts(job.outname);
    if isempty(pth)
        % {''} seems to end up as '' in SPM8, hence left hand side of ||
        if ~iscell(job.outdir) || isempty(job.outdir{1})
            pth = pwd;
        else
            pth = job.outdir{1};
        end
    end
    out = fullfile(pth, [nam ext num]);
end

if ~exist('in', 'var') || isempty(in)
    in = select('Select image'); drawnow
end
if iscellstr(in)
    in = char(in);
end
if isempty(in); error('No files selected'); end
Vi = spm_vol(in);

if ~exist('optfunc', 'var') || isempty(optfunc)
    optfunc = @opt_thr_corr;
end
if ~isa(optfunc, 'function_handle')
    % assume string, and try evaluating
    optfunc = eval(optfunc);
end

dat = spm_read_vols(Vi); dat = dat(:);

thr = optfunc(dat);

if ~exist('job', 'var') && nargin < 3
    return
elseif isempty(out)
    [pth bnm ext] = fparts(Vi.fname);
    suf = regexprep(sprintf('%g', thr), '\.', '-');
    out = fullfile(pth, [bnm '_thr' suf ext]);
elseif isempty(regexp(out, '\.(nii|img)$', 'once'))
    out = [out '.img'];
end

flags = {false, false, 0}; % (dmtx, mask, hold=nearest-neighbour)
expr = sprintf('i1>%g', thr);

Vo = struct(...
    'fname',   out,         ...
    'mat',     Vi.mat,      ...
    'descrip', 'mask image' ...
    );
if spm5
    Vo.dim = Vi.dim(1:3);
    Vo.dt = [spm_type('uint8') Vi.dt(2)];
else
    Vo.dim = [Vi.dim(1:3) spm_type('uint8')];
end

Vthr = spm_imcalc(Vi, Vo, expr, flags);

if exist('job', 'var') % called from batch system, need single output var
    out = struct('thr',thr, 'outname', {{Vthr.fname}}); thr = out;
end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thr = opt_thr_corr(img)
costfunc = @(thr) -correlation(img, img > thr);
[thr ncc] = fminbnd(costfunc, min(img), max(img));
fprintf('Maximal correlation of %g found with threshold of %g\n', ...
    -ncc, thr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = correlation(x, y)
cs = corrcoef(x, double(y));
c = cs(1, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thr = opt_thr_antimode(img) %#ok, optfunc can eval as handle to it
% See appendix B of Luo and Nichols (2003), Neuroimage 19, 1014-1032
% http://www.sph.umich.edu/ni-stat/SPMd/SPMd.pdf
% http://dx.doi.org/10.1016/S1053-8119(03)00149-6
% Implemented from the description, any mistakes are my fault - Ged

% Hartigan method
srt = sort(img);
lo = srt(floor(numel(srt) * 0.1));
hi = srt(ceil(numel(srt) * 0.9));
srt(srt <= lo) = [];
srt(srt >= hi) = [];
dfs = diff(srt);
mx = max(dfs);
k = floor(mean(find(dfs == mx))); % in case non-unique
thr = mean([srt(k) srt(k+1)]);
fprintf('Anti-mode threshold of %g by Hartigan method\n', thr);

% Histogram method
iqrange = diff(srt(round(numel(srt)*[0.25 0.75])));
binwidth = 1.595 * iqrange * numel(srt)^(-1/5);
numbins = ceil((srt(end) - srt(1)) / binwidth);
[counts bins] = hist(srt, numbins);
mn = min(counts);
k = round(mean(find(counts == mn))); % in case non-unique
thr2 = bins(k);
fprintf('Anti-mode threshold of %g by histogram method\n', thr2);

% Assume Hartigan more accurate if similar, but less reliable if different
if abs(thr - thr2) > 3 * binwidth
    thr = thr2;
end
fprintf('Anti-mode threshold of %g chosen.\n', thr);
