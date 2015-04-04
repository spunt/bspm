function tsinfo = bspm_badscan(epi, varargin)
% BSPM_BADSCAN Identifies bad scans and saves a nuisance regressor file
% 
%   USAGE: tsinfo = bspm_badscan(epi, varargin)
%
%  ARGUMENTS
%    epi = cell array of EPI filenames
%   'dvars_thresh',     2.5,  
%   'framewise_thresh', 0.5,  
%   'include_rp',       1,    
%   'maskfile',         []
%   
% 
% Bob Spunt, Caltech
% Based on a script by Donald McLaren, Ph.D.
% 2012_05_01 -- Created
% 2012_06_04 -- Added to FUNC + Presented in SCAN Lab Workshop
% 2013_03_08 -- Added option to MASK data prior to computing bad scans
% ========================================================================%
def = { 'dvars_thresh',     2.5,  ...
        'framewise_thresh', 0.5,  ...
        'include_rp',       1,    ...
        'maskfile',         []};
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end 

% | create output filename
if ischar(epi), epi = cellstr(epi); end
outfile     = sprintf('badscan_dvars%dframewise%d_%s.txt', dvars_thresh*100, framewise_thresh*100, strtrim(datestr(now,'mmm_DD_YYYY')));
nvol        = length(epi);
if nvol==1
   epi      = bspm_expand4D(epi); 
   nvol     = length(epi); 
end
epidir      = fileparts(epi{1});

% | motion parameters
rp                  = load(char(files(fullfile(epidir, 'rp_*txt'))));
rp(:,4:6)           = rp(:,4:6)*57.3;
framewise           = zeros(nvol,1);
framewise(2:end)    = max(abs(diff(rp)), [], 2);

% | get data and apply implicit masking
fprintf('\n | - Loading Data for %d Volumes', nvol); 
cfg.plot            = 0; 
if isempty(maskfile)
    cfg.vol     = bspm_read_vol(epi, 'implicit');
else
    cfg.vol     = bspm_read_vol(epi, 'mask', maskfile); 
end

% | use BRAMILA tools to compute dvars and framewise displacement
fprintf('\n | - Identifying Scans with DVARS > %2.2f and framewise displacement > %2.2f', dvars_thresh, framewise_thresh);
dvars       = bramila_dvars(cfg);

% | get indices of bad timepoints
zdvars          = oneoutzscore(dvars(2:end), 1);
badidx          = zeros(nvol, 2);
badidx(2:end,1) = zdvars > dvars_thresh; 
badidx(:,2)     = framewise > framewise_thresh;
badidx          = any(badidx, 2); 
tsinfo.pctbad   = round(100*(sum(badidx)/nvol));
fprintf('\n | - %d (%d%%) scans identified as bad', sum(badidx), tsinfo.pctbad);

% | create bad volume regressor
scrubmat = zeros(nvol, sum(badidx)); 
for i = 1:length(badidx), scrubmat(badidx(i),i) = 1; end
if include_rp==1, scrubmat = [rp scrubmat]; end

% | save as a new text file
save(fullfile(epidir, outfile), 'scrubmat', '-ascii');
fprintf('\n | - Nuisance regressors written to:  %s\n\n', outfile);

if nargout>0
    tsinfo.epidir       = epidir; 
    tsinfo.maxmotion    = max(rp) - min(rp);
    tsinfo.dvars        = dvars; 
    tsinfo.framewise    = framewise;
    tsinfo.scrubmat     = scrubmat; 
end

end
%---------------------------------------------------------------------%
% SUBFUNCTIONS
%---------------------------------------------------------------------%
function zx = oneoutzscore(x, returnas)
% ONEOUTZSCORE Perform columnwise leave-one-out zscoring
% 
% USAGE: zx = oneoutzscore(x, returnas)
% 
%   returnas: 0, signed values (default); 1, absolute values
%
if nargin<1, disp('USAGE: zx = oneoutzscore(x, returnas)'); return; end
if nargin<2, returnas = 0; end
if size(x,1)==1, x=x'; end
zx              = x; 
[nrow, ncol]    = size(x);
for c = 1:ncol
    cin         = repmat(x(:,c), 1, nrow);
    theoneout   = cin(logical(eye(nrow)))';
    theleftin   = reshape(cin(logical(~eye(nrow))),nrow-1,nrow);
    cz          = (theoneout-nanmean(theleftin))./nanstd(theleftin);
    zx(:,c)     = cz';
end
if returnas, zx = abs(zx); end
end
 
 
 
 
