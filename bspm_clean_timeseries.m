function bspm_clean_timeseries(images, TR, cutoff, covariates, outprefix) 
% BSPM_CLEAN_TIMESERIES
%
%   USAGE: bspm_filter(images, TR, cutoff, outprefix)
% 
%   images = volumes to filter
%   TR = repetition time (s)
%   cutoff = high-pass filter cutoff (s)
%   covariates = covariates to regress out
%   outprefix = prefix for output files
%

% ---------------------------- Copyright (C) 2014 ----------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<5, outprefix = 'f'; end
if nargin<4, display('bspm_clean_timeseries(images, TR, cutoff, covariates, outprefix)'); return; end
if iscell(images); images = char(images); end;

% load in raw timeseries
hdr = spm_vol(images); ts = spm_read_vols(hdr);
dims = size(ts);

% make X matrix
X = [covariates ones(length(covariates),1)];

% spm_filter stuff
K.RT = TR;
K.HParam = cutoff;
K.row = 1:dims(4);
K = spm_filter(K);

% back to the timeseries
ts_rs = reshape(ts,prod(dims(1:3)),dims(4))';
ts_rs_clean = zeros(size(ts_rs));
tmp = nanmean(ts_rs);
maskidx = find(tmp>mean(tmp)/8);
fprintf('\nCleaning Timeseries: ');
ts_rs(:,maskidx) = spm_filter(K,ts_rs(:,maskidx));
ts_rs_clean(:,maskidx) = bob_residuals(ts_rs(:,maskidx), covariates);
fprintf('Complete!\n');
% for i = 1:length(maskidx);
%     
%     ts_filt = spm_filter(K,ts_rs(:,maskidx(i)));
%     [b bint res] = regress(ts_filt,X);
%     ts_rs_clean(:,maskidx(i)) = nanmean(ts_filt) + res;
%     
% end
ts_clean = reshape(ts_rs_clean',dims);
fprintf('Writing New Volumes: ');
% write new vols
for i = 1:length(hdr)
    h = hdr(i);
    [p n e] = fileparts(h.fname);
    h.fname = [p filesep outprefix n '.nii'];
    spm_write_vol(h,ts_clean(:,:,:,i));
end
fprintf('Complete!\n\n');
    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
