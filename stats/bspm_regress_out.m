function bspm_regress_out(images, covariates, outprefix)
% BSPM_REGRESS_OUT
% 
%   USAGE: bspm_regress_out(images, covariates, outprefix)
%
%   images = volumes to filter
%   covariates = covariates to regress out
%   outprefix = prefix for output files
%

% -------------------- Copyright (C) 2014 --------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3, outprefix = 'r'; end
if nargin<2, display('USAGE: bspm_regress_out(images, covariates, outprefix)'); return; end
if iscell(images); images = char(images); end;
if iscell(covariates), covariates = char(covariates); end
if exist(covariates, 'file'), covariates = load(covariates); end

% make X matrix
X = [covariates ones(length(covariates),1)];

% load in raw timeseries
hdr = spm_vol(images); ts = spm_read_vols(hdr);
dims = size(ts);
ts_rs = reshape(ts,prod(dims(1:3)),dims(4))';
tmp = nanmean(ts_rs);
maskidx = find(tmp>mean(tmp)/8);
fprintf('\nRunning Voxelwise Regression: ');
ts_rs_clean = bspm_residuals(ts_rs, X);
fprintf('Complete!\n');
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
    
 
 
 
 
 
 
 
 
