function ts = bspm_filter_mat(ts, TR, cutoff)
% BSPM_FILTER
%
%   USAGE: bspm_filter(images, TR, cutoff, outprefix)
% 
%   images = volumes to filter
%   TR = repetition time (s)
%   cutoff = high-pass filter cutoff (s)
%   outprefix = prefix for output files
%

% ----------------- Copyright (C) 2014 -----------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 3, mfile_showhelp; return; end
dims = size(ts);

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
fprintf('\nApplying Temporal Filter: ');
for i = 1:length(maskidx);
    ts_rs_clean(:,maskidx(i)) = spm_filter(K,ts_rs(:,maskidx(i)));
end
fprintf('Complete!\n');
ts_clean = reshape(ts_rs_clean',dims);
 
 
 
 
