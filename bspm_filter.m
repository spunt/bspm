function ts = bspm_filter(images, TR, cutoff, outprefix)
% BSPM_FILTER
%
%   USAGE: bspm_filter(images, TR, cutoff, outprefix)
% 
%   images = volumes to filter
%   TR = repetition time (s)
%   cutoff = high-pass filter cutoff (s)
%   outprefix = prefix for output files
%

% ------------------- Copyright (C) 2014 -------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<4, outprefix = 'f'; end
if nargin<3, display('USAGE: bspm_filter(images, TR, cutoff, outprefix)'); return; end
if iscell(images); images = char(images); end;

% load in raw timeseries
hdr = spm_vol(images); ts = spm_read_vols(hdr);
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
fprintf('Writing New Volumes: ');
% write new vols
for i = 1:length(hdr)
    h = hdr(i);
    [p n e] = fileparts(h.fname);
    h.fname = [p filesep outprefix n '.nii'];
    spm_write_vol(h,ts_clean(:,:,:,i));
end
fprintf('Complete!\n\n');
    
 
 
 
 
