function [data, hdr] = bspm_read_vol_multi(in, reshapetag, implicitmasktag)
% BSPM_READ_VOL_MULTI  Wrapper for spm_vol/spm_read_vols + added options
%
% USAGE: [data, hdr] = bspm_read_vol_multi(in, reshapetag, implicitmasktag)
%
%   ARGUMENTS
%       in: volume filename
%       reshapetag: 
%           0: do not reshape (default)
%           1: reshape to 2D
%           2: reshape to 2D and zscore
%       implicitmasktag: 
%           0: no masking (default)
%           1: yes (NaNs voxels less than 10% of the mean)
%

% ---------------------------- Copyright (C) 2014 ----------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3, implicitmasktag = 0; end
if nargin<2, reshapetag = 0; end
if nargin<1, error('USAGE: [data, hdr] = bspm_read_vol_multi(in, reshapetag, implicitmasktag)'); end
if ischar(in), in = cellstr(in); end
hdr = spm_vol(char(in));
data = spm_read_vols(hdr);
if implicitmasktag, data(data < mean(data(:))/10) = NaN; end
if reshapetag 
    datadim = size(data);
    if length(datadim)==3, datadim(4) = 1; end
    data = reshape(data, prod(datadim(1:3)), datadim(4));
end
if reshapetag==2
    for i = 1:size(data,2), data(~isnan(data(:,i)),i) = zscore(data(~isnan(data(:,i)),i)); end
end
 
 
 
 
