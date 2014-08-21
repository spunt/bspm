function [data, hdr, info] = bspm_read_vol(in, varargin)
% BSPM_READ_VOL  Wrapper for spm_vol/spm_read_vols + added options
%
% USAGE: [data hdr info] = bspm_read_vols(in, varargin)
%
%   ARGUMENTS
%       in: filename(s) of image volumes
%
%   OPTIONAL ARGUMENTS
%       'reshape' - reshapes to 2D, rows = voxels, cols = images
%       'zscore'  - zscores the image data
%       'zscoretime'  - zscores the image data by time
%       'implicit'- NaNs voxels less than 10% of the mean
%       'mask' - (must pair with filename)
%       

% ------------------------ Copyright (C) 2014 ------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, error('USAGE: [data hdr info] = bspm_read_vols(in,option)'); end
if nargin > 1, optional = 1; else optional = 0; end
if iscell(in), in = char(in); end
hdr = spm_vol(in);
nvol = length(hdr);
data = spm_read_vols(hdr);
datadim = size(data);
if nargout==3
    info.dim = datadim;
    info.nvalid = sum(~isnan(data(:)));
    info.mean = nanmean(data(:));
    info.median = nanmedian(data(:));
    info.std = nanstd(data(:)); 
    info.max = nanmax(data(:));
    info.min = nanmin(data(:));
end
if optional
    varargin = lower(varargin);
    if ismember('implicit', varargin);
        for i = 1:nvol
            tmp = data(:,:,:,i);
            tmp(tmp < nanmean(tmp(:))/10) = NaN;
            data(:,:,:,i) = tmp;
        end
    end
    if ismember('mask', varargin)
        maskfile = varargin{find(ismember(varargin, 'mask'))+1};
        if iscell(maskfile), maskfile = char(maskfile); end
        mask = bob_reslice(maskfile, hdr(1).fname, 1, 1);
        if length(datadim)==3 datadim(4) = 1; end
        tmp = data;
        tmp = reshape(tmp, prod(datadim(1:3)), datadim(4)); 
        tmp(mask(:)~=1,:) = NaN;
        data = reshape(tmp, datadim);
    end
    if ismember('reshape', varargin);
        if length(datadim)==3 datadim(4) = 1; end
        data = reshape(data, prod(datadim(1:3)), datadim(4));
        if ismember('zscoretime', varargin);
            d = bsxfun(@minus, data', nanmean(data'));
            d = bsxfun(@rdivide, d, nanstd(d));
            data = d';
        end
        if ismember('zscore', varargin);
            for i = 1:size(data,2), data(~isnan(data(:,i)),i) = zscore(data(~isnan(data(:,i)),i)); end
        end
        return
    end
    if ismember('zscoretime', varargin);
        if nvol>1
            data = reshape(data, prod(datadim(1:3)), datadim(4));
            d = bsxfun(@minus, data', nanmean(data'));
            d = bsxfun(@rdivide, d, nanstd(d));
            data = d';
            data = reshape(data, datadim);
        end  
    end
end
 
 
 
 
