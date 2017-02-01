function overlap = bspm_save_count_map(images, heightThresh, sizeThresh, outname, pctflag)
% BSPM_SAVE_COUNT_MAP
%
% USAGE: bspm_save_count_map(images, heightThresh, sizeThresh, outname, pctflag)
%
%   ARGUMENTS
%       images:             input image name (full path if not in current directory)
%       heightThresh:   intensity threshold f(if < .10, will assume it is
%                               an alpha level and will covert to critical t)
%       sizeThresh:      extent threshold for defining clusters
%       outname:         name for file to write (default = no file written)
%       pctflag:            flag to write count as % of total N (default = 0)
%       
% Created January 1, 2013 - Bob Spunt

% ------------------------------------ Copyright (C) 2014 ------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<5, pctflag = 0; end
if nargin<4, mfile_showhelp; return; end

nim = length(images);
for i = 1:nim
    in = images{i};
    out = bspm_threshold_image(in,heightThresh,sizeThresh,1);
    overlap(:,:,:,i) = out;
end

% compute count
sum_overlap = sum(overlap,4);
if pctflag, sum_overlap = 100*(sum_overlap/nim); end

% write volume
h = spm_vol(images{i});
h.fname = outname;
h.descrip = sprintf('Count map thresholded at t > %2.3f and cluster size %d', heightThresh, sizeThresh);
spm_write_vol(h,sum_overlap);










 
 
 
 
