function [] = bspm_minT(in, outname)
% BSPM_MINT
%
%   USAGE: bspm_minT(in, outname)
%       
%       in  =  array of images to smooth (full path)
%       outname = name for output image
%

% --------------- Copyright (C) 2014 ---------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<2, disp('USAGE: bspm_minT(in, outname)'); return, end

% make sure image names are cell arrays of strings
if ischar(in)
    in = cellstr(in);
end

% read data in from images
for i = 1:length(in)
    hdr{i} = spm_vol(in{i});
    im(:,:,:,i) = spm_read_vols(hdr{i});
end

% get min at every voxel
mim = min(im,[],4);

% construct outname and write image
h = hdr{1};
h.fname = outname;
tmp = h.descrip;
tmp(regexp(tmp,'-')+2:end) = [];
h.descrip = [tmp 'Minimum Statistic Image'];
spm_write_vol(h,mim);


 
 
 
 
