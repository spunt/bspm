function bspm_dilate_roi(in, size)
% BSPM_DILATE_ROI
%
%   USAGE: bspm_dilate_roi(in, size)
%       
%       in  =  input image filename
%       size = size of dilation (default = 1);
%

% ------------ Copyright (C) 2014 ------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, disp('USAGE: bspm_dilate_roi(in, size)'); return; end
if nargin < 2, size = 1; end
if ischar(in), in = cellstr(in); end

for i = 1:length(in)
    
    %% read in image 
    hdr = spm_vol(in{i});
    img = spm_read_vols(hdr);

    %% dilate 
    kernel = cat(3,[0 0 0; 0 1 0; 0 0 0],[0 1 0; 1 1 1; 0 1 0],[0 0 0; 0 1 0; 0 0 0]);
    for s = 1:size, img = spm_dilate(img, kernel); end

    %% get extent
    extent = sum(img(:)>0);
    
    %% change name
    oldname = hdr.fname;
    [p n e] = fileparts(oldname);
    hdr.fname = [p filesep n '_D' num2str(size) '_k=' num2str(extent) '.nii'];

    %% write
    spm_write_vol(hdr, img);
    
end

 
 
 
 
