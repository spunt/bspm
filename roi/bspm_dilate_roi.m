function bspm_dilate_roi(in, ndilation, appendextent)
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

if nargin < 1, mfile_showhelp; return; end
if nargin < 2, ndilation = 1; end
if nargin < 3, appendextent = 0; end
if ischar(in), in = cellstr(in); end

for i = 1:length(in)
    
    %% read in image 
    hdr = spm_vol(in{i});
    img = spm_read_vols(hdr);

    %% dilate 
    kernel = cat(3,[0 0 0; 0 1 0; 0 0 0],[0 1 0; 1 1 1; 0 1 0],[0 0 0; 0 1 0; 0 0 0]);
    for s = 1:ndilation, img = spm_dilate(img, kernel); end

    %% get extent
    extent = sum(img(:)>0);
    
    %% change name
    oldname = hdr.fname;
    [p, n] = fileparts(oldname);
    if appendextent
       hdr.fname = fullfile(p, sprintf('%s_D%d_k=%d.nii', n, ndilation, extent)); 
    else
       hdr.fname = fullfile(p, sprintf('%s_D%d.nii', n, ndilation)); 
    end

    %% write
    spm_write_vol(hdr, img);
    
end

 
 
 
 
