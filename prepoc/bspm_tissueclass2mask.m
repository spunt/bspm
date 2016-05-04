function bspm_tissueclass2mask(cimages, outname)
% USAGE: bspm_tissueclass2mask(cimages, outname)
%


% ------------------------- Copyright (C) 2014 -------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin < 1, mfile_showhelp; return; end
if ischar(cimages), cimages = cellstr(cimages); end
if nargin < 2, outname = fullfile(fileparts(cimages{1}), 'mask_brain.nii'); end
% | Read & Sum
[im,h] = bspm_read_vol(cimages);
im = nansum(im, 4);
% | Threshold & Fill Holes
im(im < graythresh(im(im>0))) = 0; 
im = imfill(im,6,'holes');
% | Write
im = double(im > 0); 
hout        = h(1);
hout.fname  = outname;
hout.descrip = 'Mask of Brain from Tissue Classes'; 
spm_write_vol(hout, im);

end
