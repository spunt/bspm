function [img, hdr] = bnii_read(fname)
% BNII_READ Uses NII_TOOL to read NifTI images
%
%  USAGE: [img, hdr] = bnii_read(fname)  
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-04-03
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
if nargin==0, mfile_showhelp; return; end
if ischar(fname), fname = cellstr(fname); end
if length(fname) > 1
    nii     = nii_tool('cat3D', fname); 
else
    nii     = nii_tool('load', char(fname)); 
end
img = nii.img; 
hdr = nii.hdr; 


