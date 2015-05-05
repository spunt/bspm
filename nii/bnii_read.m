function [img, hdr] = bnii_read(fname, varargin)
% BNII_READ Uses NII_TOOL to read NifTI images
%
%  USAGE: [img, hdr] = bnii_read(fname, varargin)  
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-04-03
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
def = { ... 
	'reshapeflag',		0     ...
	};
vals = setargs(def, varargin);
if nargin < 1, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
if ischar(fname), fname = cellstr(fname); end
if length(fname) > 1
    nii     = nii_tool('cat3D', fname); 
else
    nii     = nii_tool('load', char(fname)); 
end
img = nii.img;
hdr = nii.hdr;
if reshapeflag, img = reshape(img, prod(hdr.dim(2:4)), hdr.dim(5)); end


