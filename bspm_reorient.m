function bspm_reorient(images, reference)
% BSPM_REORIENT  Wrapper for spm_get_space
%
% USAGE: bspm_reorient(images, reference)
%
% ARGUMENTS
%   images: image or cell array of images 
%   reference: image from which to pull transformation matrix
%       NOTE: images and reference must all be same size!
%

% -------------------- Copyright (C) 2014 --------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<2, disp('USAGE: bspm_reorient(images, reference)'); return; end
if ~iscell(images), images = cellstr(images); end
if iscell(reference), reference = char(reference); end

% check reference to ensure size is same as input images
h1 = spm_vol(reference); h2 = spm_vol(images{1});
if mean(h1.dim==h2.dim)~=1
    disp('Reference and input images do not have the same dimensions!'); 
    return; 
else
    refmat = h1.mat;
end

% apply to input images
for i = 1:length(images)
    spm_get_space(images{i},refmat);
end
 
 
 
 
