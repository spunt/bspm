function [] = bspm_reslice(image_defining_space, images_to_reslice, interpolation_method, prefix)
% BSPM_RESLICE  Just a wrapper for re-slicing images
%
% USAGE: bspm_reslice(image_defining_space, images_to_reslice, interpolation_method, prefix)
% 
% Images should be defined in cell array of string(s). 
% Options for "interpolation_method" are:
%
%   0 - Nearest neighbor
%   1 - Trilinear
%   2:7 - 2nd through 7th degree B-Spline
%

% --------------------------------------- Copyright (C) 2014 ---------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<4, prefix = 'r'; end
if nargin<3, interpolation_method = 1; end
if nargin<2, disp('USAGE: bspm_reslice(image_defining_space, images_to_reslice)'); return; end
if ischar(image_defining_space); image_defining_space = cellstr(image_defining_space); end
if ischar(images_to_reslice); images_to_reslice = cellstr(images_to_reslice); end

% run job 
spm('defaults','fmri'); spm_jobman('initcfg');         
matlabbatch{1}.spm.spatial.coreg.write.ref = image_defining_space;
matlabbatch{1}.spm.spatial.coreg.write.source = images_to_reslice;
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = interpolation_method;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = prefix;
spm_jobman('run',matlabbatch);

 
 
 
 
