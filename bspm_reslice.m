function bspm_reslice(reference, source, method, prefix)
% BSPM_RESLICE  Just a wrapper for re-slicing images
%
% USAGE: bspm_reslice(reference, source, method, prefix)
%
% Images should be defined in cell array of string(s). 
% Options for "method" are:
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
if nargin<3, method = 1; end
if nargin<2, error('USAGE: bspm_reslice(reference, source, method, prefix)'); end
if ischar(reference); reference = cellstr(reference); end
if ischar(source); source = cellstr(source); end

% run job         
matlabbatch{1}.spm.spatial.coreg.write.ref = reference;
matlabbatch{1}.spm.spatial.coreg.write.source = source;
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = method;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = prefix;

% run job
spm('defaults','fmri');
spm_jobman('run',matlabbatch);

 
 
 
 
