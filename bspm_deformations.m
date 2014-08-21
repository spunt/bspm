function [] = bspm_deformations(in, field)
% BSPM_DEFORMATIONS
%
%   USAGE: bspm_deformations(in, field)
%
%   ARGUMENTS:
%       in = cell array of images to apply deformation field to
%       field = deformation field to use
%

% --------------------- Copyright (C) 2014 ---------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<2, disp('USAGE: bspm_deformations(in, field)'); return; end

% make sure image names are cell arrays of strings
if ischar(in), in = cellstr(in); end
if iscell(field), field = char(field); end

% fix end of image filename cell array
for i = 1:length(in), in(i) = cellstr([in{i} ',1']); end

% build job
% -------------------------------------------------         
matlabbatch{1}.spm.util.defs.comp{1}.def = cellstr(field);
matlabbatch{1}.spm.util.defs.ofname = '';
matlabbatch{1}.spm.util.defs.fnames = in;
matlabbatch{1}.spm.util.defs.savedir.savesrc = 1;
matlabbatch{1}.spm.util.defs.interp = 4;


% run job
spm('defaults','fmri');   
spm_jobman('run',matlabbatch);

end

 
 
 
 
