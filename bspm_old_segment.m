function [] = bspm_old_segment(in)
% BSPM_OLD_SEGMENT
%
%   ARGUMENTS:
%       in = cell array containing paths to all images to segment
%

% ---------------------- Copyright (C) 2014 ----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1
   disp('USAGE: bspm_old_segment(in)');
   return
end

% make sure image names are cell arrays of strings
if ischar(in)
    in = cellstr(in);
end

% fix end of image filename cell array
in = cellstr([char(in) repmat(',1',length(in),1)]);

% path for tissue probability maps (in spm8/tpm) for use in segmentation 
spm_dir = which('spm');
spm_dir = regexprep(spm_dir,'spm.m','');
greymap=[spm_dir filesep '/tpm/grey.nii'];
whitemap=[spm_dir filesep '/tpm/white.nii'];
csfmap=[spm_dir filesep '/tpm/csf.nii'];

% build job
% -------------------------------------------------           
matlabbatch{1}.spm.spatial.preproc.data = in;
matlabbatch{1}.spm.spatial.preproc.output.GM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.WM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 0];
matlabbatch{1}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{1}.spm.spatial.preproc.output.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.opts.tpm{1} = greymap;
matlabbatch{1}.spm.spatial.preproc.opts.tpm{2} = whitemap;
matlabbatch{1}.spm.spatial.preproc.opts.tpm{3} = csfmap;
matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2 2 2 4];
matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;

% run job
spm('defaults','fmri'); spm_jobman('initcfg');         
spm_jobman('run',matlabbatch);

end

 
 
 
 
