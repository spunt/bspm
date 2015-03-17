function matlabbatch = bspm_dartel_norm_anat(images, flowfields, template, voxsize)
% BSPM_DARTEL_NORM_ANAT
%
%   ARGUMENTS:
%       images = anatomical files to norm (one per sub)
%       flowfields = flowfields (i.e. u_rc1*) (same length/order as images)
%       template = template image (i.e. Template_6.nii)
%       voxsize = voxel size for re-sampling (isotropic) [default = 1]
%

% ---------------------------- Copyright (C) 2014 ----------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3, disp('USAGE: bspm_dartel_norm_anat(images, flowfields, template, voxsize, fwhm)'); return; end
if nargin<4, voxsize = 1; end
fwhm = zeros(1,3); 
if length(voxsize)==1, voxsize = repmat(voxsize,1,3); end

% make sure image names are cell arrays of strings
if ischar(images), images = cellstr(images); end
if ischar(flowfields), flowfields = cellstr(flowfields); end
if ischar(template), template = cellstr(template); end
nsubs = length(images);

% Run DARTEL - Normalise
% -------------------------------------------------
matlabbatch{1}.spm.tools.dartel.mni_norm.template = template;
matlabbatch{1}.spm.tools.dartel.mni_norm.data.subjs.flowfields = flowfields; 
matlabbatch{1}.spm.tools.dartel.mni_norm.data.subjs.images{1} = images; 
matlabbatch{1}.spm.tools.dartel.mni_norm.vox = voxsize;
matlabbatch{1}.spm.tools.dartel.mni_norm.bb = [-78 -112 -50; 78 76 85];
matlabbatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
matlabbatch{1}.spm.tools.dartel.mni_norm.fwhm = fwhm; 

% run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end

end
