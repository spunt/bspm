function matlabbatch = bspm_dartel_norm_anat(images, flowfields, template, voxsize)
% BSPM_DARTEL_NORM_ANAT
%
%   ARGUMENTS:
%       images = anatomical files to norm (one per sub)
%       flowfields = flowfields (i.e. u_rc1*) (same length/order as images)
%       template = template image (i.e. Template_6.nii)
%       voxsize = voxel size for re-sampling (isotropic) [default = 1]
%       fwhm = kernel for smoothing (isotropic) [default = 0]
%

% ---------------------------- Copyright (C) 2014 ----------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3, disp('USAGE: bspm_dartel_norm_anat(images, flowfields, template, voxsize, fwhm)'); return; end
if nargin<4, voxsize = 1; end
fwhm = 0;
if length(voxsize)==1, voxsize = [voxsize voxsize voxsize]; end

% make sure image names are cell arrays of strings
if ischar(images)
    images = cellstr(images);
end
if ischar(flowfields)
    flowfields = cellstr(flowfields);
end
if ischar(template)
    template = cellstr(template);
end

% add ',1' to images
for i = 1:length(images)
    images(i) = cellstr([images{i} ',1']);
end

% number of subs
nsubs = length(images);

% Run DARTEL - Normalise
% -------------------------------------------------
matlabbatch{1}.spm.tools.dartel.mni_norm.template = cellstr(template);
for s = 1:nsubs
    matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj(s).flowfield = flowfields(s);
    matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj(s).images = images(s);
end
matlabbatch{1}.spm.tools.dartel.mni_norm.vox = voxsize;
matlabbatch{1}.spm.tools.dartel.mni_norm.bb = [-78 -112 -50; 78 76 85];
matlabbatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
matlabbatch{1}.spm.tools.dartel.mni_norm.fwhm = [fwhm fwhm fwhm];

% run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end


end
