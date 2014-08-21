function bspm_dartel_norm_func(images, flowfields, template, voxsize, fwhm)
% BSPM_DARTEL_NORM_FUNC
%
%   ARGUMENTS:
%       images = array of cell arrays containing epi files to norm (length = nsubs)
%       flowfields = flowfields (i.e. u_rc1*) (same length/order as images)
%       template = template image (i.e. Template_6.nii)
%       voxsize = voxel size for re-sampling (isotropic) [default = 3]
%       fwhm = kernel for smoothing (isotropic) [default = 6]
%

% ------------------------------- Copyright (C) 2014 -------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3, error('USAGE: bspm_dartel_norm_func(images, flowfields, template, voxsize, fwhm)'); end
if nargin<4, voxsize = 3; end
if nargin<5, fwhm = 8; end
if length(fwhm)==1, fwhm = repmat(fwhm,1,3); end
if length(voxsize)==1, voxsize = repmat(voxsize,1,3); end
if ischar(flowfields), flowfields = cellstr(flowfields); end
if ischar(template), template = cellstr(template); end
if ischar(images), images = cellstr(images); end
if ischar(images{1}), images = {images}; end
nsubs = length(images);

% Run DARTEL - Normalise
% -------------------------------------------------
matlabbatch{1}.spm.tools.dartel.mni_norm.template = cellstr(template);
for s = 1:nsubs
    cim = images{s};
    for i = 1:length(cim)
        cim(i) = cellstr([cim{i} ',1']);
    end
    matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj(s).flowfield = flowfields(s);
    matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj(s).images = cim;
end                     
matlabbatch{1}.spm.tools.dartel.mni_norm.vox = voxsize;
matlabbatch{1}.spm.tools.dartel.mni_norm.bb = [-78 -112 -50; 78 76 85];
matlabbatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
matlabbatch{1}.spm.tools.dartel.mni_norm.fwhm = fwhm;
spm_jobman('initcfg');    
global defaults
spm_get_defaults('normalise.write.prefix',sprintf('w%d',voxsize*100));
spm_get_defaults('smooth.prefix',sprintf('s%d',fwhm));
spm_jobman('run',matlabbatch);

end

 
 
 
 
