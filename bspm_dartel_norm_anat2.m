function bspm_dartel_norm_anat2(input)
% BSPM_DARTEL_NORM_ANAT
%
%   ARGUMENTS:
%       images = anatomical files to norm (one per sub)
%       flowfields = flowfields (i.e. u_rc1*) (same length/order as images)
%       template = template image (i.e. Template_6.nii)
%

% --------------------------- Copyright (C) 2014 ---------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, error('No input!'); end
fn = {'images' 'flowfields' 'template'};
[status, msg] = checkfields(input, fn);
if ~status, error(msg); end
voxsize = [1];
fwhm = [0];

images = input.images;
flowfields = input.flowfields;
template = input.template;

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
spm('defaults','fmri'); spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

end

 
 
 
 
