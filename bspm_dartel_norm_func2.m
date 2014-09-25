function bspm_dartel_norm_func2(input)
% BSPM_DARTEL_NORM_FUNC
%
%   FIELDS:
%       epipat = cell array of patterns for grabbing epi images to norm
%       flowfields = flowfields (i.e. u_rc1*) (same length/order as images)
%       template = template image (i.e. Template_6.nii)
%       voxsize = voxel size for re-sampling (isotropic) 
%       fwhm = kernel for smoothing (isotropic) 
%

% --------------------------- Copyright (C) 2014 ---------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, error('No input!'); end
fn = {'epipat' 'flowfields' 'template' 'voxsize' 'fwhm'};
[status, msg] = checkfields(input, fn);
if ~status, error(msg); end
template = input.template;
flowfields = input.flowfields;
fwhm = input.fwhm;
voxsize = input.voxsize;
epipat = input.epipat;
if ischar(epipat), epipat = cellstr(epipat); end
images = cell(length(epipat),1);
for i = 1:length(epipat)
    images{i} = files(epipat{i});
end
if length(fwhm)==1, fwhm = repmat(fwhm,1,3); end
if length(voxsize)==1, voxsize = repmat(voxsize,1,3); end
if ischar(flowfields), flowfields = cellstr(flowfields); end
if ischar(template), template = cellstr(template); end

% number of subs
nsubs = length(images);

% Run DARTEL - Normalise
% --------------------------------------------------------
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
global defaults; 
defaults.smooth.prefix = sprintf('s%d',fwhm(1));
defaults.normalise.write.prefix = sprintf('w%d', voxsize(1)); 
spm_jobman('run',matlabbatch);

end

 
 
 
 
