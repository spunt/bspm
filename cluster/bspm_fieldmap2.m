function bspm_fieldmap2(input)
% BSPM_FIELDMAP2
%
%   input structure must contain the following fields
%
%       .epipat - first epi of each run
%       .magimg - the magnitude image (if 2 entered, mean is taken)
%       .phaseimg - the presubtracted phase image
%       .echotimes - array of short and long echo times (e.g., [2.55 5.01])
%       .epirot - echo spacing * # of lines of data acquired
%       .blip - blip direction (-1 or 1) [DEFAULT = -1]
%       .jacob - use jacobian modulation? (0 = NO, 1 = YES) [DEFAULT = 0]
%       .method - unwrapping method ('Mark2D', 'Mark3D') [DEFAULT = 'Mark3D']
%

% ---------------------------- Copyright (C) 2014 ----------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, error('No input!'); end
if ~isfield(input, 'blip'), input.blip = -1; end
if ~isfield(input, 'method'), input.method = 'Mark3D'; end
if ~isfield(input, 'jabob'), input.jacob = 0; end
fn = {'epipat' 'magimg' 'phaseimg' 'echotimes' 'epirot' 'blip' 'jacob' 'method'};
[status, msg] = checkfields(input, fn);
if ~status, error(msg); end
phase_image = input.phaseimg;
magnitude_image = input.magimg;
epiimg = files(input.epipat);
epiimg = epiimg(1);

% grab spm directory
spm_dir = which('spm');
spm_dir = regexprep(spm_dir,'spm.m','');

% if two magnitude images, take average
if length(magnitude_image)==2
    hdr = spm_vol(char(magnitude_image));
    imdata = spm_read_vols(hdr);
    immean = mean(imdata,4);
    magnitude_image = [fileparts(hdr(1).fname) filesep 'mean_mag.img'];
    hdr(1).fname = magnitude_image;
    spm_write_vol(hdr(1), immean);
end
    
% build job variable
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.phase = cellstr(phase_image);
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.magnitude = cellstr(magnitude_image);
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.et = input.echotimes;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.maskbrain = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.blipdir = input.blip;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.tert = input.epirot;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.epifm = 0;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.ajm = input.jacob;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.method = input.method;   % Mark2D or Mark3D
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.fwhm = 10;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.pad = 0;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.ws = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.template = cellstr([spm_dir filesep 'templates/T1.nii']);
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.nerode = 2;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.reg = 0.02;
if length(epiimg) > 1
    for r = 1:length(epiimg)
        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(r).epi = epiimg(r);
    end
else
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(1).epi = epiimg;
end
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchvdm = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.sessname = 'run';
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.writeunwarped = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.anat = '';
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchanat = 0;

% run job
spm('defaults','fmri');      
spm_jobman('run',matlabbatch);

end

 
 
 
 
