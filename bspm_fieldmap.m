function [] = bspm_fieldmap(magnitude_image, phase_image, epi_images, echo_times, total_epi_readout_time, blip, jacobTAG, method)
% BSPM_FIELDMAP
%
%   ARGUMENTS:
%       magnitude image = the magnitude image (if 2 entered, mean is taken)
%       phase image = the presubtracted phase image
%       epi_images = first epi of each run
%       echo_times = array of short and long echo times (e.g., [2.55 5.01])
%       total_epi_readout_time = echo spacing * # of lines of data acquired
%       blip = blip direction (-1 or 1) [DEFAULT = 1]
%       jacbobTAG = use jacobian modulation? (0 = NO, 1 = YES) [DEFAULT = 0]
%       method = unwrapping method ('Mark2D', 'Mark3D')
%

% ------------------------------------------------------- Copyright (C) 2014 -------------------------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<5
   disp('Give me arguments! bspm_fieldmap(magnitude_image, phase_image, epi_images, echo_times, total_epi_readout_time');
   return
end

% default blip and jacobian modulation
if nargin<7
    blip = -1;  
    jacobTAG = 0; 
end

% make sure image names are cell arrays of strings
if ischar(phase_image)
    phase_image = cellstr(phase_image);
end
if ischar(magnitude_image)
    magnitude_image = cellstr(magnitude_image);
end
if ischar(epi_images)
    epi_images = cellstr(epi_images);
end

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
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.et = echo_times;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.maskbrain = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.blipdir = blip;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.tert = total_epi_readout_time;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.epifm = 0;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.ajm = jacobTAG;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.method = method;   % Mark2D or Mark3D
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.fwhm = 10;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.pad = 0;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.ws = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.template = cellstr([spm_dir filesep 'templates/T1.nii']);
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.nerode = 2;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.reg = 0.02;
if length(epi_images)>1
    for r = 1:length(epi_images)
        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(r).epi = epi_images(r);
    end
else
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(1).epi = epi_images;
end
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchvdm = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.sessname = 'run';
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.writeunwarped = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.anat = '';
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchanat = 0;

% run job
spm('defaults','fmri'); spm_jobman('initcfg');         
spm_jobman('run',matlabbatch);

end

 
 
 
 
