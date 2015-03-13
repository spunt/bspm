function matlabbatch = bspm_fieldmap(mag, phase, epi, epi_rot, varargin)
% BSPM_FIELDMAP
%
%   USAGE: matlabbatch = bspm_fieldmap(mag, phase, epi, epi_rot, varargin)
%
%   ARGUMENTS:
%       mag     = the magnitude image (if 2 entered, mean is taken)
%       phase   = the presubtracted phase image
%       epi     = first epi of each run
%       epi_rot = echo spacing * # of lines of data acquired
%
%   VARARGIN: 
%       ets     = array of short and long echo times (e.g., [2.55 5.01])
%       (will attempt to determine this automatically from descrip field of hdr of mag images)
%       blip    = blip direction (-1 or 1) [DEFAULT = -1]
%       jacbob  = use jacobian modulation? (0 = NO, 1 = YES) [DEFAULT = 0]
%       method  = unwrapping method ('Mark2D', 'Mark3D')
%

% ------------------------------------------------------- Copyright (C) 2014 -------------------------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<4, disp('USAGE: matlabbatch = bspm_fieldmap(mag, phase, epi, epi_rot, varargin)'); return; end
def = { 'ets',          [], ...
        'blip',         -1, ...
        'jacob',        0, ...
        'method',       'Mark3D'};
bspm_setdefaults(def, varargin); 
if ischar(phase), phase = cellstr(phase); end
if ischar(mag), mag = cellstr(mag); end
if ischar(epi), epi = cellstr(epi); end
spm_dir = fileparts(which('spm'));
fm_dir  = fullfile(spm_dir, 'toolbox', 'FieldMap');
if length(mag)==2
    hdr     = spm_vol(char(mag));
    if isempty(ets)
        for i = 1:2
            tmp = hdr(i).descrip;
            [i1,i2] = regexp(tmp, 'TE=.+ms'); 
            ets(i) = str2double(tmp(i1+3:i2-2)); 
        end
    end
    imdata  = spm_read_vols(hdr);
    mag     = fullfile(fileparts(hdr(1).fname), 'mean_mag.nii'); 
    hdr(1).fname = mag;
    spm_write_vol(hdr(1), nanmean(imdata,4));
end
% build job variable
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.phase = cellstr(phase);
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.magnitude = cellstr(mag);
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.et = ets;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.maskbrain = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.blipdir = blip;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.tert = epi_rot;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.epifm = 0;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.ajm = jacob;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.method = method;   % Mark2D or Mark3D
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.fwhm = 10;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.pad = 0;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.ws = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.template = cellstr(fullfile(fm_dir, 'T1.nii'));
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.nerode = 2;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.reg = 0.02;
if length(epi)>1
    for r = 1:length(epi)
        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(r).epi = epi(r);
    end
else
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(1).epi = epi;
end
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchvdm = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.sessname = 'run';
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.writeunwarped = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.anat = '';
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchanat = 0;

% run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end


end
