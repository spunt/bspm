function matlabbatch = bspm_clinical_segment_normalize(t1, lesion, varargin)
% BSPM_CLINICAL_SEGMENT_NORMALIZE
%
%   USAGE: matlabbatch = bspm_clinical_segment_normalize(t1, lesion, varargin)
% ____________________________________________________________________________
% REQUIRED ARGUMENTS
%   t1:       T1-weighted image filename(s)
%   lesion:   lesion map(s) [registered with t1 image(s)]
% ____________________________________________________________________________
% OPTIONAL ARGUMENTS (name-value pairs)
%   template:           template to normalize to ('older' or 'younger')
%   cleanuplevel:       tissue cleanup (0=none, 1=light, 2=thorough)
%   voxelsize:          resolution for resampling normalized images
%   boundingbox:        bounding box for normalized images
%   skullstripthresh:   prob threshold for skull stripping
%   delintermediate:    delete intermediate images (0 or 1)
%   enantiomorphic:     enantiomorphic normalization (see below)
%   autosetorigin:      attempt to auto-orient to ac/pc origin
% ____________________________________________________________________________
% If you choose 'Enantiomorphic normalization' the lesion will be replaced
% by tissue from the intact hemisphere, if you set this to false a
% traditional lesion-masked normalization will occur. In general
% enantiomorphic is better for large unilateral lesions, lesion masking is
% required for bilateral lesions and these methods perform similarly for
% smaller lesions (see Nachev et al., 2008).
% 

% ---------------------------------- Copyright (C) 2015 ----------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: May_26_2015
def = { 'template',             'older'                     , ...
        'cleanuplevel',         2                           , ...
        'voxelsize',            [1 1 1]                     , ...
        'boundingbox',          [-78 -112 -50; 78 76 85]    , ...
        'skullstripthresh'      0.005                       , ...                                                                                         
        'delintermediate',      1                           , ...
        'enantiomorphic',       1                           , ...
        'autosetorigin',        1                           , ...
        };
vals = setargs(def, varargin);
if nargin < 2, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals); 

% | CHECK INPUTS
if ischar(t1), t1 = cellstr(t1); end
if ischar(lesion), lesion = cellstr(lesion); end
templatetag = 1; 
if strcmpi(template, 'younger'), templatetag = 0; end

% | BUILD JOB
matlabbatch{1}.spm.tools.MRI.MRnormseg.anat             = strcat(t1, ',1');
matlabbatch{1}.spm.tools.MRI.MRnormseg.les              = strcat(lesion, ',1');
matlabbatch{1}.spm.tools.MRI.MRnormseg.t2               = '';
matlabbatch{1}.spm.tools.MRI.MRnormseg.clinicaltemplate = templatetag;
matlabbatch{1}.spm.tools.MRI.MRnormseg.clean            = cleanuplevel;
matlabbatch{1}.spm.tools.MRI.MRnormseg.bb               = boundingbox;
matlabbatch{1}.spm.tools.MRI.MRnormseg.vox              = voxelsize;
matlabbatch{1}.spm.tools.MRI.MRnormseg.ssthresh         = skullstripthresh;
matlabbatch{1}.spm.tools.MRI.MRnormseg.DelIntermediate  = delintermediate;
matlabbatch{1}.spm.tools.MRI.MRnormseg.Enantiomorphic   = enantiomorphic;
matlabbatch{1}.spm.tools.MRI.MRnormseg.AutoSetOrigin    = autosetorigin;

% | RUN IF NO OUTPUT ARGS SPECIFIED
if nargout==0, bspm_runbatch(matlabbatch); end

end




