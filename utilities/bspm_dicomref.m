function ref = bspm_dicomref(filepat)
% BSPM_DICOMREF
%
%   USAGE: ref = bspm_dicomref(filepat)
%
%         amy_2pt5
%         conte
%         lois
%         surf
% 
if nargin==0, disp('USAGE: ref = bspm_dicomref(filepat)'); end
bdir    = fileparts(whichdir('spm.m')); 
ref     = files(fullfile(bdir, 'dicom_refs', filepat));
end
 
