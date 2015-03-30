function ref = bspm_dicomref(filepat)
% BSPM_DICOMREF
%
%   USAGE: ref = bspm_dicomref(filepat)
%
pdir    = pwd; 
[a,b]   = regexpi(pdir, 'Matlab');
pdir    = pdir(1:b); 
ref     = files(fullfile(pdir, 'dicom_refs', filepat));
end
 
