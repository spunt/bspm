function matlabbatch = bspm_smooth(in, fwhm, implicitmaskTAG)
% BSPM_SMOOTH
%
%   USAGE: matlabbatch = bspm_smooth(in, fwhm, prefix)
%
%       in              =  array of images to smooth (full path)
%       fwhm            =  smoothing kernel
%       implicitmaskTAG = 0: no implicit mask, 1: implicit mask
%

% --------------- Copyright (C) 2014 ---------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<2, disp('USAGE: matlabbatch = bspm_smooth(in, fwhm, implicitmaskTAG)'); return, end
if nargin<3, implicitmaskTAG = 0; end
if length(fwhm)==1, fwhm = [fwhm fwhm fwhm]; end
if ischar(in), in = cellstr(in); end
in = strcat(in, ',1');

% | Build job variable
% | =======================================================================
matlabbatch{1}.spm.spatial.smooth.data  = cellstr(in);
matlabbatch{1}.spm.spatial.smooth.fwhm  = fwhm;
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im    = implicitmaskTAG;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';

% | Run job (only if no output arguments requested)
% | =======================================================================
if nargout==0, spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end

end
