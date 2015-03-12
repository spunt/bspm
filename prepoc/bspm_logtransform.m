function matlabbatch = bspm_logtransform(epi)
% BSPM_LOGTRANSFORM
%
%   USAGE: matlabbatch = bspm_logtransform(epi)
%

% -------------------------------- Copyright (C) 2014 --------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<1, disp('USAGE: matlabbatch = bspm_logtransform(epi)'); end
if ischar(epi), epi = cellstr(epi); end

% | build job variable
matlabbatch{1}.spm.tools.logtransform.data{1}           = strcat(epi, ',1');
matlabbatch{1}.spm.tools.logtransform.scaling.ascale    = 1;
matlabbatch{1}.spm.tools.logtransform.clipneg           = false;
matlabbatch{1}.spm.tools.logtransform.dtype             = 0; 

% | run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end

end
