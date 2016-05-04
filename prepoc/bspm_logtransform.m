function matlabbatch = bspm_logtransform(epi, is4D)
% BSPM_LOGTRANSFORM
%
%   USAGE: matlabbatch = bspm_logtransform(epi, is4D)
%

% -------------------------------- Copyright (C) 2014 --------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<1, mfile_showhelp; return; end
if nargin<2, is4D = 0; end; 
if ischar(epi), epi = cellstr(epi); end
if is4D, epi = bspm_expand4D(char(epi)); 
else epi = strcat(epi, ',1'); end

% | build job variable
matlabbatch{1}.spm.tools.logtransform.data{1}           = epi; 
matlabbatch{1}.spm.tools.logtransform.scaling.ascale    = 1;
matlabbatch{1}.spm.tools.logtransform.clipneg           = false;
matlabbatch{1}.spm.tools.logtransform.dtype             = 0; 

% | run job
if nargout==0, bspm_runbatch(matlabbatch); end

end
