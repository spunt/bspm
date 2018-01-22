function matlabbatch = bspm_dartel_existing_template(rc1, rc2, templatedir)
% BSPM_DARTEL_CREATE_TEMPLATE
%
%   USAGE: matlabbatch = bspm_dartel_existing_template(rc1, rc2, templatedir)
%   ARGUMENTS:
%

% ------------------------ Copyright (C) 2014 ------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<1, mfile_showhelp; return; end
if iscell(templatedir), templatedir = char(templatedir); end
if ischar(rc1), rc1 = cellstr(rc1); end
if ischar(rc2), rc2 = cellstr(rc2); end
if length(rc1)~=length(rc2), error('Number of rc1 and rc2 images must be equal!'); end

% | Parameters
rparam = [  4 2 1e-06; ...
            2 1 1e-06; ...
            1 0.5 1e-06; ...
            0.5 0.25 1e-06; ...
            0.25 0.125 1e-06; ...
            0.25 0.125 1e-06 ...
            ]; 
K   = [0 0 1 2 4 6];
its = [3 3 3 3 3 3];

% Build job for DARTEL (Create Template)
% -------------------------------------------------
matlabbatch{1}.spm.tools.dartel.warp1.images{1} = strcat(rc1, ',1');
matlabbatch{1}.spm.tools.dartel.warp1.images{2} = strcat(rc2, ',1');
matlabbatch{1}.spm.tools.dartel.warp1.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.its = 3;
for i = 1:6
    matlabbatch{1}.spm.tools.dartel.warp1.settings.param(i).its = its(i);
    matlabbatch{1}.spm.tools.dartel.warp1.settings.param(i).rparam = rparam(i,:); 
    matlabbatch{1}.spm.tools.dartel.warp1.settings.param(i).K = K(i);
    matlabbatch{1}.spm.tools.dartel.warp1.settings.param(i).template = {fullfile(templatedir, sprintf('Template_%d.nii', i))};
end

% run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end


end
