%_______________________________________________________________________
%
% Initialize matlabbatch interface for Variability Toolbox
%_______________________________________________________________________
%
% This file is part of the Variability Toolbox for SPM
% Published by the Lifespan Neural Dynamics Group
% Provided under the GNU General Public License
% See LICENSE for terms and conditions

function tbx_def_Variability

matlabbatch{1}.spm.tools.variability = struct();

% Set modality
spm('defaults', 'FMRI')

% Avoid warning when launched directly from Matlab shell
spm_jobman('initcfg');

% Open the SPM Batch Editor with this module
spm_jobman('interactive', matlabbatch)
