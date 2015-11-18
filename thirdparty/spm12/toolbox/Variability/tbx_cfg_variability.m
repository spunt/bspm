%_______________________________________________________________________
%
% Generate GUI for SPM's matlabbatch interface
%_______________________________________________________________________
%
% Output
%
% variability | matlabbatch interface (struct)
%
%_______________________________________________________________________
%
% This file is part of the Variability Toolbox for SPM
% Published by the Lifespan Neural Dynamics Group
% Provided under the GNU General Public License
% See LICENSE for terms and conditions

function variability = tbx_cfg_variability

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','Variability')); end

%--------------------------------------------------------------------------
% modeltype Model Type
%--------------------------------------------------------------------------

modeltype        = cfg_menu;
modeltype.tag    = 'modeltype';
modeltype.name   = 'Model type';
modeltype.help   = {
'Select a model type to compute variabilty'
''
'When choosing ''Boxcar'' the condition-wise variability is computed across all blocks of all sessions.'
''
};
modeltype.labels = {'Boxcar'}';
modeltype.values = {'boxcar'};
modeltype.val    = {'boxcar'};

% ---------------------------------------------------------------------
% modelmat Model File
% ---------------------------------------------------------------------
modelmat         = cfg_files;
modelmat.tag     = 'modelmat';
modelmat.name    = 'Model file';
modelmat.help    = {
'Select the MAT-file with the first-level model specification.'
''
'This file can be generated with the SPM batch editor by setting up a model for first-level analysis and saving the specification by selecting File > Save Batch from the menu bar of the batch editor. When processing multiple subjects this file could be used as a template for scripting.'
''
'The Boxcar model uses the following information of the first-level specification:'
''
'* Timing parameters > Units for design'
'* Timing parameters > Interscan Interval'
'* Data & Design > Subject/Session > Scans'
'* Data & Design > Subject/Session > Conditions > Condition > Name'
'* Data & Design > Subject/Session > Conditions > Condition > Onsets'
'* Data & Design > Subject/Session > Conditions > Condition > Durations'
'* Data & Design > Subject/Session > Regressors > Regressor > Name'
'* Data & Design > Subject/Session > Regressors > Regressor > Value'
'* Data & Design > Explicit Mask'
}';

modelmat.filter  = 'mat';
modelmat.ufilter = '.*\.mat$';
modelmat.num     = [1 1];


%--------------------------------------------------------------------------
% metric Varibility Metric
%--------------------------------------------------------------------------

metric        = cfg_menu;
metric.tag    = 'metric';
metric.name   = 'Variability function';
metric.help   = {
''
'VAR: variance'
'SD: standard deviation'
'MSSD: mean square successive difference'
'SQRT(MSSD): square root of MSSD'
''
'For boxcar models VAR and SD will be detrended prior to the variability computation.'
''
};
metric.labels = {'VAR' 'SD' 'MSSD' 'SQRT(MSSD)'}';
metric.values = {'var' 'sd' 'mssd' 'sqrt_mssd'};
metric.val 		= {'sd'};

% ---------------------------------------------------------------------
% resultprefix Result prefix
% ---------------------------------------------------------------------

resultprefix         = cfg_entry;
resultprefix.tag     = 'resultprefix';
resultprefix.name    = 'Filename prefix';
resultprefix.help    = {
'By default the result files are named after the condition and prefixed with the string entered here (e.g. sd_1-back.nii). To avoid overwriting result files when processing multiple subjects or metrics adapt this prefix or choose a different result directory for each subject and metric. '
};
resultprefix.strtype = 's';
resultprefix.val     = {'sd'};
resultprefix.num     = [0 Inf];

% ---------------------------------------------------------------------
% resultdir Result directory
% ---------------------------------------------------------------------

resultdir         = cfg_files;
resultdir.tag     = 'resultdir';
resultdir.name    = 'Result directory';
resultdir.help    = {
'Select the directory where the results should be saved.'
};
resultdir.filter  = 'dir';
resultdir.ufilter = '.*';
resultdir.num     = [1 1];

%--------------------------------------------------------------------------
% Variability Calculation of variability in series of volumes
%--------------------------------------------------------------------------
variability				= cfg_exbranch;
variability.tag		= 'variability';
variability.name	= 'Variability';
variability.val		= {modeltype modelmat metric resultprefix resultdir};
variability.help	= {
''
'Variability Toolbox for SPM'
'Version 0.1'
''
'This toolbox measures voxel-based temporal variability in fMRI data, and will output Nifti files. The toolbox structure is intended to be similar to a standard SPM first-level analysis. You can then proceed with those level-1 variability-based images to a level-2 analysis in order to model effects of interest. However, the nii output files could also be used within other statistics programs of your choice (e.g., FSL, AFNI, PLS).'
''
'The toolbox currently supports modeling block designs with a boxcar model, and computes variability using measures such as: detrended standard deviation (SD), variance (VAR), mean squared successive difference (MSSD) and SQRT(MSSD). We plan to add additional modeling approaches and variability measures in future releases.'
''
'To use the toolbox, you would specify your model in the first-level analysis as usual. You then save the model, which is used as input for the toolbox to specify all sessions, conditions, onsets, durations, and nuisance regressors.'
''
'For further information visit'
'http://douglasdgarrett.com/software'
''
'To contribute, check out the current release at'
'http://github.com/stefanschmidt/VariabilityToolbox'
''
'Published by:'
'Lifespan Neural Dynamics Group (PI: Douglas Garrett; Developer: Stefan Schmidt),'
'Max Planck UCL Centre for Computational Psychiatry and Ageing Research,'
'Max Planck Institute for Human Development, Berlin, Germany.'
''
'Provided under the GNU General Public License, see LICENSE for terms and conditions.'
''
'Caveats'
''
'- data transfer may be slow when accessing SMB network shares with OS X'
''
'Release History'
''
'0.1 (2014-12-20)'
'- initial release'
};
variability.prog	= @tbx_run_variability;
variability.vout	= @vout;

%--------------------------------------------------------------------------
% validation of input data
%--------------------------------------------------------------------------

function dep = vout(varargin)
	dep(1)            = cfg_dep;
	dep(1).sname      = 'Condition-wise variability of a series of volumes';
	dep(1).src_output = substruct('.','variability');
	dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

end
