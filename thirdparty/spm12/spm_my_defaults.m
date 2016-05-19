function spm_my_defaults
% SPM_MY_DEFAULTS
global defaults

% Mask defaults
%==========================================================================
defaults.mask.thresh    = 0.8; % default is 0.8
% change defaults.mask.thresh to -inf to turn off implicit masking (this is
% very useful if you have one or two subjects who have a lot of signal
% loss, since you can use an explicit mask and all voxels in the explicit
% mask will be included, including those with signal loss)

% Stats defaults
%==========================================================================
defaults.stats.maxmem      = 2^33;  % maximum amount of RAM to use (2^33 = 8GB)
defaults.stats.maxres      = 128;
defaults.stats.resmem      = true;  % flag to store temporary files on disk (false) or in memory (true)
defaults.stats.fmri.ufp    = 0.001;  % Upper tail F-probability
% change defaults.stats.fmri.ufp to relax the threshold for defining voxels
% for variance estimation at the ReML stage. By default, SPM does the
% variance calculations on only those voxels that survive an F-test at
% .001. If there are no voxels which survive this threshold (which is
% commonly the case when estimating models within small volumes), the
% estimation will fail. In this case, you can relax the threshold for
% finding the voxels to do the estimation on.