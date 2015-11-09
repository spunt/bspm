function masking = tbx_cfg_masking
% MATLABBATCH Configuration file for toolbox 'Mask Creation'
% This code has been automatically generated.
addpath(fullfile(spm('Dir'), 'toolbox', 'Masking'))
% ---------------------------------------------------------------------
% innames Input Images
% ---------------------------------------------------------------------
innames         = cfg_files;
innames.tag     = 'innames';
innames.name    = 'Input Images';
innames.filter = 'image';
innames.ufilter = '.*';
innames.num     = [1 Inf];
% ---------------------------------------------------------------------
% avgexpr Average Expression
% ---------------------------------------------------------------------
avgexpr         = cfg_menu;
avgexpr.tag     = 'avgexpr';
avgexpr.name    = 'Average Expression';
avgexpr.help    = {'Type of average to compute'};
avgexpr.def = @(val)tbx_def_masking('makeavg.avgexpr',val{:});
avgexpr.labels = {
                  'Arithmetic mean'
                  'Median'
                  }';
avgexpr.values = {
                  'mean(X)'
                  'median(X)'
                  }';
% ---------------------------------------------------------------------
% outname Output Filename
% ---------------------------------------------------------------------
outname         = cfg_entry;
outname.tag     = 'outname';
outname.name    = 'Output Filename';
outname.help    = {
                   'The output image is saved with this name. '
                   'If a path name is given here, the output directory setting will be ignored. If there is no path and no output directory, then the current working directory is used'
                   }';
outname.def = @(val)tbx_def_masking('makeavg.outname',val{:});
outname.strtype = 's';
outname.num     = [1  Inf];
% ---------------------------------------------------------------------
% outdir Output Directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output Directory';
outdir.help    = {'Output files will be written into this directory. If no directory is given, images will be written to current working directory. If both output filename and output directory contain a directory, then output filename takes precedence.'};
outdir.def = @(val)tbx_def_masking('makeavg.outdir',val{:});
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [0 1];
% ---------------------------------------------------------------------
% makeavg Make Average
% ---------------------------------------------------------------------
makeavg         = cfg_exbranch;
makeavg.tag     = 'makeavg';
makeavg.name    = 'Make Average';
makeavg.val     = {innames avgexpr outname outdir };
makeavg.help    = {
                   'The Masking toolbox contains code for creating a thresholded binary image that can be specified as an explicit analysis mask when setting up a statistical model in SPM. '
                   'This user-interface implements a strategy for automatically computing a threshold for an image which is typically an average of the scans to be analysed.'
                   'This module first creates an average from the selected images.'
                   }';
makeavg.prog = @make_average;
makeavg.vout = @make_average_vout;
% ---------------------------------------------------------------------
% inname Input Image
% ---------------------------------------------------------------------
inname         = cfg_files;
inname.tag     = 'inname';
inname.name    = 'Input Image';
inname.filter = 'image';
inname.ufilter = '.*';
inname.num     = [1 1];
% ---------------------------------------------------------------------
% optfunc Optimality Criterion
% ---------------------------------------------------------------------
optfunc         = cfg_menu;
optfunc.tag     = 'optfunc';
optfunc.name    = 'Optimality Criterion';
optfunc.help    = {'Choice of objective. See help on Optimal Thresholding.'};
optfunc.def = @(val)tbx_def_masking('optthr.optfunc',val{:});
optfunc.labels = {
                  'Corr(img, img > thr)'
                  'Luo-Nichols anti-mode'
                  }';
optfunc.values = {
                  '@opt_thr_corr'
                  '@opt_thr_antimode'
                  }';
% ---------------------------------------------------------------------
% outname Output Filename
% ---------------------------------------------------------------------
outname         = cfg_entry;
outname.tag     = 'outname';
outname.name    = 'Output Filename';
outname.help    = {
                   'The output image is saved with this name. '
                   'If a path name is given here, the output directory setting will be ignored. If there is no path and no output directory, then the current working directory is used'
                   }';
outname.def = @(val)tbx_def_masking('optthr.outname',val{:});
outname.strtype = 's';
outname.num     = [1  Inf];
% ---------------------------------------------------------------------
% outdir Output Directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output Directory';
outdir.help    = {'Output files will be written into this directory. If no directory is given, images will be written to current working directory. If both output filename and output directory contain a directory, then output filename takes precedence.'};
outdir.def = @(val)tbx_def_masking('optthr.outdir',val{:});
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [0 1];
% ---------------------------------------------------------------------
% optthr Optimal Thresholding
% ---------------------------------------------------------------------
optthr         = cfg_exbranch;
optthr.tag     = 'optthr';
optthr.name    = 'Optimal Thresholding';
optthr.val     = {inname optfunc outname outdir };
%%
optthr.help    = {
                  'The Masking toolbox contains code for creating a thresholded binary image that can be specified as an explicit analysis mask when setting up a statistical model in SPM. '
                  'This user-interface implements a strategy for automatically computing a threshold for an image which is typically an average of the scans to be analysed.'
                  'This module optimally thresholds an image, according to one of two objective functions.'
                  'The first is the one from the paper by Ridgway et al, '
                  ' http://dx.doi.org/10.1016/j.neuroimage.2008.08.045 '
                  'It is based on maximising the correlation between the original and thresholded images. This is equivalent to finding the threshold that gives the largest two-sample t-statistic between the above and below threshold parts of the data, which is similar to the classic Otsu thresholding criterion, based on maximum between-class separation with the histrogram. '
                  ' http://en.wikipedia.org/wiki/Otsu''s_method'
                  'Please cite Ridgway et al. if you find the toolbox useful.'
                  ''
                  'The second is described in appendix B of the paper '
                  ' Diagnosis and exploration of massively univariate neuroimaging models'
                  ' Luo & Nichols (2003) Neuroimage 19:1014-1032 '
                  ' http://www.sph.umich.edu/ni-stat/SPMd/SPMd.pdf '
                  ' http://dx.doi.org/10.1016/S1053-8119(03)00149-6 '
                  'Please cite both Ridgway et al. and Luo & Nichols if you use this method.'
                  }';
%%
optthr.prog = @opt_thresh;
optthr.vout = @opt_thresh_vout;
% ---------------------------------------------------------------------
% masking Mask Creation
% ---------------------------------------------------------------------
masking         = cfg_repeat;
masking.tag     = 'masking';
masking.name    = 'Mask Creation';
%%
masking.help    = {
                   'The Masking toolbox contains code for creating a thresholded binary image that can be specified as an explicit analysis mask when setting up a statistical model in SPM. '
                   'This user-interface implements a strategy for automatically computing a threshold for an image which is typically an average of the scans to be analysed.'
                   'Command-line functions are available for some other approaches,for example make_majority_mask.m. '
                   'These strategies, and related aspects, are described in the paper'
                   ' Issues with threshold masking in Voxel Based Morphometry of atrophied brains. '
                   ' Ridgway, G.R.; Omar, R.; Ourselin, S.; Hill, D.L.G.; Warren, J.D.; and Fox, N.C. (in press) NeuroImage'
                   ' http://dx.doi.org/10.1016/j.neuroimage.2008.08.045 '
                   ' pre-print available from http://eprints.ucl.ac.uk/13060/'
                   'Please cite this paper if you find the toolbox useful.'
                   }';
%%
masking.values  = {makeavg optthr };
masking.num     = [0 Inf];
