% Bramila pipeline demo for task-related group data

% NOTE: full 4D EPI data must fit into memory
% If you have multiple subjects, you should set 'group_cfg.save_memory = 1'
% Otherwise EPI data from all subjects are kept in memory

% group_cfg = group-level parameters, these are copied for all subjects
% cfg = subject-level parameters that will overwrite group-level parameters
% you should only modify subject-wise EPI+motion+tissue masks, all other
% parameters should be shared between all subjects

clc;
close all;
clear all;

% REQUIRED group parameters
group_cfg.bramilapath = pwd; % bramila path
addpath(group_cfg.bramilapath);
group_cfg.TR = 2.0;

% OPTIONAL group parameters (see "bramila_checkparameters_group.m" for defaults)
group_cfg.save_memory = 1; % save intermediate EPI data into disk instead of keeping all in the memory
group_cfg.do_temporal_tissue_ISC = 1; % remove shared tissue signal (1=only if significant with p<0.05, 2=always)
% group_cfg.do_spatial_ISC = 1; % remove synchronized voxels from nuisance masks (not recommended for unsmoothed data)
% group_cfg.ISC_mask = 'ISC_mask_fdr005.nii'; % if you already have ISC (or any other) mask
group_cfg.temporal_tissue_ISC_method = 'pca'; % use 'full' for no regularization (full, pls, pca)
group_cfg.do_spatial_smooth=1; % use 0 if already smoothed
group_cfg.smooth_FWHM = 6; % in millimeters
group_cfg.smooth_method = 'FSLgauss'; % 'SPM' (simple gaussian), 'FSL' (susan, non-gaussian), 'FSLgauss' (gaussian) 'AFNI' (iterative gaussian)
group_cfg.FSLDIR = 'fsl5.0-'; % path to your FSL binaries OR a proper command prefix
%group_cfg.save_path = '/triton/becs/scratch/braindata/kauttoj2/Deren_part2/TEST_GROUP2/group';

% OPTIONAL parameters shared between all subjects (see "bramila_checkparameters.m" for defaults)
% Note! These parameters are copied into all subject-wise csf files
group_cfg.white_mask_th = 0.90; % probability threshold for white matter tissue
group_cfg.csf_mask_th = 0.90;   % probability threshold for ventricles (csf tissue)
group_cfg.motion_reg_type = 'friston'; % motion regression type
% group_cfg.voxelsize=[2,2,2];        % voxelsize
% group_cfg.mask = 'MNI152_ENLARGED';  % initial EPI mask
% group_cfg.tissue_derivatives = 0;   % tissue regressor derivative order
group_cfg.min_tissue_var_expl = 50; % cut-off variance percentage for tissue PCA regressors (dynamic upper limit)
group_cfg.max_tissue_pca_count = 8; % maximum upper limit for tissue PCA regressors (despite variance explained)
% group_cfg.remove_global = 0;    % remove global signal
% group_cfg.mot_derivatives = 0;  % motion regressor derivatives
group_cfg.detrend_type='polynomial-nodemean';   % detrending type (linear or polynomial)
group_cfg.filter_type = 'butter';   % temporal filter type
group_cfg.filter_limits = [0,0.008,0.1,0.11];   % filter limits in Hz
group_cfg.write = 0;    % write all intermediate EPI's (otherwise only the final network)
% group_cfg.network_nodes = 'mynodes.mat'; % custom roi definition file that includes proper 'rois' variable

% REQUIRED subject parameters
%cfg{1}.infile
%cfg{1}.motionparam

% OPTIONAL subject parameters (otherwise use templates)

%warning! setting other parameters will override above group parameters!

k=1;
cfg{k}.infile = '/folder1/folder2/my_fmri_data/subj1.nii';
cfg{k}.motionparam = '/folder1/folder2/subj1_motion.txt';
%cfg{k}.fileID = 'Subject_1'; (OPTIONAL)
%cfg{k}.csf_mask = '/folder1/folder2/my_fmri_data/subj1_CSF_mask.nii';
%cfg{k}.wm_mask = '/folder1/folder2/my_fmri_data/subj1_WM_mask.nii';
cfg{k}.save_path =  '/folder1/folder2/my_fmri_data/results/';

k=2;
cfg{k}.infile = '/folder1/folder2/my_fmri_data/subj2.nii';
cfg{k}.motionparam = '/folder1/folder2/my_fmri_data/subj2_motion.txt';
%cfg{k}.fileID = 'Subject_2'; (OPTIONAL)
%cfg{k}.csf_mask = '/folder1/folder2/my_fmri_data/subj2_CSF_mask.nii';
%cfg{k}.wm_mask = '/folder1/folder2/my_fmri_data/subj2_WM_mask.nii';
cfg{k}.save_path =  '/folder1/folder2/my_fmri_data/results/';

%---- RUN BRAMILA CODE (do not edit, order of functions is fixed)

% (1) check parameters, mask, detrend and motion+tissue nuisance regression
[cfg,group_cfg] = bramila_clean_signal(cfg,group_cfg);
% (2) compute some diagnostics (optional)
[cfg,group_cfg] = bramila_diagnostics(cfg,group_cfg);
% (3) temporal filtering and spatial smoothing
cfg = bramila_spatiotemporal_filtering(cfg,group_cfg);
% (4) create network (optional)
[cfg,adj,tsdata]=bramila_makenet(cfg,group_cfg);

fprintf('\nAll Finished!\n')
