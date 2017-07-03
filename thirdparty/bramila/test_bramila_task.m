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
group_cfg.bramilapath = '/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/bramila'; % bramila path
group_cfg.TR = 2.015;

% OPTIONAL group parameters (see "bramila_checkparameters_group.m" for defaults)
% group_cfg.save_memory = 1; % save intermediate EPI data into disk instead of keeping all in the memory
% group_cfg.do_temporal_tissue_ISC = 1; % remove shared tissue signal
% group_cfg.do_spatial_ISC = 1; % remove synchronized voxels from nuisance tissue masks
% group_cfg.ISC_mask = 'ISC_mask_fdr005.nii'; % if you already have ISC (or any other) mask
% group_cfg.temporal_tissue_ISC_method = 'pca'; % use 'full' for no regularization (full, pls, pca)
% group_cfg.save_path = '/triton/becs/scratch/braindata/kauttoj2/code/bramila/testdata';

% OPTIONAL parameters shared between all subjets (see "bramila_checkparameters.m" for defaults)
% Note! These parameters are copied into all subject-wise csf files
 group_cfg.white_mask_th = 0.80; % probability threshold for white matter tissue
 group_cfg.csf_mask_th = 0.80;   % probability threshold for ventricle (csf tissue)
 group_cfg.motion_reg_type = 'friston'; % motion regression type
% group_cfg.voxelsize=[2,2,2];        % voxelsize
% group_cfg.mask = 'MNI152_ENLARGED';  % initial EPI mask
% group_cfg.tissue_derivatives = 0;   % tissue regressor derivative order
 group_cfg.min_tissue_var_expl = 50; % minimum variance percentage for tissue PCA regressors
 group_cfg.max_tissue_pca_count = 15; % upper limit for tissue PCA regressors
% group_cfg.remove_global = 0;    % remove global signal
% group_cfg.mot_derivatives = 0;  % motion regressor derivatives
% group_cfg.detrend_type='linear-nodemean';   % detrending type
% group_cfg.filter_type = 'butter';   % temporal filter type
% group_cfg.filter_limits = [0,0.01,0.08,0.09];   % filter limits in Hz
% group_cfg.write = 0;    % write all intermediate EPI's (otherwise only the final network)
% group_cfg.network_nodes = 'mynodes.mat'; % custom roi definition file that includes proper 'rois' variable

% REQUIRED subject parameters
%cfg{1}.infile
%cfg{1}.motionparam

% OPTIONAL subject parameters (otherwise usetemplates)
%cfg{1}.csf_mask
%cfg{1}.wm_mask
%warning! setting other parameters will override above group parameters!

k=1;
cfg{k}.infile = '/scratch/braindata/kauttoj2/Tulitikku/new_fdpa/FunImgNormDetrendedArtrepSmoothedNii/S3/4D.nii';
cfg{k}.motionparam = '/scratch/braindata/kauttoj2/Tulitikku/new_fdpa/RealignParameter/S3/rp_epi_0003.txt';
cfg{k}.wm_mask = '/scratch/braindata/kauttoj2/Tulitikku/new_fdpa/T1ImgSegment/S3/wc2coSPGR.img';
cfg{k}.csf_mask = '/scratch/braindata/kauttoj2/Tulitikku/new_fdpa/T1ImgSegment/S3/wc3coSPGR.img';
cfg{k}.save_path = '/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/bramila/testdata/S3';

k=2;
cfg{k}.infile = '/scratch/braindata/kauttoj2/Tulitikku/new_fdpa/FunImgNormDetrendedArtrepSmoothedNii/S5/4D.nii';
cfg{k}.motionparam = '/scratch/braindata/kauttoj2/Tulitikku/new_fdpa/RealignParameter/S5/rp_epi_0003.txt';
cfg{k}.wm_mask = '/scratch/braindata/kauttoj2/Tulitikku/new_fdpa/T1ImgSegment/S5/wc2coSPGR.img';
cfg{k}.csf_mask = '/scratch/braindata/kauttoj2/Tulitikku/new_fdpa/T1ImgSegment/S5/wc3coSPGR.img';
cfg{k}.save_path = '/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/bramila/testdata/S5';

k=3;
cfg{k}.infile = '/scratch/braindata/kauttoj2/Tulitikku/new_fdpa/FunImgNormDetrendedArtrepSmoothedNii/S6/4D.nii';
cfg{k}.motionparam = '/scratch/braindata/kauttoj2/Tulitikku/new_fdpa/RealignParameter/S6/rp_epi_0003.txt';
cfg{k}.wm_mask = '/scratch/braindata/kauttoj2/Tulitikku/new_fdpa/T1ImgSegment/S6/wc2coSPGR.img';
cfg{k}.csf_mask = '/scratch/braindata/kauttoj2/Tulitikku/new_fdpa/T1ImgSegment/S6/wc3coSPGR.img';
cfg{k}.save_path = '/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/bramila/testdata/S6';

[adj,roits,cfg,group_cfg] = bramila_taskFCpipeline(cfg,group_cfg);