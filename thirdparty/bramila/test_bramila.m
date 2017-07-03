% Bramila pipeline demo for one subject
% 
% NOTE: full 4D EPI data must fit into memory

clc;
close all;
clear cfg;

% REQUIRED processing parameters
cfg.bramilapath = '/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/bramila'; % bramila path
cfg.TR = 2.015;

% OPTIONAL parameters (otherwise use defaults)
% Modify only if needed (see "bramila_checkparameters.m" for defaults)
% cfg.motion_reg_type = 'friston'; % motion regression type
% cfg.voxelsize=[2,2,2];        % voxelsize
% cfg.mask = [];                % initial EPI mask
% cfg.tissue_derivatives = 0;   % tissue regressor derivative order
% cfg.min_tissue_var_expl = 75; % minimum variance percentage for tissue PCA regressors
% cfg.max_tissue_pca_count = 7; % upper limit for tissue PCA regressors
% cfg.remove_global = 0;        % remove global signal
% cfg.mot_derivatives = 1;  % motion regressor derivatives
% cfg.white_mask_th = 0.90; % probability threshold for white matter tissue
% cfg.csf_mask_th = 0.90;   % probability threshold for ventricle tissue
% cfg.detrend_type='linear-nodemean';   % detrending type
% cfg.filter_type = 'butter';   % temporal filter type
% cfg.filter_limits = [0,0.01,0.08,0.09];   % filter limits in Hz
% cfg.write = 0;    % write all intermediate EPI's (otherwise only the final network)
% cfg.network_nodes = []; % custom roi definition file that include proper 'rois' variable
% cfg.save_path = '/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/bramila/out'; % custom save path

% subject 1 EPI + motion + (tissues)
cfg.infile = '/triton/becs/scratch/braindata/kauttoj2/code/bramila/sampledata/p18/P18_EPI_4D.nii';
cfg.motionparam = '/triton/becs/scratch/braindata/kauttoj2/code/bramila/sampledata/p18/P18_movement.txt';
%cfg.white_mask = '/scratch/braindata/kauttoj2/code/bramila/sampledata/p15/P15_tissue_wm.img';
%cfg.csf_mask = '/scratch/braindata/kauttoj2/code/bramila/sampledata/p15/P15_tissue_csf.img';
[adj,roits,cfg] = bramila_restFCpipeline(cfg);

% custom save path subject 1 EPI + motion
cfg.infile = '/triton/becs/scratch/braindata/kauttoj2/code/bramila/sampledata/p18/P18_EPI_4D.nii';
cfg.motionparam = '/triton/becs/scratch/braindata/kauttoj2/code/bramila/sampledata/p18/P18_movement.txt';
cfg.save_path = '/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/bramila/out';
[adj,roits,cfg] = bramila_restFCpipeline(cfg);
