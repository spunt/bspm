function P = fcParams
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%% Copyright (C) 2014,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.

%%% Choose which opperations to perform and the order in which they are to be executed.
P.PO = []; %%% Specify which Steps to run and in which order.
P.PO{1} = 'nfn = drop_vols(fn)';
P.PO{2} = 'nfn = slice_time(nfn);';
P.PO{3} = 'nfn = realign(nfn);';
P.PO{4} = 'nfn = normalize(nfn);';
P.PO{5} = 'nfn = smooth(nfn);';
P.PO{6} = 'nfn = motion_regress(nfn);';
P.PO{7} = 'nfn = filter_data(nfn);';
P.PO{8} = 'nfn = nuisance_regress(nfn);';
P.PO{9} = 'Bad_Vols(nfn)';

P.TR = 3;  % Specify your TR here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drop Specified Volumes
P.dv.do = 1;
P.dv.dropVols = [1 2 3 4]; % All volumes specified here will be removed before analysis takes place.
P.dv.prefix = 'dv_';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slice Timing Paramters: See SPM8 documentation for more info.
P.st.do = 1;
P.st.sliceorder = [1:2:47 2:2:46]; % This is the order in which the slices were acquired
P.st.refslice = 1;
P.st.timing = [(P.TR/numel(P.st.sliceorder)) (P.TR/numel(P.st.sliceorder))-((P.TR/numel(P.st.sliceorder))/numel(P.st.sliceorder))];
P.st.prefix = 'st_';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Realign and Reslice Parameters
P.rr.do = 1;
%%% Realign: See SPM8 documentation for more info.

% Realign and Reslice Parameters
P.rr.do = 1;  %%% 1 =  do;  0 =  skip;
P.rr.RealignmentType = 2;  % 1 = use inria_realign;         2 = use spm_realign

%%% inria_realign specific parameters
P.rr.RealignParsI.quality = .95;  % Quality
P.rr.RealignParsI.fwhm    = 5;  % pre alignment smoothing kernel in mm.
P.rr.RealignParsI.sep     = 3;  % not 100% sure what this does.
P.rr.RealignParsI.interp  = 5;  % order of interpolation function for resampling.
P.rr.RealignParsI.wrap    = [0 0 0]; %% No wrapping.
P.rr.RealignParsI.rho_func = 'geman';
P.rr.RealignParsI.cutoff = 2.5;
P.rr.RealignParsI.hold = -9;

%%% spm_realign realignment paramters
P.rr.RealignParsS.quality = .50;  % Quality
P.rr.RealignParsS.fwhm    = 5;  % pre alignment smoothing kernel in mm.
P.rr.RealignParsS.sep     = 3;  % not 100% sure what this does.
P.rr.RealignParsS.interp  = 5;  % order of interpolation function for resampling.
P.rr.RealignParsS.wrap    = [0 0 0]; %% No wrapping.
P.rr.RealignParsS.rtm     = 0;  % no Register to mean (second pass alignment).
P.rr.RealignParsS.PW      = ''; % Weigthing?

%%% Reslice: See SPM8 documentation for more info.
P.rr.ReslicePars.mask   = 1; %% if one timepoint in a voxel is 0 all  timepoints in that voxel are set to 0;
P.rr.ReslicePars.mean   = 1; %% Write a Mean image.
P.rr.ReslicePars.interp = 5; %% interpolation order for resampling;
P.rr.ReslicePars.which  = 0; %% 0= No Reslice; 2 = Reslice all images including the first image
P.rr.ReslicePars.prefix = 'rr_';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization Parameters: See SPM8 documentation for more info.
P.nn.do = 1;
P.nn.source = 'mean*.nii';
tmp = which('spm');
ind = find(tmp==filesep);
whe = tmp(1:ind(end)-1);
P.nn.template = [whe '/templates/epi.nii'];

%%% Normalize data to MNI epi template: See SPM8 documentation for more info.
P.nn.NormPars.smosrc = 8;
P.nn.NormPars.smoref = 0;
P.nn.NormPars.regtype = 'mni';
P.nn.NormPars.cutoff = 25;
P.nn.NormPars.nits = 16;
P.nn.NormPars.reg = 1;
%%% Warping Paramters:  See SPM8 documentation for more info.
P.nn.rflags.skipThis = 0; %Set this equal to 1 to continue preprocessing in native space.
P.nn.rflags.preserve = 0;
P.nn.rflags.bb = [-78 -112 -70; 78 76 90]; %%[[-90 -126 -72];[ 90 90 108]] % Buckner BB
P.nn.rflags.vox = [3 3 3];
P.nn.rflags.interp = 5; % Buckner uses 7 order interp
P.nn.rflags.wrap = [0 0 0];
P.nn.rflags.prefix =  'nn_';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smoothing Parameters % 1 = Perform Smoothing. 0 = Skip it.
P.ss.do = 1;
P.ss.prefix = 'ss_';
P.ss.kernel = [6 6 6]; % 3D smoothing kernel in mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P.bv.do = 1;

P.bv.CheckGlobalSNR = 1;
P.bv.SNR_prefix = 'Global_SNR_Report';
P.bv.SNR_suffix = '.txt';
P.bv.SNRthresh = 115;       %%% The global SNR threshold that must be met.  If a run has less than this value it will be junked.
                            %%% liberal value = 115  moderate value = 175 conservative value =  235
% P.bv.RealignPrefix = 'realignment_';
P.bv.RealignPrefix = 'rp_';
P.bv.GlobalSigThresh = 2.5; %% Units are Standard Deviations to assess outliers of global signal velocity
P.bv.MoveThresh = .75;      %% Units are mm/TR
P.bv.RotThresh = 1.5;       %% Units are degrees/TR
P.bv.BadVolThresh = 20;     %% The maximum number of Bad Volumes allowed before the entire run is junked.
P.bv.MeanMovementThresh = .2; %% average change in position between TRs
P.bv.TotalMovementThresh = 5; %% maximum difference in Position in mm
P.bv.TotalRotationThresh = 5; %% maximum difference in Rotation in degrees
P.bv.LogFileName = 'AnalysisLog.txt';
P.bv.BadVolRegsName = 'ExtraRegressors.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motion Regression Params
P.mot.do = 1;
P.mot.MotionDeriv = 1;
P.mot.prefix = 'motRes_';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtering Parameters
P.ff.do = 1; % 1 = Perform Filtering. 0 = Skip it.
P.ff.Detrend = 0; % 1 = Perform Detrending. 0 = Skip it.
P.ff.HighCut = 0.08;
P.ff.LowCut = 0.01;
P.ff.FilterOrder = 4;
P.ff.prefix = 'ffBPS_';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nuisance Regressor Parameters
P.res.do = 1; % 1 = Perform the Nuisance variable regression. 0 = Skip it.
P.res.Motion = 0;  % 1 = Remove Variance attributable to motion. 0 = Do NOT remove variance attributable to motion.
P.res.MotionDeriv = 0; % 1 = Also regress out the first dirivitive of motion. 0 = Do NOT regress out the first dirivitive of motion.
P.res.Masks = {...
               which('avg152T1_brain_MNI.nii') ... % Masks used for computing physiological signals to be regressed out. Add as many as you would like.
               which('avg152T1_ventricles_MNI.nii') ...
               which('avg152T1_WM_MNI.nii') ...
               };
P.res.MaskDeriv = 1; % 1 = Also regress out the first dirivitive of physio vars. 0 = Do NOT regress out the first dirivitive of physio vars.
P.res.prefix = 'res_';

