function P = MixedDes_Block_Params(dirName)
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

if nargin==0
    dirName = 'BlockAna_Standard';
end
P.PO = []; %%% Specify which Steps to run and in which order.
P.PO{end+1} = 'Bad_Vols(fn)';         % 1
P.PO{end+1} = 'Create_Model(fn);';    % 2
P.PO{end+1} = 'Make_Cons(fn);';       % 3

P.root = pwd;
P.DestFold = ['functional/First_Level_Models/' dirName];
P.TR = 2;  % Specify your TR here.
P.List = 'RunIDs.txt';
P.BadRuns = 'bad_runs.txt';
P.Runs = [1 2 3 4 5 6];
P.nVols = [127 127 127 127 127 127];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Screening Parameters.
P.bv.do = 1;
P.bv.SourceFold = 'functional/Standard_Preproc/3_8mm_Smoothed';
P.bv.SourcePrefix = 'ss_';
P.bv.SourceSuffix = '.nii';

P.bv.CheckGlobalSNR = 1;
P.bv.SNR_SourceFold = 'functional/Standard_Preproc/2_Normalize_To_MNI_epi/SNR_Images';
P.bv.SNR_prefix = 'Global_SNR_Report';
P.bv.SNR_suffix = '.txt';
P.bv.SNRthresh = 115;       %%% The global SNR threshold that must be met.  If a run has less than this value it will be junked.
                            %%% liberal value = 115  moderate value = 175 conservative value =  235
P.bv.RealignDir = 'functional/Standard_Preproc/1_SliceTimeCorrected';
P.bv.RealignPrefix = 'realignment_';
P.bv.GlobalSigThresh = 2.5; %% Units are Standard Deviations to assess outliers of global signal velocity
P.bv.MoveThresh = .75;      %% Units are mm/TR
P.bv.RotThresh = 1.5;       %% Units are degrees/TR
P.bv.BadVolThresh = 20;     %% The maximum number of Bad Volumes allowed before the entire run is junked.
P.bv.TotalMovementThresh = 5; %% maximum difference in Position in mm
P.bv.TotalRotationThresh = 5; %% maximum difference in Rotation in degrees
P.bv.LogFileName = 'AnalysisLog.txt';
P.bv.BadVolRegsName = 'ExtraRegressors.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Model and Estimate
P.cm.do = 1;
P.cm.SourceFold = 'functional/Standard_Preproc/3_8mm_Smoothed';
P.cm.SourcePrefix = 'ss_';
P.cm.SourceSuffix = '.nii';
P.cm.DestFold = 'functional/First_Level_Models/BlockAna_Standard/SPM_ana';

P.cm.Units = 'secs';
P.cm.MicrotimeRes = 30;    %% Changed this from 16 to 30 after fMRI course in Abq. new value = number of slices
P.cm.MicrotimeOnset = 15;  %% Changed this from 1 to 15 after fMRI course in Abq. new value = middle slice;
P.cm.basis_func = 'hrf';   %% Other options include 'hrf (with time derivative)' and 'hrf (with time and dispersion derivatives)'

P.cm.Volterra = 1;  %% 1 is not modeling volterra interactions

%%% Put together SPM.xGX (This has to do with global normalization which we are not using).
P.cm.iGXcalc = 'None';
P.cm.sGXcalc = 'mean voxel value';
P.cm.sGMsca = 'session specific';

%%% Put together SPM.xVi This has to do with using and AR(1) process to remove serial correlations.  We are not using this.
P.cm.AR = 'none';
% P.cm.AR = 'AR(1)';

%%% Put together SPM.xX.  This is the specification of the high pass filter. 128 is the reccomended default value.
%%% We use 260 since the experiment is also set up as a block design with
%%% 130 seconds between blocks of the same kind.  This requires a higher HPF.
P.cm.HP_filt = [260 260 260 260 260 260];
P.cm.OnsetsFile = 'RunInfo.csv';

%%% This has to do with temporal modulation.  We are not using this.                
P.cm.TempMod = 'none';

P.cm.addBadVolRegs = 1;
P.cm.BadVolRegs = 'ExtraRegressors.mat';

P.cm.DisableThresholdMasking = 1;
P.cm.ExplicitMask(1) = spm_vol(which('brainmask.nii'));

P.cm.addMotionRegressors = 0;
P.cm.MotionFold = 'functional/Standard_Preproc/1_SliceTimeCorrected';
P.cm.MotionPrefix = 'realignment_';

P.cm.ScreenTaskCorrMot = 0;
P.cm.TaskCorrThresh = .25;  %%% R^2 values for cutoff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Contrasts.

P.ct.do = 1;
P.ct.MoveContrasts = 1;
P.ct.GroupConFold = ['GroupContrasts/' dirName];
P.ct.MinEvents = 1;  %%% the minimum number of events required for a contrast to be generated.

c = 0;

c = c+1; % 1
P.ct.con(c).name = 'Novel_gt_F';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

c = c+1; % 2
P.ct.con(c).name = 'Repeated_gt_F';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Repeated'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

c = c+1; % 3
P.ct.con(c).name = 'Novel_gt_Repeated';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right.mid = 'Repeated'; 
P.ct.con(c).right.pre = ''; 
P.ct.con(c).right.post = ''; 

c = c+1; % 4
P.ct.con(c).name = 'All_gt_F';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = '(Novel|Repeated)'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

%%

c = c+1; % 1
P.ct.con(c).name = 'Repeated_gt_F_R1';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Repeated'; 
P.ct.con(c).left.pre = '^Sn\([1]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

c = c+1; % 1
P.ct.con(c).name = 'Repeated_gt_F_R2';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Repeated'; 
P.ct.con(c).left.pre = '^Sn\([2]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];


c = c+1; % 1
P.ct.con(c).name = 'Repeated_gt_F_R3';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Repeated'; 
P.ct.con(c).left.pre = '^Sn\([3]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];


c = c+1; % 1
P.ct.con(c).name = 'Repeated_gt_F_R4';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Repeated'; 
P.ct.con(c).left.pre = '^Sn\([4]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];


c = c+1; % 1
P.ct.con(c).name = 'Repeated_gt_F_R5';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Repeated'; 
P.ct.con(c).left.pre = '^Sn\([5]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];


c = c+1; % 1
P.ct.con(c).name = 'Repeated_gt_F_R6';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Repeated'; 
P.ct.con(c).left.pre = '^Sn\([6]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];
%%

c = c+1; % 1
P.ct.con(c).name = 'Novel_gt_F_R1';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = '^Sn\([1]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

c = c+1; % 1
P.ct.con(c).name = 'Novel_gt_F_R2';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = '^Sn\([2]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];


c = c+1; % 1
P.ct.con(c).name = 'Novel_gt_F_R3';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = '^Sn\([3]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];


c = c+1; % 1
P.ct.con(c).name = 'Novel_gt_F_R4';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = '^Sn\([4]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];


c = c+1; % 1
P.ct.con(c).name = 'Novel_gt_F_R5';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = '^Sn\([5]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];


c = c+1; % 1
P.ct.con(c).name = 'Novel_gt_F_R6';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = '^Sn\([6]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];
%%

c = c+1; % 1
P.ct.con(c).name = 'Novel_gt_Repeated_R1';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = '^Sn\([1]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right.mid = 'Repeated'; 
P.ct.con(c).right.pre = '^Sn\([1]).*'; 
P.ct.con(c).right.post = ''; 

c = c+1; % 1
P.ct.con(c).name = 'Novel_gt_Repeated_R2';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = '^Sn\([2]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right.mid = 'Repeated'; 
P.ct.con(c).right.pre = '^Sn\([2]).*'; 
P.ct.con(c).right.post = ''; 


c = c+1; % 1
P.ct.con(c).name = 'Novel_gt_Repeated_R3';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = '^Sn\([3]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right.mid = 'Repeated'; 
P.ct.con(c).right.pre = '^Sn\([3]).*'; 
P.ct.con(c).right.post = '';  


c = c+1; % 1
P.ct.con(c).name = 'Novel_gt_Repeated_R4';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = '^Sn\([4]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right.mid = 'Repeated'; 
P.ct.con(c).right.pre = '^Sn\([4]).*'; 
P.ct.con(c).right.post = ''; 


c = c+1; % 1
P.ct.con(c).name = 'Novel_gt_Repeated_R5';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = '^Sn\([5]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right.mid = 'Repeated'; 
P.ct.con(c).right.pre = '^Sn\([5]).*';  
P.ct.con(c).right.post = ''; 


c = c+1; % 1
P.ct.con(c).name = 'Novel_gt_Repeated_R6';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'Novel'; 
P.ct.con(c).left.pre = '^Sn\([6]).*'; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right.mid = 'Repeated'; 
P.ct.con(c).right.pre = '^Sn\([6]).*'; 
P.ct.con(c).right.post = ''; 
% cd /autofs/space/lyrica_003/users/R01_Long_SPM8/RLB008C4_MA_6m/functional/First_Level_Models/BlockAna_Standard/SPM_ana2
% tmp = dir_wfp('beta*.img'); tmp = tmp(contains('Repeated',SPM.xX.name)); m = FastRead(tmp,VOI
%%
%%%%%%% 
% P.ct.do = 1;
% P.ct.MoveContrasts = 1;
% P.ct.GroupConFold = 'GroupContrasts/MixedDes_Block';
% P.ct.WeightContrasts = 15;
% P.ct.MinEvents = 1;  %%% the minimum number of events required for a contrast to be generated.
% 
% P.ct.name{1} = 'Novel_gt_F';
% P.ct.Con{1}.left  = {'Novel'};
% P.ct.Con{1}.right = {'none'};
% P.ct.ConPrefix{1} = [];
% P.ct.ConTrail{1}  = [];
% 
% P.ct.name{2} = 'F_gt_Novel';
% P.ct.Con{2}.left  = {'none'};
% P.ct.Con{2}.right = {'Novel'};
% P.ct.ConPrefix{2} = [];
% P.ct.ConTrail{2}  = [];
% 
% P.ct.name{3} = 'Repeated_gt_F';
% P.ct.Con{3}.left  = {'Repeated'};
% P.ct.Con{3}.right = {'none'};
% P.ct.ConPrefix{3} = [];
% P.ct.ConTrail{3}  = [];
% 
% P.ct.name{4} = 'F_gt_Repeated';
% P.ct.Con{4}.left  = {'none'};
% P.ct.Con{4}.right = {'Repeated'};
% P.ct.ConPrefix{4} = [];
% P.ct.ConTrail{4}  = [];
% 
% P.ct.name{5} = 'Novel_gt_Repeated';
% P.ct.Con{5}.left  = {'Novel'};
% P.ct.Con{5}.right = {'Repeated'};
% P.ct.ConPrefix{5} = [];
% P.ct.ConTrail{5}  = [];
% 
% P.ct.name{6} = 'Repeated_gt_Novel';
% P.ct.Con{6}.left  = {'Repeated'};
% P.ct.Con{6}.right = {'Novel'};
% P.ct.ConPrefix{6} = [];
% P.ct.ConTrail{6}  = [];
% 
% P.ct.name{7} = 'All_gt_F';
% P.ct.Con{7}.left  = {'Novel' 'Repeated'};
% P.ct.Con{7}.right = {'none'};
% P.ct.ConPrefix{7} = [];
% P.ct.ConTrail{7}  = [];
% 
% P.ct.name{8} = 'F_gt_All';
% P.ct.Con{8}.left  = {'none'};
% P.ct.Con{8}.right = {'Novel' 'Repeated'};
% P.ct.ConPrefix{8} = [];
% P.ct.ConTrail{8}  = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
