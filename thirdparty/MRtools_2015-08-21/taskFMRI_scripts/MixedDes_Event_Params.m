function P = MixedDes_Event_Params(dirName)
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

if nargin == 0;
    dirName = 'EventAna_Standard';
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
P.cm.DestFold = ['functional/First_Level_Models/' dirName '/SPM_Ana'];

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

%%%%%%%%
P.ct.do = 1;
P.ct.MoveContrasts = 1;
P.ct.GroupConFold = ['GroupContrasts/' dirName];
P.ct.MinEvents = 10;  %%% the minimum number of events required for a contrast to be generated.

%%%%%%%%%  For createVec_aps by APS
c = 0;

c = c+1; % 1
P.ct.con(c).name = 'HCH_gt_F';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'hch'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

c = c+1; % 2
P.ct.con(c).name = 'LCH_gt_F';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'lch'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

c = c+1; % 3
P.ct.con(c).name = 'LCM_gt_F';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'lcm'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

c = c+1; % 4
P.ct.con(c).name = 'HCM_gt_F';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'hcm'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

c = c+1; % 5
P.ct.con(c).name = 'Rep_gt_F';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'rep'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

c = c+1; % 6
P.ct.con(c).name = 'HCH_gt_LCH';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'hch'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right.mid = 'lch';
P.ct.con(c).right.pre = '';
P.ct.con(c).right.post = '';

c = c+1; % 7
P.ct.con(c).name = 'AllHits_gt_F';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = '(hch|lch)'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

c = c+1; % 8
P.ct.con(c).name = 'AllMiss_gt_F';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = '(hcm|lcm)'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

c = c+1; % 9
P.ct.con(c).name = 'All_gt_F';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = '(hch|hcm|lch|lcm|rep)'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

c = c+1; % 10
P.ct.con(c).name = 'AllNovel_gt_F';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = '(hch|hcm|lch|lcm)'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right = [];

c = c+1; % 11
P.ct.con(c).name = 'AllNov_gt_Rep';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = '(hch|hcm|lch|lcm)'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right.mid = 'rep';
P.ct.con(c).right.pre = '';
P.ct.con(c).right.post = '';

c = c+1; % 12
P.ct.con(c).name = 'AllHits_gt_AllMiss';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = '(hch|lch)'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right.mid = '(hcm|lcm)';
P.ct.con(c).right.pre = '';
P.ct.con(c).right.post = '';

c = c+1; % 13
P.ct.con(c).name = 'HCH_gt_AllMiss';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'hch'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right.mid = '(hcm|lcm)';
P.ct.con(c).right.pre = '';
P.ct.con(c).right.post = '';

c = c+1; % 14
P.ct.con(c).name = 'HCH_gt_Rep';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = 'hch'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right.mid = 'rep';
P.ct.con(c).right.pre = '';
P.ct.con(c).right.post = '';

c = c+1; % 15
P.ct.con(c).name = 'AllHits_gt_Rep';
P.ct.con(c).WeightWithin = 0;
P.ct.con(c).BlockThresh = 15;
P.ct.con(c).left.mid = '(hch|lch)'; 
P.ct.con(c).left.pre = ''; 
P.ct.con(c).left.post = ''; 
P.ct.con(c).right.mid = 'rep';
P.ct.con(c).right.pre = '';
P.ct.con(c).right.post = '';


%%%%%%%%%  For createVec by DGM
% P.ct.do = 1;
% P.ct.MoveContrasts = 1;
% P.ct.GroupConFold = 'GroupContrasts/MixedDes_Event';
% P.ct.WeightContrasts = 15;
% P.ct.MinEvents = 10;  %%% the minimum number of events required for a contrast to be generated.
% 
% c = 0;
% 
% c = c+1; % 1
% P.ct.name{c} = 'HCH_gt_F';
% P.ct.Con{c}.left = {'hch'};
% P.ct.Con{c}.right = {'none'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 2
% P.ct.name{c} = 'F_gt_HCH';
% P.ct.Con{c}.left = {'none'};
% P.ct.Con{c}.right = {'hch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 3
% P.ct.name{c} = 'LCH_gt_F';
% P.ct.Con{c}.left = {'lch'};
% P.ct.Con{c}.right = {'none'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 4
% P.ct.name{c} = 'F_gt_LCH';
% P.ct.Con{c}.left = {'none'};
% P.ct.Con{c}.right = {'lch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 5
% P.ct.name{c} = 'LCM_gt_F';
% P.ct.Con{c}.left = {'lcm'};
% P.ct.Con{c}.right = {'none'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 6
% P.ct.name{c} = 'F_gt_LCM';
% P.ct.Con{c}.left = {'none'};
% P.ct.Con{c}.right = {'lcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 7
% P.ct.name{c} = 'HCM_gt_F';
% P.ct.Con{c}.left = {'hcm'};
% P.ct.Con{c}.right = {'none'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 8
% P.ct.name{c} = 'F_gt_HCM';
% P.ct.Con{c}.left = {'none'};
% P.ct.Con{c}.right = {'hcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 9
% P.ct.name{c} = 'Rep_gt_F';
% P.ct.Con{c}.left = {'rep'};
% P.ct.Con{c}.right = {'none'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 10
% P.ct.name{c} = 'F_gt_Rep';
% P.ct.Con{c}.left = {'none'};
% P.ct.Con{c}.right = {'rep'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 11
% P.ct.name{c} = 'HCH_gt_HCM';
% P.ct.Con{c}.left = {'hch'};
% P.ct.Con{c}.right = {'hcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 12
% P.ct.name{c} = 'HCM_gt_HCH';
% Con(c).STAT = 'T';
% P.ct.Con{c}.left = {'hcm'};
% P.ct.Con{c}.right = {'hch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 13
% P.ct.name{c} = 'HCH_gt_LCM';
% P.ct.Con{c}.left = {'hch'};
% P.ct.Con{c}.right = {'lcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 14
% P.ct.name{c} = 'LCM_gt_HCH';
% P.ct.Con{c}.left = {'lcm'};
% P.ct.Con{c}.right = {'hch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 15
% P.ct.name{c} = 'HCH_gt_LCH';
% P.ct.Con{c}.left = {'hch'};
% P.ct.Con{c}.right = {'lch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 16
% P.ct.name{c} = 'LCH_gt_HCH';
% Con(c).STAT = 'T';
% P.ct.Con{c}.left = {'lch'};
% P.ct.Con{c}.right = {'hch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 17
% P.ct.name{c} = 'HCM_gt_LCH';
% P.ct.Con{c}.left = {'hcm'};
% P.ct.Con{c}.right = {'lch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 18
% P.ct.name{c} = 'LCH_gt_HCM';
% P.ct.Con{c}.left = {'lch'};
% P.ct.Con{c}.right = {'hcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 19
% P.ct.name{c} = 'HCM_gt_LCM';
% P.ct.Con{c}.left = {'hcm'};
% P.ct.Con{c}.right = {'lch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 20
% P.ct.name{c} = 'LCM_gt_HCM';
% P.ct.Con{c}.left = {'lcm'};
% P.ct.Con{c}.right = {'hcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 21
% P.ct.name{c} = 'LCH_gt_LCM';
% P.ct.Con{c}.left = {'lch'};
% P.ct.Con{c}.right = {'lcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 22
% P.ct.name{c} = 'LCM_gt_LCH';
% P.ct.Con{c}.left = {'lcm'};
% P.ct.Con{c}.right = {'lch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 23
% P.ct.name{c} = 'AllHits_gt_F';
% P.ct.Con{c}.left = {'hch' 'lch'};
% P.ct.Con{c}.right = {'none'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 24
% P.ct.name{c} = 'F_gt_AllHits';
% P.ct.Con{c}.left = {'none'};
% P.ct.Con{c}.right = {'hch' 'lch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 25
% P.ct.name{c} = 'AllMiss_gt_F';
% P.ct.Con{c}.left = {'hcm' 'lcm'};
% P.ct.Con{c}.right = {'none'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 26
% P.ct.name{c} = 'F_gt_AllMiss';
% P.ct.Con{c}.left = {'none'};
% P.ct.Con{c}.right = {'hcm' 'lcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 27
% P.ct.name{c} = 'All_gt_F';
% P.ct.Con{c}.left = {'hch' 'hcm' 'lch' 'lcm' 'rep'};
% P.ct.Con{c}.right = {'none'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 28
% P.ct.name{c} = 'F_gt_All';
% P.ct.Con{c}.left = {'none'};
% P.ct.Con{c}.right = {'hch' 'hcm' 'lch' 'lcm' 'rep'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 29
% P.ct.name{c} = 'AllNov_gt_F';
% P.ct.Con{c}.left = {'hch' 'hcm' 'lch' 'lcm'};
% P.ct.Con{c}.right = {'none'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 30
% P.ct.name{c} = 'F_gt_AllNov';
% P.ct.Con{c}.left = {'none'};
% P.ct.Con{c}.right = {'hch' 'hcm' 'lch' 'lcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 31
% P.ct.name{c} = 'AllNov_gt_Rep';
% P.ct.Con{c}.left = {'hch' 'hcm' 'lch' 'lcm'};
% P.ct.Con{c}.right = {'rep'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 32
% P.ct.name{c} = 'Rep_gt_AllNov';
% Con(c).STAT = 'T';
% P.ct.Con{c}.left = {'rep'};
% P.ct.Con{c}.right = {'hch' 'hcm' 'lch' 'lcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 33
% P.ct.name{c} = 'AllHits_gt_AllMiss';
% P.ct.Con{c}.left = {'hch' 'lch'};
% P.ct.Con{c}.right = {'hcm' 'lcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 34
% P.ct.name{c} = 'AllMiss_gt_AllHits';
% Con(c).STAT = 'T';
% P.ct.Con{c}.left = {'hcm' 'lcm'};
% P.ct.Con{c}.right = {'hch' 'lch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 35
% P.ct.name{c} = 'HCH_gt_AllMiss';
% P.ct.Con{c}.left = {'hch'};
% P.ct.Con{c}.right = {'hcm' 'lcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 36
% P.ct.name{c} = 'AllMiss_gt_HCH';
% P.ct.Con{c}.left = {'hcm' 'lcm'};
% P.ct.Con{c}.right = {'hch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 37
% P.ct.name{c} = 'LCH_gt_AllMiss';
% P.ct.Con{c}.left = {'lch'};
% P.ct.Con{c}.right = {'hcm' 'lcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 38
% P.ct.name{c} = 'AllMiss_gt_LCH';
% P.ct.Con{c}.left = {'hcm' 'lcm'};
% P.ct.Con{c}.right = {'lch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 39
% P.ct.name{c} = 'HCM_gt_AllHits';
% P.ct.Con{c}.left = {'hcm'};
% P.ct.Con{c}.right = {'hch' 'lch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 40
% P.ct.name{c} = 'AllHits_gt_HCM';
% P.ct.Con{c}.left = {'hch' 'lch'};
% P.ct.Con{c}.right = {'hcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 41
% P.ct.name{c} = 'LCM_gt_AllHits';
% P.ct.Con{c}.left = {'lcm'};
% P.ct.Con{c}.right = {'hch' 'lch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 42
% P.ct.name{c} = 'AllHits_gt_LCM';
% P.ct.Con{c}.left = {'hch' 'lch'};
% P.ct.Con{c}.right = {'lcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 43
% P.ct.name{c} = 'HCH_gt_Rep';
% P.ct.Con{c}.left = {'hch'};
% P.ct.Con{c}.right = {'rep'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 44
% P.ct.name{c} = 'Rep_gt_HCH';
% P.ct.Con{c}.left = {'rep'};
% P.ct.Con{c}.right = {'hch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 45
% P.ct.name{c} = 'HCM_gt_Rep';
% P.ct.Con{c}.left = {'hcm'};
% P.ct.Con{c}.right = {'rep'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 46
% P.ct.name{c} = 'Rep_gt_HCM';
% P.ct.Con{c}.left = {'rep'};
% P.ct.Con{c}.right = {'hcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 47
% P.ct.name{c} = 'LCH_gt_Rep';
% P.ct.Con{c}.left = {'lch'};
% P.ct.Con{c}.right = {'rep'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 48
% P.ct.name{c} = 'Rep_gt_LCH';
% P.ct.Con{c}.left = {'rep'};
% P.ct.Con{c}.right = {'lch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 49
% P.ct.name{c} = 'LCM_gt_Rep';
% P.ct.Con{c}.left = {'lcm'};
% P.ct.Con{c}.right = {'rep'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 50
% P.ct.name{c} = 'Rep_gt_LCM';
% P.ct.Con{c}.left = {'rep'};
% P.ct.Con{c}.right = {'lcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 51
% P.ct.name{c} = 'AllHits_gt_Rep';
% P.ct.Con{c}.left = {'hch' 'lch'};
% P.ct.Con{c}.right = {'rep'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 52
% P.ct.name{c} = 'Rep_gt_AllHits';
% P.ct.Con{c}.left = {'rep'};
% P.ct.Con{c}.right = {'hch' 'lch'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 53
% P.ct.name{c} = 'AllMiss_gt_Rep';
% P.ct.Con{c}.left = {'hcm' 'lcm'};
% P.ct.Con{c}.right = {'rep'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];
% 
% c = c+1; % 54
% P.ct.name{c} = 'Rep_gt_AllMiss';
% P.ct.Con{c}.left = {'rep'};
% P.ct.Con{c}.right = {'hcm' 'lcm'};
% P.ct.ConPrefix{c} = [];
% P.ct.ConTrail{c}  = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%