function P = Preproc_Params
%%% Written by Aaron Schultz on Feb. 5th, 2011 (aschultz@martinos.org).
%%% Choose which opperations to perform and the order in which they are to be executed.
%%%
%%% All filenames must end with an accurate run number.
%%%
%%% Orgamization Tips:
%%% This script makes three assumptions
%%% 1. All Sessions are located within a single folder.
%%% 2. All Runs for each Session are 4D Nifti files that are all within one
%%%    folder.
%%% 3. The name of the nifti files end with _# (e.g. _1.nii , _01.nii, 
%%%    _001.nii, etc.).  
%%%
%%% As a general rule of thumb you will want the numberinf convention to be 
%%% such that the matlab dir commmand will return the file names in the
%%% correct run order.  If this is not the case the script may fail or
%%% possibly mix up the run orders.
%%%
%%% Before you can run this script you must manually create the P.DestFold
%%% folder and create the P.List file with the Run Identifiers.
%%% Example RunIDs.txt:
%%% Run_1
%%% Run_2
%%% Run_3
%%% Run_4
%%% Run_5
%%% Run_6
%%%
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

P.PO = []; %%% Specify which Steps to run and in which order.
P.PO{end+1} = 'nfn = slice_time(nfn);';         % 1
P.PO{end+1} = 'nfn = realign(nfn);';            % 2
P.PO{end+1} = 'nfn = normalize(nfn);';          % 3
P.PO{end+1} = 'nfn = smooth(nfn);';             % 4

P.TR = 2.0;  % Specify your TR here.
P.root = pwd;
P.DestFold = 'functional/Standard_Preproc';
P.List = 'RunIDs.txt';
P.Runs = [1 2 3 4 5 6];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slice Timing Paramters: See SPM8 documentation for more info.
P.st.do = 1;  %%% 1 =  do slice timing;  0 =  skip slice timing;
P.st.SourceFold = 'functional/Orig';
P.st.SourcePrefix = '*un';
P.st.SourceSuffix = '.nii';
P.st.DestFold = [P.DestFold '/1_SliceTimeCorrected'];

P.st.sliceorder = [2:2:30 1:2:29]; % This is the order in which the slices were acquired
% P.st.sliceorder = [1:2:47 2:2:46]; % This is the order in which the slices were acquired
P.st.refslice = 2;
P.st.timing = [(P.TR/numel(P.st.sliceorder)) (P.TR/numel(P.st.sliceorder))-((P.TR/numel(P.st.sliceorder))/numel(P.st.sliceorder))];
P.st.prefix = 'st_';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Realign and Reslice Parameters
P.rr.do = 1;  %%% 1 =  do;  0 =  skip;
P.rr.SourceFold = P.st.DestFold;
P.rr.SourcePrefix = 'st_*un';
P.rr.SourceSuffix = '.nii';
P.rr.DestFold = P.st.DestFold;
P.rr.RealignSeparate = 1;  % 1 = do separate realignments;  2 = relaign first images then align each run to to first image of each run
P.rr.RealignmentType = 1;  % 1 = use inria_realign;         2 = use spm_realign

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
P.rr.RealignParsS.quality = .95;  % Quality
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
P.rr.ReslicePars.which  = 3; %% Reslice all images including the first image
P.rr.ReslicePars.prefix = 'rr_';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization Parameters: See SPM8 documentation for more info.
P.nn.do = 1; %%% 1 =  do;  0 =  skip;
P.nn.SourceFold = P.st.DestFold;
P.nn.SourcePrefix = 'st_*un';
P.nn.SourceSuffix = '.nii';
P.nn.DestFold = [P.DestFold '/2_Normalize_To_MNI_epi'];

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
P.nn.rflags.preserve = 0;
P.nn.rflags.bb = [-78 -112 -70; 78 76 90]; %%[[-90 -126 -72];[ 90 90 108]] % Buckner BB
P.nn.rflags.vox = [3 3 3];
P.nn.rflags.interp = 5; % Buckner uses 7 order interp
P.nn.rflags.wrap = [0 0 0];
P.nn.rflags.prefix =  'nn_';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smoothing Parameters % 1 = Perform Smoothing. 0 = Skip it.
P.ss.do = 1; %%% 1 =  do;  0 =  skip;

P.ss.SourceFold = P.nn.DestFold;
P.ss.SourcePrefix = 'nn_';
P.ss.SourceSuffix = '.nii';
P.ss.DestFold = [P.DestFold '/3_8mm_Smoothed'];

P.ss.prefix = 'ss_';
P.ss.kernel = [8 8 8]; % 3D smoothing kernel in mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
