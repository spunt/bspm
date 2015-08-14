%% Editable GIPC variables

% NOTE! All time is measured in seconds unless stated otherwise.

% FDR (False Discovery Rate) Corrected or uncorrected
% If uncorrected the p-threshholds are uncorrected.
% (If group has too few subjects FDR-correction may sort away 
% all results and therefore is better switched to uncorrected to get any results)
gIpc.bFdr = true;

% Seconds per TR
gIpc.dTr = 1.5; %seconds

% Factor to interpolate TR more syncronized with seconds 
% Higher factor makes TR more interpolated to seconds and will find better
% match to onsets.
% However there are factors that interpolates TR to match the onsets exactly.
gIpc.dSyncFact = 20;

% Units used in the onset file
gIpc.bOnsetInTr0Sec1 = 1;

% Memory splitter - higher integers makes ipctb run with less memory (and slower).
gIpc.nSplit = 1;

% Number of clusters for cluster images
gIpc.nClusters = 5;

% Brain slices in Z direction
gIpc.ronSloverSlices = 72:-6:-36;

% Onsets of 3 different types for 2 sessions/runs
% Symbol used in graph, Line type and color (':k'=dotted black), 
% filename of tab/commaseparated onset file in the onset directory of GIPC
gIpc.nSess(1).susOnset = {'S',':k','soa_standards_run1.txt';...
                          'N',':k','soa_novels_run1.txt';...
                          'T',':r','soa_targets_run1.txt'};
gIpc.nSess(2).susOnset = {'S',':k','soa_standards_run2.txt';...
                          'N',':k','soa_novels_run2.txt';...
                          'T',':r','soa_targets_run2.txt'};    

% Timelock Graph titles
gIpc.sTitleGrp1 = 'HC';
gIpc.sTitleGrp2 = 'SZ';
% Onset type to lock to
gIpc.sOnsLockType = 'T';
% Time after lock that will be graphed
gIpc.dLockTime = 10; %seconds

% Cluster colors (first has to be 0 0 0).
gIpc.cmClust = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; .75 .5 1; .25 .5 .25];

% Defaults for p-val threshholds
gIpc.dFdrQ1Samp = 0.005; %.00005;
gIpc.dFdrQ2Samp = 0.9865; %.005;

% Default directory for output data
if isunix
    gIpc.sDirSave='/myGipcOutput/'; % unix Directory where you'd like to save your new data.
else
    gIpc.sDirSave='C:\myGipcOutput\'; %PC Directory where you'd like to save your new data.
end

% 2 sample T-test order
% if b2SampT_Grp1Min2 = false it means grp2 minus grp1 instead
gIpc.b2SampT_Grp1Min2 = true;
% Filename for 2 sample ttest for Grp1 > Grp1
if gIpc.b2SampT_Grp1Min2
    gIpc.sFileT2 = 'Grp1Minus2T_Thresh';
    gIpc.sFileT2Unthresh = 'Grp1Minus2T_UnthreshUncorr';
else
    gIpc.sFileT2 = 'Grp2Minus1T_Thresh';
    gIpc.sFileT2Unthresh = 'Grp2Minus1T_UnthreshUncorr';
end

% Names of resulting correlational image files
gIpc.sGrp1 = 'Group1T';
gIpc.sGrp2 = 'Group2T';

% Cluster file names including start of time courses
gIpc.sClust1 = 'ClustGrp1';
gIpc.sClust2 = 'ClustGrp2';  

% Default output directory for images
gIpc.sImageDir = 'Images';

% Directory for masks
gIpc.sMaskDir = 'mask';

% If all correlation maps are already saved 
% you may copy these files to an output directory 
% (put the correlation maps in their sub folders Grp1 and Grp1)
% plus you need to copy the settings.mat to the new output folder.
% Note! The source nifti-files has to remain at location they existed at
% original run!
% Now you can change all settings in this defaults file
% and p values in the user interface slashing processing time to 
% a fraction of what it was before
gIpc.bSkipStep1 = false;

%% Automatically set variables

% path to application
gIpc.sPathApp = '';
if strcmp(gIpc.sPathApp, '')
    gIpc.sPathApp = fileparts(mfilename('fullpath'));
end
gIpc.sPathApp = ipctb_backslash(gIpc.sPathApp); %run this even for user defined

% Directory for onset files (relative to main directory)
gIpc.sDirOnset = [gIpc.sPathApp ipctb_backslash('onsets')];

% Masks. Note if imported images have other dimensions 
% than 63*53*46 voxels the following supplied mask will have to be replaced
gIpc.sThreshMaskFile = [gIpc.sPathApp ipctb_backslash(gIpc.sMaskDir) 'sampleEPI.img'];
% gIpc.sThreshMaskFile = [gIpc.sPathApp ipctb_backslash(gIpc.sMaskDir) 'sampleEPI_alt.img'];

% AnatTemp
gIpc.sAnatTemp_dir = [gIpc.sPathApp ipctb_backslash(gIpc.sMaskDir) 'nsingle_subj_T1_2_2_3.img'];
% gIpc.sAnatTemp_dir = [gIpc.sPathApp ipctb_backslash(gIpc.sMaskDir) 'nsingle_subj_T1_2_2_3_alt.img'];

% SPM5 Directory
gIpc.sDirSpm = [gIpc.sPathApp ipctb_backslash('ipctb_spm5_files')];

% Set paths
sPath = which('gipc.m');
sPath = ipctb_backslash(fileparts(sPath));
addpath(sPath, '-end');
sPath2 = [sPath ipctb_backslash('ipctb_ica')];
addpath(sPath2, '-end');
sPath2 = [sPath ipctb_backslash('ipctb_spm5_files')];
addpath(sPath2, '-end');



%% Defaults originating from GIFT 102207
% colormap info
global COLORMAP_FILE;
global COLORLIST;

% conventions for naming components, timecourses and structural
global COMPONENT_NAMING;
global TIMECOURSE_NAMING;
global STRUCTURAL_NAMING;

% Matlab files
global PARAMETER_INFO_MAT_FILE; %Holds information for group session parameters

% File Indices
global SUBJECT_ICA_INDEX;
global MEAN_INDEX;
global TMAP_INDEX;
global STD_INDEX;
global MEAN_ALL_INDEX;  % Index for mean over different sessions

% Screen Color Defaults
global BG_COLOR; % Figure background color
global BG2_COLOR; % User Interface background color except pushbutton
global FG_COLOR;
global AXES_COLOR;
global FONT_COLOR; % User Interface foreground color except pushbutton
global LEGENDCOLOR;
global BUTTON_COLOR; % BUTTON background color
global BUTTON_FONT_COLOR; % BUTTON foreground color
global HELP_FONT_COLOR;

% Pictures
global COMPOSITE_VIEWER_PIC;
global COMPONENT_EXPLORER_PIC;
global ORTHOGONAL_VIEWER_PIC;

% display defaults
global SORT_COMPONENTS;
global IMAGE_VALUES;
global CONVERT_Z;
global THRESHOLD_VALUE;
global IMAGES_PER_FIGURE;
global COMPLEX_IMAGES_PER_FIGURE; % for complex data
global ANATOMICAL_PLANE;

% Slice Range defaults
global USE_DEFAULT_SLICE_RANGE;
global SLICE_RANGE; 

global ASPECT_RATIO_FOR_SQUARE_FIGURE;
global WS;
global MIN_SCREEN_DIM_IN_PIXELS;
global PERCENT_SCREEN_OCCUPIED;
global S0;

DIRS_TO_BE_STORED = {'C:\MATLAB6p5\work\Visuomotor_Data\visomot', 'C:\MATLAB6p5\work\Example Subjects'};

BG_COLOR = [.353 .471 .573]; % Figure background color
BG2_COLOR = [.353 .471 .573]; % User Interface controls background color except push button
FG_COLOR = [0 0 0];
FONT_COLOR = [1 1 1]; % User Interface controls foreground color except push button
AXES_COLOR = [0 0 0];
LEGENDCOLOR = [1 1 1];

%% DEBUG Settings for GIPC
dbstop if error;    %ce042408Takeaway
if isunix           %ce042408Takeaway
    bMinUi = true;  %ce042408Takeaway
else                %ce042408Takeaway
    bMinUi = false; %ce042408Takeaway
end                 %ce042408Takeaway