function filename = bramila_compute_ISC_mask(cfg,group_cfg)
% bramila_compute_ISC_mask - ISC toolbox wrapper
%   - Usage:
%   filename = bramila_compute_ISC_mask(cfg,group_cfg)
%   - Input:
%       cfg = a valid subject-wise cfg file
%       group_cfg = a valid group_cfg file
%   - Output:
%       filename = computed ISC maskfile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis destination directory. This must be same as
% specified in the parameter file.
if ~exist(group_cfg.ISC_path,'dir')
    mkdir(group_cfg.ISC_path);
end
analysisDestinationPath = group_cfg.ISC_path;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory of the templates:
templatePath = group_cfg.ISC_toolbox_path;
% Directory of the codes for running the ISC analysis:
ISCtoolboxPath = group_cfg.ISC_toolbox_path;
% Directory of the codes to process nifti -files:
niftiPath =  [group_cfg.ISC_toolbox_path,'niftitools/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not edit this part:
if exist(templatePath) ~= 7
    warning(['MNI templates not found in ' templatePath '.'])
end
if exist(ISCtoolboxPath) ~= 7
    error(['ISC toolbox not found in ' ISCtoolboxPath '.'])
end
if exist(templatePath) ~= 7
    warning(['Destination directory ' analysisDestinationPath ' does not exist.'])
end
curDir = cd;
cd(ISCtoolboxPath) 
setISCToolboxPath;
cd(curDir)

% set path:
disp(' ');disp(' ')
P=path;
P=path(P,ISCtoolboxPath);
P=path(P,analysisDestinationPath);
P=path(P,niftiPath);
P=path(P,templatePath);

disp(' ');disp(' ')
% load parameters from the destination directory. 
% Then you can run analysis by typing: runAnalysis(Params).
try
  load([analysisDestinationPath 'Tag'])
  load([analysisDestinationPath Tag '.mat'])
  [Params,flag] = paramsValidation(Params);
  if ~flag
    disp('Current parameters loaded to workspace.');disp(' ')
    disp('If you fill in parameters manually, validate parameters by typing: "Params=paramsValidation(Params);"');disp(' ')
    disp('Type "ISCanalysis" to use GUI to set parameters.');disp(' ')    
  end            
catch
    warning(['No parameter struct found in ' analysisDestinationPath]);disp(' ')
    disp('Default parameters loaded to workspace.');
    Params = initParams(templatePath,analysisDestinationPath,group_cfg,cfg);
    [Params,flag] = paramsValidation(Params);    
    if ~flag
        disp('If you fill in parameters manually, Type "Params.PublicParams" to see them.');disp(' ')
        disp('Validate parameters by typing: "Params=paramsValidation(Params);"');disp(' ')
        disp('Type "ISCanalysis" to use GUI to set parameters!')
    end
end

disp(' ')

if ~flag
    error('Cannot proceed to analysis!')
end

runAnalysis(Params);

load([analysisDestinationPath,'memMaps.mat']);
corvals = (memMaps.resultMap.whole.band0.Session1.cor.Data.xyz);
corvals(isnan(corvals))=-1;
load([analysisDestinationPath,'/stats/Thband0Session1win0.mat'],'Th','Th_info');
fdr_005 = Th(3);
mask = 0*corvals;
mask(corvals>=fdr_005)=1;
mask = bramila_clusterextend(mask,20);
save_nii(make_nii(mask,[2,2,2]),[analysisDestinationPath,'ISC_mask_fdr005_clustermin20.nii']);
filename=[analysisDestinationPath,'ISC_mask_fdr005_clustermin20.nii'];

fprintf('ISC mask has %i voxels\n',nnz(mask));

end

function Params = initParams(atlasPath,destinationPath,group_cfg,cfg)
% This function initializes parameters for intersubject 
% synchronization analysis. User must specify all fields in 
% PublicParams before performing the analysis (see runAnalysis.m).
%
% inputs:
% atlasPath - directory containing templates and  mask
% destinationPath - directory where analysis results will be stored

% Example usage:
% Set fields in runInit and initParams and then call:
% runInit;Params = initParams;
%
% See also: RUNANALYSIS

% Last updated: 5.8.2013 by Jukka-Pekka Kauppi
% University of Helsinki
% Department of Computer Science
% e-mail: jukka-pekka.kauppi@helsinki.fi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PUBLIC PARAMETERS
% These parameters must be set correctly before running the analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
% General settings:


% description/tag of the project:
PublicParams.dataDescription = 'Bramila_tools_ISC';

PublicParams.atlasPath = atlasPath;
PublicParams.dataDestination = destinationPath;

% Next, set source data file names for subjectSource -field. Full path 
% names must be given inside curly brackets (cell-format), where rows 
% present subjects and columns sessions,
% e.g (2 subjects and 2 scanning sessions):
%  PublicParams.subjectSource = 
%  [ {/fullpathname/subj1_ses1.nii}, {/fullpathname/subj2_ses1.nii}; 
%    {/fullpathname/subj1_ses2.nii}, {/fullpathname/subj2_ses2.nii} ]

for i=1:length(cfg)
    PublicParams.subjectSource{1,i}  = cfg{i}.infile;
end  

% set subject source file format either to 'mat' or 'nii':
PublicParams.fileFormatSubj = 'nii';

%%%%%%%%%%%%%%%%%%
% Settings regarding inter-subject synchronization:

% Select intersubject similarity criteria that are calculated:
PublicParams.ssiOn = 0; % sign-similarity index
PublicParams.nmiOn = 0; % mutual information based index
PublicParams.corOn = 1; % Pearson's correlation based index (recommended!)
PublicParams.kenOn = 0; % Kendall's W based index

% Determine parts of the analysis that will be performed:
PublicParams.calcStandard = 1; % standard analysis
PublicParams.calcStats = 0; % median, quartile, t and variance ISC maps
PublicParams.calcCorMatrices = 0; % save full correlation matrices
PublicParams.calcPhase = 0;

% use MNI template (FSL nifti)
PublicParams.useTemplate = 1;

PublicParams.permutSessionComp = 25000;
PublicParams.sessionCompOn = 0;

%%%%%%%%%%%%%%%%%%
% Frequency settings:

% set data sampling frequency:
PublicParams.samplingFrequency = group_cfg.TR;

% Set number of frequency subbands:
PublicParams.nrFreqBands = 0;
PublicParams.permutFreqComp = 25000;
PublicParams.freqCompOn = 1;
%%%%%%%%%%%%%%%%%% 
% Time-window settings:

% time-windowing on/off:
PublicParams.winOn = 0;
% Set window size and step size in samples. 
% Note: if windowSize exceeds time-series length,
% windowing is not used.
PublicParams.windowSize = 30;
PublicParams.windowStep = 30;

%%%%%%%%%%%%%%%%%%
% Permutation test settings. 
% To perform the following analysis, corOn must be set to 1.

% Set permutation parameters for frequency-specific ISC.
% To save time and memory, calculation can be done parallel
%  using several CPU:s (nrPermutationSets). Total 
% number of realizations in the "null distribution"
% will be: nrPermutationSets * nrPermutations (e.g. 10e7).
PublicParams.nrPermutationSets = 1;
PublicParams.nrPermutations = 1e6;

% Set permutation parameters for frequency-band comparison statistic. 
% Total number of realization in the "null distribution" will be
% nrPermutationsZPF * numberOfBrainVoxels
%PublicParams.nrPermutationsZPF = 25000; % e.g. 25000

%%%%%%%%%%%%%%%%
% Grid computation settings:
PublicParams.disableGrid = true;

%%%%%%%%%%%%%%%%
% Memmap removal after analysis settings:
PublicParams.removeMemmaps = false;
PublicParams.removeFiltermaps = true;

Params.PublicParams = PublicParams;

end
