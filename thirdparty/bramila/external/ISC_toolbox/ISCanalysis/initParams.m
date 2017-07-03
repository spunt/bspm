function Params = initParams(atlasPath,destinationPath)
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
PublicParams.dataDescription = 'Test';

if nargin < 2
    c = cd;
    if isunix
        sl = '/';
    else
        sl = '\';
    end
    c = [c sl];
    path0 = cd;
    path1 = [c 'templates'];
    PublicParams.atlasPath = path1;
    if nargin < 1
        path2 = [c 'results'];
        PublicParams.dataDestination = path2;
    else
        PublicParam.dataDestination = destinationPath;
    end
end
    if nargin == 2
        PublicParams.atlasPath = atlasPath;
        PublicParams.dataDestination = destinationPath;
    end
    if nargin > 2
        error('Number of inputs must be less than three!')
    end

% Next, set source data file names for subjectSource -field. Full path 
% names must be given inside curly brackets (cell-format), where rows 
% present subjects and columns sessions,
% e.g (2 subjects and 2 scanning sessions):
%  PublicParams.subjectSource = 
%  [ {/fullpathname/subj1_ses1.nii}, {/fullpathname/subj2_ses1.nii}; 
%    {/fullpathname/subj1_ses2.nii}, {/fullpathname/subj2_ses2.nii} ]

PublicParams.subjectSource{1,1}  = 'I1.mat';
PublicParams.subjectSource{1,2}  = 'I2.mat';    

% set subject source file format either to 'mat' or 'nii':
PublicParams.fileFormatSubj = 'mat';

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
PublicParams.sessionCompOn = 1;

%%%%%%%%%%%%%%%%%%
% Frequency settings:

% set data sampling frequency:
PublicParams.samplingFrequency = NaN;

% Set number of frequency subbands:
PublicParams.nrFreqBands = 2;
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
PublicParams.nrPermutations = 10e4;

% Set permutation parameters for frequency-band comparison statistic. 
% Total number of realization in the "null distribution" will be
% nrPermutationsZPF * numberOfBrainVoxels
%PublicParams.nrPermutationsZPF = 25000; % e.g. 25000

%%%%%%%%%%%%%%%%
% Grid computation settings:
PublicParams.disableGrid = false;

%%%%%%%%%%%%%%%%
% Memmap removal after analysis settings:
PublicParams.removeMemmaps = false;
PublicParams.removeFiltermaps = false;

Params.PublicParams = PublicParams;
