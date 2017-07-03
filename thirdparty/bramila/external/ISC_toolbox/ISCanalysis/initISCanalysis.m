function handles = initISCanalysis(handles)

Params = initParams;

% description/tag of the project:
% PublicParams.dataDescription = 'movieData1';
% 
% % set paths:
% PublicParams.maskPath = '/ISCtoolbox/templates/';
% PublicParams.atlasPath = '/ISCtoolbox/templates/';
% PublicParams.dataDestination = '/ISCresults/';
% 
% % Next, set source data file names for subjectSource -field. Full path 
% % names must be given inside curly brackets (cell-format), where rows 
% % present subjects and columns sessions,
% % e.g (2 subjects and 2 scanning sessions):
% %  PublicParams.subjectSource = 
% %  [ {/fullpathname/subj1_ses1.nii}, {/fullpathname/subj1_ses2.nii}; 
% %    {/fullpathname/subj1_ses2.nii}, {/fullpathname/subj2_ses2.nii} ]
% 
% PublicParams.subjectSource{1,1} = '/fMRIdata/Session1/preprocessed1.nii';
% PublicParams.subjectSource{1,2} = '/fMRIdata/Session1/preprocessed2.nii';
% 
% % set data size for each scanning session in the form 
% % [x y z t]. If this field is commented away, sizes 
% % will be calculated automatically (in this case, files
% % are temporally loaded into the Matlab to obtain size information).
% %PublicParams.dataSize = [91 109 91 244];%91 109 91 382];
% 
% % set source file format either to 'mat' or 'nii':
% PublicParams.fileFormatSubj = 'nii';
% 
% %%%%%%%%%%%%%%%%%%
% % Settings regarding inter-subject synchronization:
% 
% % Select intersubject similarity criteria that are calculated:
% PublicParams.ssiOn = 0; % sign-similarity index
% PublicParams.nmiOn = 0; % mutual information based index
% PublicParams.corOn = 1; % Pearson's correlation based index (recommended!)
% PublicParams.kenOn = 0; % Kendall's W based index
% 
% % Determine parts of the analysis that will be performed:
% PublicParams.calcStandard = 1; % standard analysis
% PublicParams.calcStats = 0; % median, quartile, t and variance ISC maps
% PublicParams.calcCorMatrices = 0; % save full correlation matrices
% PublicParams.calcPhase = 0;
% 
% %%%%%%%%%%%%%%%%%%
% % Frequency settings:
% 
% % set data sampling frequency:
% PublicParams.samplingFrequency = (1/3.4);
% 
% % Set number of frequency subbands:
% PublicParams.nrFreqBands = 3;
% 
% %%%%%%%%%%%%%%%%%% 
% % Time-window settings:
% 
% % time-windowing on/off:
% PublicParams.winOn = 0;
% % Set window size and step size in samples. 
% % Note: if windowSize exceeds time-series length,
% % windowing is not used.
% PublicParams.windowSize = 36;
% PublicParams.windowStep = 9;
% 
% %%%%%%%%%%%%%%%%%%%%
% % Template settings:
% PublicParams.useTemplate = 1;
% 
% %%%%%%%%%%%%%%%%%%
% % Permutation test settings. 
% % To perform the following analysis, corOn must be set to 1.
% 
% % Set permutation parameters for frequency-specific ISC.
% % To save time and memory, calculation can be done parallel
% %  using several CPU:s (nrPermutationSets). Total 
% % number of realizations in the "null distribution"
% % will be: nrPermutationSets * nrPermutations (e.g. 10e7).
% PublicParams.nrPermutationSets = 1;
% PublicParams.nrPermutations = 10e5;
% 
% % Set permutation parameters for frequency-band comparison statistic. 
% % Total number of realization in the "null distribution" will be
% % nrPermutationsZPF * numberOfBrainVoxels
% PublicParams.nrPermutationsZPF = 25000; % e.g. 25000

handles.Pub = Params.PublicParams;
handles.validFlag = false;
handles.Priv = [];
handles.ParamsValid = false;