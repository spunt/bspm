function runAnalysis(Params)

% This is the main function that performs ISC analysis for 
% preprocessed fMRI data. Analysis parameters (see initParams.m) 
% must be set by the user before running the analysis.
%
% Input fMRI data must be either in "nii" or "mat" -format 
% (see initParams.m). Output data is mapped into the memory and 
% can be quickly accessed through memory pointer objects.
%
% COMMAND LINE ANALYSIS:
% 
% After running the analysis, memory pointer objects can be loaded 
% to the Matlab's workspace by typing:
% load memMaps % objects are saved in analysis destination folder.
%
% GRAPHICAL USER INTERFACE (GUI):
% 
% Specific GUI has been designed to visualize and analyze the results.
% After running the analysis, GUI can be launched by typing:
% load Tag % Tag is saved in the analysis destination folder
% ISCtool(Tag); % launch GUI from your GUI folder
% 
% Note: when you launch GUI make sure you have cleared all memory
% map pointer objects from Matlab workspace.
% 
% See also: ISCANALYSIS, MEMMAPDATA, INITPARAMS

% Updated: 7.11.2013 by Jukka-Pekka Kauppi
% University of Helsinki
% Department of Computer Science
% email: jukka-pekka.kauppi@helsinki.fi
%
% Updated: 5.8.2013 by Juha Pajula
% Tampere University of Technology
% Department of Signal Processing
% e-mail: juha.pajula@tut.fi
%creating the log folder, this should be moved to params struct (and edited
%to setPrivParams.m as well as setDataPaths.m

total_time = tic;
if ~exist([Params.PublicParams.dataDestination,'scripts'],'dir')
    mkdir([Params.PublicParams.dataDestination,'scripts']); %create folder for scripts
end
log_path=[Params.PublicParams.dataDestination,'scripts'];

%saving the Log from command window:
diary([log_path, Params.PublicParams.dataDestination(end), 'main_log.txt']);
disp(datestr(now))
% Init run for grid
gridOff=Params.PublicParams.disableGrid;
if gridOff
    grid_type=''; %disable grid computing;
    disp('Grid computing disabled')
else
    grid_type=testGrid; %test if cluster environments (slurm/SGE) are here
end
if isempty(grid_type)
    deleteLockFiles(Params);
end


pauset = 5; %check interval for waitGrid in seconds


% STAGE 1:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize parameters. Note: user must set all public parameters
% in initParams.m before running the analysis.
%Params = initParams;
% initialize data matrices into the disk and create pointers
% using dynamic memory mapping:

if(~isempty(grid_type))
    [~,outpt]=gridParser('memMapData',Params,{},grid_type,'ISC_1');
    disp(outpt);
else
    memMapData(Params);
end
%if submitted to grid, waiting that all processes are finished before next
%Batch
if(~isempty(grid_type))
    waitGrid(grid_type,pauset,log_path)
end

%reload the params after memMapping
load([Params.PublicParams.dataDestination, Params.PublicParams.dataDescription])
Priv = Params.PrivateParams;
Pub = Params.PublicParams;

% STAGE 2:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter all data through stationary wavelet filter bank.
% Computationally faster way is to distribute the
% calculations across several processors.
disp(' ')
disp('Filtering:')
for nrSubject = 1:Priv.nrSubjects
  for nrSession = 1:Priv.nrSessions
    disp(['Subject: ' num2str(nrSubject) ' , Session: ' num2str(nrSession) ':'])
    if(~isempty(grid_type))
        [~,outpt]=gridParser('filterData',Params,{nrSubject,nrSession},grid_type,'ISC_2');
        disp(outpt);
    else
        filterData(Params,nrSubject,nrSession);
    end
  end
end

%if submitted to grid, waiting that all processes are finished before next
%Batch
if(~isempty(grid_type))
    waitGrid(grid_type,pauset,log_path)
end

% STAGE 3:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate intersubject synchronization maps.
% Computationally faster way is to distribute the calculations over several processors.
% Here, nrBand 0 stands for full frequency band.
disp(' ')
disp('Calculating inter-subject synchronization maps:')
for nrBand = 0:Pub.nrFreqBands
  for nrSession = 1:Priv.nrSessions
    disp(['Band: ' num2str(nrBand) ', Session: ' num2str(nrSession)])
    if(~isempty(grid_type))
        [~,outpt]=gridParser('calculateSimilarityMaps',Params,{nrBand,nrSession},grid_type,'ISC_3_1');
        disp(outpt);
    else
        calculateSimilarityMaps(Params,nrBand,nrSession);
    end
  end
end

%if submitted to grid, waiting that all processes are finished before next
if(~isempty(grid_type))
    waitGrid(grid_type,pauset,log_path)
end

% Compute resampling distribution for statistical si:
% Generate permutation (null) distributions:
disp(' ')
disp('Calculating permutation distributions for each ISC map:')
for nrSession = 1:Priv.nrSessions
  for permSetIdx = 1:Priv.nrPermutationSets
    if(~isempty(grid_type))
        [~,outpt1]=gridParser('permutationTest',Params,{nrSession,permSetIdx,0},grid_type,'ISC_3_2a');
        [~,outpt2]=gridParser('permutationTest',Params,{nrSession,permSetIdx,1},grid_type,'ISC_3_2b');
        disp(outpt1);
        disp(outpt2);
    else
        permutationTest(Params,nrSession,permSetIdx,0); % across session
        permutationTest(Params,nrSession,permSetIdx,1); % time-windows
    end
  end
end

%if submitted to grid, waiting that all processes are finished before next
if(~isempty(grid_type))
    waitGrid(grid_type,pauset,log_path)
end


% GRID
% Calculate session comparison maps (sum ZPF statistic).
disp(' ')
disp('Calculating sum ZPF statistic for session comparisons:')
% Get total number of session comparisons:
sessComps = ((Priv.nrSessions)^2-(Priv.nrSessions))/2;
%  Calculate null distributions for all comparisons:
for nrBand = 0:Priv.maxScale + 1
    for sessComp = 1:sessComps
        PearsonFilonAcrossSessions(Params,nrBand,sessComp);
    end
end

% GRID
% Assess statistical significance for session differences through
% permutation testing.
disp(' ')
disp('Calculating permutation distributions for session comparisons:')
for nrBand = 0:Priv.maxScale + 1
  for sessComp = 1:sessComps
    permutationPFAcrossSessions(Params,nrBand,sessComp);
  end
end




% Calculate frequency comparison maps (sum ZPF statistic).
% Get total number of frequency band comparisons:
disp(' ')
disp('Calculating sum ZPF statistic for frequency-band comparisons:')
freqComps = ((Priv.maxScale+2)^2-(Priv.maxScale+2))/2;
%  Calculate null distributions for all comparisons indexed from
% 1 to freqComps. To distribute computations, you can also call
% function  by giving subblocks as input, e.g. call function
% separately with comparisons 1:3, 4:6, 7:9, 10:12, and 13:15.
for nrSession = 1:Priv.nrSessions
  for freqComp = 1:freqComps
      if(~isempty(grid_type))
          [~,outpt]=gridParser('PearsonFilon',Params,{nrSession,freqComp},grid_type,'ISC_3_3');
          disp(outpt);
      else
          PearsonFilon(Params,nrSession,freqComp);
      end
  end
end

%if submitted to grid, waiting that all processes are finished before next
if(~isempty(grid_type))
    waitGrid(grid_type,pauset,log_path)
end

% Assess statistical significance between frequency bands through
% permutation testing:
disp(' ')
disp('Calculating permutation distributions for frequency band comparisons:')
for nrSession = 1:Priv.nrSessions
  for freqComp = 1:freqComps
    if(~isempty(grid_type))
        [~,outpt]=gridParser('permutationPF',Params,{nrSession,freqComp},grid_type,'ISC_3_4');
        disp(outpt);
    else
        permutationPF(Params,nrSession,freqComp)
    end
  end
end

%if submitted to grid, waiting that all processes are finished before next
if(~isempty(grid_type))
    waitGrid(grid_type,pauset,log_path)
end


% Calculate full intersubject correlation matrices. These calculations
% are required if subject-wise correlations are investigated 
% (average ISC is calculated using the function 
% calculateSimilarityMaps.m but it saves only the mean correlation values).
disp(' ')
disp('Calculating ISC matrices:')
for nrBand = 0:Pub.nrFreqBands
  for nrSession = 1:Priv.nrSessions
    disp(['Band: ' num2str(nrBand) ', Session: ' num2str(nrSession)])
    if(~isempty(grid_type))
        [~,outpt]=gridParser('calculateCorMats',Params,{nrBand,nrSession},grid_type,'ISC_3_5');
        disp(outpt);
    else
        calculateCorMats(Params,nrBand,nrSession);
    end
  end
end

%if submitted to grid, waiting that all processes are finished before next
if(~isempty(grid_type))
    waitGrid(grid_type,pauset,log_path)
end

% calculate other ISC maps (t-stats, median ISC map, percentile ISC maps):
disp(' ')
disp('Calculating other ISC maps:')
for nrBand = 0:Pub.nrFreqBands
  for nrSession = 1:Priv.nrSessions
    disp(['Band: ' num2str(nrBand) ', Session: ' num2str(nrSession)])
    if(~isempty(grid_type))
        [~,outpt]=gridParser('calculateStatsMaps',Params,{nrBand,nrSession},grid_type,'ISC_3_6');
        disp(outpt);
    else
        calculateStatsMaps(Params,nrBand,nrSession);
    end
  end
end

%if submitted to grid, waiting that all processes are finished before next
if(~isempty(grid_type))
    waitGrid(grid_type,pauset,log_path)
end

disp(' ')                                                                                                      
disp('Calculating inter-subject phase synchronization maps:')                                                           
for nrBand = 0:Pub.nrFreqBands                                                                                      
  for nrSession = 1:Priv.nrSessions                                                                                 
    disp(['Band: ' num2str(nrBand) ', Session: ' num2str(nrSession)])                                               
    if(~isempty(grid_type))
        [~,outpt]=gridParser('calculatePhaseSynch',Params,{nrBand,nrSession},grid_type,'ISC_3_7');
        disp(outpt);
    else
        calculatePhaseSynch(Params,nrBand,nrSession);
    end
  end                                                                                                               
end

%if submitted to grid, waiting that all processes are finished before next
if(~isempty(grid_type))
    waitGrid(grid_type,pauset,log_path)
end

% STAGE 4:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate ISC map thresholds using the resampling distribution and FDR:
disp(' ')
disp('Calculating critical thresholds:')
for nrSession = 1:Priv.nrSessions
    if(~isempty(grid_type))
        [~,outpt1]=gridParser('calculateThresholds',Params,{nrSession,0},grid_type,'ISC_4_1a');
        [~,outpt2]=gridParser('calculateThresholds',Params,{nrSession,1},grid_type,'ISC_4_1b');
        disp(outpt1);
        disp(outpt2);
    else
        calculateThresholds(Params,nrSession,0); % across session
        calculateThresholds(Params,nrSession,1); % time-windows
    end
end
%if submitted to grid, waiting that all processes are finished before next
if(~isempty(grid_type))
    waitGrid(grid_type,pauset,log_path)
end

% GRID
% calculate session comparison thresholds according to maximal statistic:
calculateThresholdsPFsessions(Params);

% calculate frequency comparison thresholds according to maximal statistic:
if(~isempty(grid_type))
    [~,outpt]=gridParser('calculateThresholdsPF',Params,{},grid_type,'ISC_4_2');
    disp(outpt);
else
    calculateThresholdsPF(Params);
end
%if submitted to grid, waiting that all processes are finished before next
if(~isempty(grid_type))
    waitGrid(grid_type,pauset,log_path)
end

% STAGE 5:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate intersubject synchronization curves from time-windowed data:
disp(' ')
disp('Calculating inter-subject synchronization curves:')
for nrBand = 0:Pub.nrFreqBands
  for nrSession = 1:Priv.nrSessions
    disp(['Band: ' num2str(nrBand) ', Session: ' num2str(nrSession)])
    if(~isempty(grid_type))
        [~,outpt]=gridParser('calculateSynchCurves',Params,{nrBand,nrSession},grid_type,'ISC_5');
        disp(outpt);
    else
        calculateSynchCurves(Params,nrBand,nrSession);
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if submitted to grid, waiting that all processes are finished before next
if(~isempty(grid_type))
    waitGrid(grid_type,pauset,log_path)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STAGE 6 (experimental):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results to Nii-files
disp(' ')
disp('Saving Nii')

saveNiiResults(Params)

% BATCH 7 Memmapped sourcedata removal:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% removing the memmapped source data to save disk space, if this Batch is
% selected the data must be memmapped again to be able to run the analysis
disp(' ')

removeMemmapData(Params)




disp(' ')
disp('Finished!!')

toc(total_time)
disp(datestr(now))
diary off