function icatb_postprocess_timecourses(param_file)
%% Write FNC and spectra information in *_postprocess_results.mat.
% FNC correlations (transformed to fisher z-scores) are saved as fnc_corrs_all with
% dimensions subjects x sessions x components x components. Spectra is
% saved as spectra_tc_all of dimensions subjects x sessions x spectral length.
%

%% Load defaults
icatb_defaults;
global EXPERIMENTAL_TR;
global TIMECOURSE_POSTPROCESS;
global DETRENDNUMBER;
global PARAMETER_INFO_MAT_FILE;

filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', filterP);
    drawnow;
end

if (isempty(param_file))
    error('ICA parameter file is not selected');
end

if (ischar(param_file))
    load(param_file);
    
    if ~exist('sesInfo', 'var')
        error('Not a valid parameter file');
    end
    outputDir = fileparts(param_file);
else
    sesInfo = param_file;
    outputDir = sesInfo.outputDir;
end

if (isempty(outputDir))
    outputDir = pwd;
end

try
    modalityType = sesInfo.modality;
catch
    modalityType = icatb_get_modality;
end

if (~strcmpi(modalityType, 'fmri'))
    warning('!!!Timecourse post-processing will be done only for fMRI modality');
    return;
end

writeInfo = 0;
try
    writeInfo = TIMECOURSE_POSTPROCESS.write;
catch
end

if (~writeInfo)
    return;
end

if (isempty(EXPERIMENTAL_TR))
    warning('!!!Experimental TR variable (EXPERIMENTAL_TR) in seconds is missing.');
    return;
end

TR = EXPERIMENTAL_TR;

% Spectra info
tapers = [3, 5];
sampling_frequency = 1/TR;
frequency_band = [0, 1/(2*TR)];

% FNC (despike timecourses and High freq cutoff in Hz)
despike_tc = 1;
cutoff_frequency = 0.15;

%% Write results
if (writeInfo)
    
    % SEPCTRA PARAMS (tapers, sampling_freq, frequency_band)
    try
        tapers = TIMECOURSE_POSTPROCESS.spectra.tapers;
    catch
    end
    
    try
        sampling_frequency = TIMECOURSE_POSTPROCESS.spectra.sampling_frequency;
    catch
    end
    
    try
        frequency_band = TIMECOURSE_POSTPROCESS.spectra.frequency_band;
    catch
    end
    
    % FNC PARAMs (Despike timecourses and high frequency cutoff in Hz)
    try
        despike_tc = TIMECOURSE_POSTPROCESS.fnc.despike_tc;
    catch
    end
    
    try
        cutoff_frequency = TIMECOURSE_POSTPROCESS.fnc.cutoff_frequency;
    catch
    end
    
    outputFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_postprocess_results.mat']);
    
    %% Uncompress files
    subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess, 'flagTimePoints', ...
        sesInfo.flagTimePoints);
    sesInfo.outputDir = outputDir;
    fileIn = dir(fullfile(outputDir, [sesInfo.calibrate_components_mat_file, '*.mat']));
    filesToDelete = {};
    if (length(fileIn) ~= sesInfo.numOfSub*sesInfo.numOfSess)
        disp('Uncompressing subject component files ...');
        filesToDelete = icatb_unZipSubjectMaps(sesInfo, subjectICAFiles);
    end
    
    %% Spectra
    spectra_params = struct('tapers', tapers, 'Fs', sampling_frequency, 'fpass', frequency_band);
    countS = 0;
    fprintf('\nComputing spectra of all subjects and sessions components ...');
    for nSub = 1:sesInfo.numOfSub
        for nSess = 1:sesInfo.numOfSess
            countS = countS + 1;
            timecourses = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'subjects', nSub, 'sessions', nSess, 'vars_to_load', 'tc', 'truncate_tp', 1, 'detrend_no', ...
                DETRENDNUMBER, 'subject_ica_files', subjectICAFiles);
            timecourses = timecourses(1:min(sesInfo.diffTimePoints), :);
            [temp_spectra, freq] = icatb_get_spectra(timecourses', TR, spectra_params);
            temp_spectra = temp_spectra./repmat(sum(temp_spectra, 2), [1, size(temp_spectra, 2)]);
            temp_spectra = temp_spectra';
            if (countS == 1)
                spectra_tc_all = zeros(sesInfo.numOfSub, sesInfo.numOfSess, size(temp_spectra, 1), sesInfo.numComp);
            end
            spectra_tc_all(nSub, nSess, :, :) = temp_spectra;
        end
    end
    
    if (sesInfo.numOfSub*sesInfo.numOfSess == 1)
        spectra_tc_all = squeeze(spectra_tc_all);
    end
    
    save(outputFile, 'spectra_tc_all', 'freq');
    
    %% FNC
    fprintf('\nComputing FNC correlations ...\n');
    
    if (despike_tc)
        disp('Timecourses will be despiked ...');
    end
    
    if (cutoff_frequency > 0)
        disp(['Timecourses will be filtered using HF cutoff of ', num2str(cutoff_frequency), ' Hz ...']);
    end
    
    countS = 0;
    fprintf('\n');
    for nSub = 1:sesInfo.numOfSub
        for nSess = 1:sesInfo.numOfSess
            countS = countS + 1;
            timecourses = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'subjects', nSub, 'sessions', nSess, 'vars_to_load', 'tc', 'truncate_tp', 1, 'detrend_no', ...
                DETRENDNUMBER, 'subject_ica_files', subjectICAFiles);
            for nC = 1:sesInfo.numComp
                tmp = timecourses(:, nC);
                if (despike_tc)
                    tmp = icatb_despike_tc(tmp, TR);
                end
                if (cutoff_frequency > 0)
                    tmp = icatb_filt_data(tmp, TR, cutoff_frequency);
                end
                timecourses(:, nC) = tmp;
            end
            if (countS == 1)
                fnc_corrs_all = zeros(sesInfo.numOfSub, sesInfo.numOfSess, sesInfo.numComp, sesInfo.numComp);
            end
            c = icatb_corr(timecourses);
            c(1:size(c, 1) + 1:end) = 0;
            c = icatb_r_to_z(c);
            fnc_corrs_all(nSub, nSess, :, :) = c;
        end
    end
    
    if (sesInfo.numOfSub*sesInfo.numOfSess == 1)
        fnc_corrs_all = squeeze(fnc_corrs_all);
    end
    
    save(outputFile, 'fnc_corrs_all', '-append');
    
    fprintf('Done\n\n');
    disp(['File ', outputFile, ' contains spectra and FNC correlations']);
    disp('1. spectra_tc_all - Timecourses spectra. spectra_tc_all variable is of dimensions subjects x sessions x spectral length x components');
    disp('2. fnc_corrs_all - FNC correlations transformed to fisher z-scores. fnc_corrs_all variable is of dimensions subjects x sessions x components x components');
    
    
    if (exist('filesToDelete', 'var') && ~isempty(filesToDelete))
        icatb_cleanupFiles(filesToDelete, outputDir);
    end
    
end
