function icatb_run_mancovan(mancovanInfo, step)
%% Run Mancovan
%
% Inputs:
% 1. mancovanInfo - Mancovan information
%


icatb_defaults;
global MANCOVA_DEFAULTS;

if (~exist('mancovanInfo', 'var') || isempty(mancovanInfo))
    mancovanInfo = icatb_selectEntry('title', 'Select Mancovan Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', '*mancovan.mat');
    if (isempty(mancovanInfo))
        error('Mancovan parameter file is not selected');
    end
    drawnow;
end

if (ischar(mancovanInfo))
    mancovaParamFile = mancovanInfo;
    clear mancovanInfo;
    load(mancovaParamFile);
    outputDir = fileparts(mancovaParamFile);
    if (isempty(outputDir))
        outputDir = pwd;
    end
    mancovanInfo.outputDir = outputDir;
else
    mancovanInfo.outputDir = mancovanInfo.userInput.outputDir;
end

if (~exist('step', 'var'))
    step = 1;
end


desCriteria = 'mancova';
try
    desCriteria = mancovanInfo.designCriteria;
catch
end

mancovanInfo.designCriteria = desCriteria;

cd(mancovanInfo.outputDir);

logFile = fullfile(mancovanInfo.outputDir, [mancovanInfo.prefix, '_mancovan_results.log']);
tic;
diary(logFile);
mancovanInfo = mancovan_err_chk(mancovanInfo);
%mancovanInfo.numOfSub = mancovanInfo.userInput.numOfSub;
mancovanInfo.comp = mancovanInfo.userInput.comp;
doEst = 0;
if (isfield(mancovanInfo.userInput, 'doEstimation'))
    doEst = mancovanInfo.userInput.doEstimation;
end
mancovanInfo.doEstimation = doEst;
mancovanInfo.numOfPCs = mancovanInfo.userInput.numOfPCs;
mancovanInfo.features = mancovanInfo.userInput.features;
if (length(mancovanInfo.numOfPCs) == 1)
    mancovanInfo.numOfPCs = repmat(mancovanInfo.numOfPCs, 1, length(mancovanInfo.features));
    mancovanInfo.userInput.numOfPCs = mancovanInfo.numOfPCs;
end
TR = mancovanInfo.userInput.TR;
mancovanInfo.prefix = mancovanInfo.userInput.prefix;
try
    mancovanInfo.modelInteractions = mancovanInfo.userInput.modelInteractions;
catch
end
mancovanInfo.comps = [mancovanInfo.comp.value];
mancovanInfo.comps = mancovanInfo.comps(:)';
step_P = mancovanInfo.userInput.p_threshold;
if (~isfield(mancovanInfo.userInput, 'feature_params'))
    feature_params = icatb_mancovan_feature_options('tr', TR, 'mask_dims', mancovanInfo.userInput.HInfo(1).dim(1:3));
    out = icatb_OptionsWindow(feature_params.inputParameters, feature_params.defaults, 'off', 'title', 'Feature Options');
    feature_params.final = out.results;
    clear out;
else
    feature_params = mancovanInfo.userInput.feature_params;
end

def_mask_stat = 't';
def_mask_z_thresh = 1;
try
    def_mask_stat = lower(feature_params.final.stat_threshold_maps);
    def_mask_z_thresh = feature_params.final.z_threshold_maps;
catch
end

def_mask_t_std = 4;
try
    def_mask_t_std = MANCOVA_DEFAULTS.sm.t_std;
catch
end

% Only positive voxels
t_image_values = 2;

try
    t_image_values = MANCOVA_DEFAULTS.sm.t_image_values;
catch
end

% Filter cutoff frequency
filter_cutoff = 0;
if (ischar(feature_params.final.fnc_tc_filter))
    if (strcmpi(feature_params.final.fnc_tc_filter, 'yes'))
        filter_cutoff = 0.15;
    end
else
    filter_cutoff = feature_params.final.fnc_tc_filter;
end

covariates = 'none';
try
    covariates = mancovanInfo.userInput.feature_params.final.fnc_tc_covariates;
catch
end


covariate_files = '';
scansToInclude = [];
if (~strcmpi(covariates, 'none'))
    try
        covariate_files = mancovanInfo.userInput.feature_params.final.fnc_tc_covariates_userdata.filesList;
        scansToInclude = mancovanInfo.userInput.feature_params.final.fnc_tc_covariates_userdata.file_numbers;
        disp('Variance associated with the covariates will be removed from the timecourses ... ');
    catch
    end
end


load(mancovanInfo.userInput.ica_param_file);
subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess, 'flagTimePoints', sesInfo.flagTimePoints);
outDir = fileparts(mancovanInfo.userInput.ica_param_file);
sesInfo.outputDir = outDir;
fileIn = dir(fullfile(outDir, [sesInfo.calibrate_components_mat_file, '*.mat']));
filesToDelete = {};
if (length(fileIn) ~= sesInfo.numOfSub*sesInfo.numOfSess)
    disp('Uncompressing subject component files ...');
    filesToDelete = icatb_unZipSubjectMaps(sesInfo, subjectICAFiles);
end

if (~strcmpi(desCriteria, 'mancova'))
    good_sub_inds = (1:sesInfo.numOfSub);
else
    good_sub_inds = mancovanInfo.good_sub_inds;
end

fprintf('\n');
if ((step == 1) || (step == 2))
    %% Compute features
    %% Load subject components and average components across runs
    
    tp = min(sesInfo.diffTimePoints);
    
    comp_inds = mancovanInfo.comps;
    
    features = mancovanInfo.features;
    Vol = sesInfo.HInfo.V(1);
    mancovanInfo.HInfo = Vol;
    
    outputFiles = repmat(struct('feature_name', '', 'filesInfo', ''), 1, length(features));
    
    disp('Computing features ...');
    fprintf('\n');
    
    for nF = 1:length(features)
        
        cF = features{nF};
        
        outputFiles(nF).feature_name = cF;
        
        if (strcmpi(features{nF}, 'spatial maps'))
            %% Spatial maps
            
            if (~strcmpi(feature_params.final.sm_mask, 'default'))
                disp(['Loading mask ', feature_params.final.sm_mask_userdata, ' ...']);
                userMask = icatb_loadData(feature_params.final.sm_mask_userdata);
                userMask = find(abs(userMask) > eps);
                [dd, mask_rel_inds] = intersect(sesInfo.mask_ind, userMask);
            end
            
            countComp = 0;
            result_files = repmat({''}, 1, length(comp_inds));
            tmap_files = repmat({''}, 1, length(comp_inds));
            
            dirName = 'sm_stats';
            
            if (exist(fullfile(mancovanInfo.outputDir, dirName)) ~= 7)
                mkdir(mancovanInfo.outputDir, dirName);
            end
            
            for ncomps = comp_inds
                
                countComp = countComp + 1;
                disp(['Loading subject spatial maps of component ', num2str(ncomps), ' ...']);
                
                
                if (~strcmpi(desCriteria, 'paired t-test'))
                    SM = icatb_loadComp(sesInfo, ncomps, 'vars_to_load', 'ic', 'subjects', good_sub_inds, 'average_runs', 1, ...
                        'subject_ica_files', subjectICAFiles);
                    
                    % average runs
                    meanmap = mean(SM);
                    
                else
                    SM = icatb_loadComp(sesInfo, ncomps, 'vars_to_load', 'ic', 'subjects', good_sub_inds, 'average_runs', 0, ...
                        'subject_ica_files', subjectICAFiles);
                    SM = reshape(SM, sesInfo.numOfSub, sesInfo.numOfSess, length(sesInfo.mask_ind));
                    meanmap = squeeze(mean(mean(SM, 2)));
                    meanmap = meanmap(:)';
                end
                
                %% adjust magnitude
                stdterm = norm(meanmap) / sqrt(length(meanmap) - 1);
                meanmap = meanmap / stdterm;
                SM = SM / stdterm;
                
                offset = 0;
                
                if (strcmpi(feature_params.final.sm_center, 'yes'))
                    disp('Centering component spatial maps ...');
                    %% recenter
                    [meanmap, offset] = icatb_recenter_image(meanmap);
                end
                
                SM = SM - offset;
                if (~strcmpi(desCriteria, 'paired t-test'))
                    tmap = meanmap*sqrt(size(SM, 1)) ./ std(SM);
                else
                    tmap = meanmap*sqrt(size(SM, 1)) ./ std(squeeze(mean(SM, 2)));
                end
                
                sm_params.offset = offset;
                sm_params.std = stdterm;
                
                cutoff = 0;
                
                if (strcmpi(feature_params.final.sm_mask, 'default'))
                    %% determine cutoff
                    if (strcmpi(def_mask_stat, 't'))
                        disp(['Using tmap statistics (t > mean + ', num2str(def_mask_t_std), '*std) to compute default mask ...']);
                        [y, cutoff, fitPARAMS, P0, FH] = icatb_fitggmix(tmap, def_mask_t_std);
                        close(FH);
                        clear y;
                        
                        if (t_image_values == 1)
                            % positive and negative voxels
                            mask_rel_inds = (abs(tmap) >= abs(cutoff(2)));
                        elseif (t_image_values == 2)
                            % positive voxels
                            mask_rel_inds = (tmap >= abs(cutoff(2)));
                        else
                            % Negative voxels
                            mask_rel_inds = (tmap <= -abs(cutoff(2)));
                        end
                        
                    else
                        disp(['Using Z threshold of ', num2str(def_mask_z_thresh), ' on mean map to compute default mask ...']);
                        mask_rel_inds = (abs(meanmap./std(meanmap)) >= def_mask_z_thresh);
                    end
                end
                
                mask_ind = sesInfo.mask_ind(mask_rel_inds);
                
                if (isempty(mask_ind))
                    fprintf('\n');
                    errorMsg = ['No significant voxels found in component ', num2str(ncomps), '.\nTry using user specified mask or use default mask with z-statistic.'];
                    error('Error:Mancovan', errorMsg);
                    fprintf('\n');
                end
                
                if (~strcmpi(desCriteria, 'paired t-test'))
                    SM = SM(:, mask_rel_inds);
                    if (~strcmpi(desCriteria, 'mancova'))
                        datasetNo = [mancovanInfo.ttestOpts.t.val{:}];
                        SM = SM(datasetNo, :);
                    end
                else
                    SM = SM(:, :, mask_rel_inds);
                    datasetNo = [mancovanInfo.ttestOpts.t.val{:}];
                    SM = convertDataTo2D(SM, datasetNo, sesInfo);
                end
                
                tmpData = zeros(Vol.dim(1:3));
                tmpData(mask_ind) = tmap(mask_rel_inds);
                
                
                tmapName = fullfile(dirName, [mancovanInfo.prefix, '_tmap_', icatb_returnFileIndex(ncomps), '.img']);
                Vol.fname = fullfile(mancovanInfo.outputDir, tmapName);
                Vol.n(1) = 1;
                
                icatb_write_vol(Vol, tmpData);
                
                comp_est = 0;
                
                if (strcmpi(desCriteria, 'mancova'))
                    if (mancovanInfo.doEstimation)
                        %comp_est = icatb_estimate_dimension(SM, (tmpData ~= 0));
                        comp_est = order_selection(SM);
                    else
                        comp_est = mancovanInfo.numOfPCs(strmatch(cF, lower(mancovanInfo.features), 'exact'));
                    end
                end
                
                if (length(mask_ind) < comp_est)
                    fprintf('\n');
                    error('Error:Mancovan', ['No of voxels (', num2str(length(mask_ind)), ') in the mask is less than the number of desired PCs (', num2str(comp_est), ')', ...
                        '\nTry using user specified mask or use default mask with z-statistic.']);
                end
                
                resultsFile = fullfile(dirName, [mancovanInfo.prefix, '_results_sm_', icatb_returnFileIndex(ncomps), '.mat']);
                comp_number = ncomps;
                disp(['Saving file ', resultsFile, ' ...']);
                
                result_files{countComp} = resultsFile;
                tmap_files{countComp} = tmapName;
                
                icatb_save(fullfile(mancovanInfo.outputDir, resultsFile), 'cutoff', 'tmapName', 'comp_est', 'mask_ind', 'comp_number', 'sm_params', 'def_mask_t_std', ...
                    'SM');
                if (step == 1)
                    fprintf('\n');
                    disp(['Running Mancovan on ', cF, ' ...']);
                    disp('');
                    Stepwise_options = {};
                    if (strcmpi(desCriteria, 'mancova'))
                        Stepwise_options = {'reduced', mancovanInfo.modelInteractions.types{:}, 'SVD', 'FIXED', ['FIXED_' num2str(comp_est)]};
                    end
                    
                    
                    %                     [MULT, UNI] = run_model(mancovanInfo, SM, Stepwise_options, step_P);
                    %                     disp(['Saving file ', resultsFile, ' ...']);
                    %                     icatb_save(fullfile(mancovanInfo.outputDir, resultsFile), 'MULT', 'UNI', '-append');
                    %                     clear UNI MULT;
                    %                     fprintf('\n');
                    
                    
                    computeMancova(mancovanInfo, SM, Stepwise_options, step_P, fullfile(mancovanInfo.outputDir, resultsFile));
                    
                end
                
                clear mask_ind cutoff tmapName comp_est mask_ind;
                fprintf('\n');
                
            end
            outputFiles(nF).filesInfo.tmap_files = tmap_files;
        elseif (strcmpi(features{nF}, 'timecourses spectra'))
            %% Timecourses spectra
            
            dirName = 'spectra_stats';
            
            if (exist(fullfile(mancovanInfo.outputDir, dirName)) ~= 7)
                mkdir(mancovanInfo.outputDir, dirName);
            end
            
            spectra_params = struct('tapers', feature_params.final.spectra_tapers, 'Fs', feature_params.final.spectra_sampling_freq, 'fpass', feature_params.final.spectra_freq_band);
            
            countComp = 0;
            result_files = repmat({''}, 1, length(comp_inds));
            for ncomps = comp_inds
                
                countComp = countComp + 1;
                
                disp(['Loading subject timecourses of component ', num2str(ncomps), ' ...']);
                
                
                timecourses = icatb_loadComp(sesInfo, ncomps, 'vars_to_load', 'tc', 'subjects', good_sub_inds, 'truncate_tp', 1, ...
                    'subject_ica_files', subjectICAFiles, 'detrend_no', feature_params.final.spectra_detrend);
                
                timecourses = reshape(timecourses, size(timecourses, 1), sesInfo.numOfSess, size(timecourses, length(size(timecourses))));
                timecourses = timecourses(:, :, 1:min(sesInfo.diffTimePoints));
                
                disp('Doing multi-taper spectral estimation ...');
                
                for nSubjects = 1:size(timecourses, 1)
                    for nSessions = 1:sesInfo.numOfSess
                        tempTC = squeeze(timecourses(nSubjects, nSessions, :));
                        [temp_spectra, freq] = icatb_get_spectra(tempTC, mean(TR), spectra_params);
                        if ((nSubjects == 1) && (nSessions == 1))
                            spectra_tc = zeros(size(timecourses, 1), sesInfo.numOfSess, length(temp_spectra));
                        end
                        spectra_tc(nSubjects, nSessions, :) = temp_spectra;
                        clear tempTC;
                    end
                end
                
                clear timecourses;
                
                if (strcmpi(feature_params.final.spectra_normalize_subs, 'yes'))
                    disp('Using fractional amplitude ...');
                    spectra_tc = spectra_tc./repmat(sum(spectra_tc,3), [1, 1, size(spectra_tc, 3)]);
                end
                
                
                spectra_tc_all = spectra_tc;
                
                if (~strcmpi(desCriteria, 'paired t-test'))
                    spectra_tc = squeeze(mean(spectra_tc, 2));
                    if (~strcmpi(desCriteria, 'mancova'))
                        datasetNo = [mancovanInfo.ttestOpts.t.val{:}];
                        spectra_tc = spectra_tc(datasetNo, :);
                    end
                else
                    datasetNo = [mancovanInfo.ttestOpts.t.val{:}];
                    spectra_tc = convertDataTo2D(spectra_tc, datasetNo, sesInfo);
                end
                
                resultsFile = fullfile(dirName, [mancovanInfo.prefix, '_results_spectra_', icatb_returnFileIndex(ncomps), '.mat']);
                icatb_save(fullfile(mancovanInfo.outputDir, resultsFile), 'spectra_tc', 'spectra_tc_all');
                
                clear spectra_tc_all;
                
                if (strcmpi(feature_params.final.spectra_transform, 'yes'))
                    disp('Applying log transform to spectra ...');
                    %% normalizing transform
                    spectra_tc = log(spectra_tc);
                end
                
                comp_est = 0;
                
                if (strcmpi(desCriteria, 'mancova'))
                    if (mancovanInfo.doEstimation)
                        comp_est = order_selection(spectra_tc);
                    else
                        comp_est = mancovanInfo.numOfPCs(strmatch(cF, lower(mancovanInfo.features), 'exact'));
                    end
                end
                
                comp_number = ncomps;
                disp(['Saving file ', resultsFile, ' ...']);
                
                result_files{countComp} = resultsFile;
                
                icatb_save(fullfile(mancovanInfo.outputDir, resultsFile), 'comp_est', 'comp_number', 'freq', '-append');
                
                if (step == 1)
                    
                    fprintf('\n');
                    disp(['Running Mancovan on ', cF, ' ...']);
                    disp('');
                    %Stepwise_options = {'reduced', mancovanInfo.modelInteractions.types{:}, 'SVD', 'FIXED', ['FIXED_' num2str(comp_est)]};
                    Stepwise_options = {};
                    if (strcmpi(desCriteria, 'mancova'))
                        Stepwise_options = {'reduced', mancovanInfo.modelInteractions.types{:}, 'SVD', 'FIXED', ['FIXED_' num2str(comp_est)]};
                    end
                    
                    
                    computeMancova(mancovanInfo, spectra_tc, Stepwise_options, step_P, fullfile(mancovanInfo.outputDir, resultsFile));
                    
                end
                
                
                clear spectra_tc cutoff tmapName comp_est mask_ind freq;
                
                fprintf('\n');
                
            end
            
        else
            %% FNC correlations;
            
            dirName = 'fnc_stats';
            
            if (exist(fullfile(mancovanInfo.outputDir, dirName)) ~= 7)
                mkdir(mancovanInfo.outputDir, dirName);
            end
            
            tmpTR = TR;
            
            if (length(tmpTR) == 1)
                tmpTR = repmat(tmpTR(1), 1, sesInfo.numOfSub);
            else
                if (length(tmpTR) ~= sesInfo.numOfSub)
                    error(['TRs length must match the number of subjects (', num2str(sesInfo.numOfSub), ')']);
                end
            end
            
            tmpTR = tmpTR(good_sub_inds);
            
            disp('Loading subject timecourses of components ...');
            % %             timecourses = icatb_loadComp(sesInfo, comp_inds, 'vars_to_load', 'tc', 'subjects', good_sub_inds, 'truncate_tp', 1, ...
            % %                 'subject_ica_files', subjectICAFiles, 'detrend_no', feature_params.final.fnc_tc_detrend, 'covariates', ...
            % %                 covariate_files, 'scansToInclude', scansToInclude);
            % %
            % %             timecourses = reshape(timecourses, size(timecourses, 1), sesInfo.numOfSess, size(timecourses, length(size(timecourses)) - 1), length(comp_inds));
            % %             timecourses = timecourses(:, :, 1:min(sesInfo.diffTimePoints), :);
            
            if (strcmpi(feature_params.final.fnc_tc_despike, 'yes') && (filter_cutoff > 0))
                disp('Despiking and filtering timecourses ...');
            else
                if (strcmpi(feature_params.final.fnc_tc_despike, 'yes'))
                    disp('Despiking timecourses ...');
                elseif (filter_cutoff > 0)
                    disp('Filtering timecourses ...');
                end
            end
            
            numSubjects = length(good_sub_inds);
            for nSub = 1:numSubjects
                for nSess = 1:sesInfo.numOfSess
                    timecourses = icatb_loadComp(sesInfo, comp_inds, 'vars_to_load', 'tc', 'subjects', good_sub_inds(nSub), 'sessions', nSess, ...
                        'subject_ica_files', subjectICAFiles, 'detrend_no', feature_params.final.fnc_tc_detrend, 'covariates', ...
                        covariate_files, 'scansToInclude', scansToInclude);
                    
                    timecourses = squeeze(timecourses);
                    
                    for nComp = 1:length(comp_inds)
                        if (strcmpi(feature_params.final.fnc_tc_despike, 'yes'))
                            timecourses(:, nComp) = icatb_despike_tc(timecourses(:, nComp), tmpTR(nSub));
                        end
                        if (filter_cutoff > 0)
                            timecourses(:, nComp) = icatb_filt_data(timecourses(:, nComp), tmpTR(nSub), filter_cutoff);
                        end
                    end
                    c = icatb_corr(timecourses);
                    c(1:size(c, 1) + 1:end) = 0;
                    c = icatb_mat2vec(icatb_r_to_z(c));
                    if ((nSub == 1) && (nSess == 1))
                        fnc_corrs = zeros(numSubjects, sesInfo.numOfSess, numel(c));
                    end
                    fnc_corrs(nSub, nSess, :) = c(:)';
                end
            end
            
            fnc_corrs_all = fnc_corrs;
            
            if (~strcmpi(desCriteria, 'paired t-test'))
                fnc_corrs = squeeze(mean(fnc_corrs, 2));
                if (~strcmpi(desCriteria, 'mancova'))
                    datasetNo = [mancovanInfo.ttestOpts.t.val{:}];
                    fnc_corrs = fnc_corrs(datasetNo, :);
                end
            else
                datasetNo = [mancovanInfo.ttestOpts.t.val{:}];
                fnc_corrs = convertDataTo2D(fnc_corrs, datasetNo, sesInfo);
            end
            
            clear timecourses;
            
            comp_est = 0;
            
            if (strcmpi(desCriteria, 'mancova'))
                if (mancovanInfo.doEstimation)
                    comp_est = order_selection(fnc_corrs);
                else
                    comp_est = mancovanInfo.numOfPCs(strmatch(cF, lower(mancovanInfo.features), 'exact'));
                end
            end
            
            resultsFile = fullfile(dirName, [mancovanInfo.prefix, '_results_fnc.mat']);
            comp_number = comp_inds;
            disp(['Saving file ', resultsFile, ' ...']);
            
            result_files{1} = resultsFile;
            
            icatb_save(fullfile(mancovanInfo.outputDir, resultsFile), 'comp_est', 'comp_number', 'fnc_corrs', 'fnc_corrs_all');
            
            if (step == 1)
                
                fprintf('\n');
                disp(['Running Mancovan on ', cF, ' ...']);
                disp('');
                %Stepwise_options = {'reduced', mancovanInfo.modelInteractions.types{:}, 'SVD', 'FIXED', ['FIXED_' num2str(comp_est)]};
                Stepwise_options = {};
                if (strcmpi(desCriteria, 'mancova'))
                    Stepwise_options = {'reduced', mancovanInfo.modelInteractions.types{:}, 'SVD', 'FIXED', ['FIXED_' num2str(comp_est)]};
                end
                
                computeMancova(mancovanInfo, fnc_corrs, Stepwise_options, step_P, fullfile(mancovanInfo.outputDir, resultsFile));
                
            end
            
            
            clear fnc_corrs cutoff tmapName comp_est mask_ind freq;
            fprintf('\n');
            
        end
        
        outputFiles(nF).filesInfo.result_files = result_files;
        
        clear result_files;
        
    end
    
    mancovanInfo.outputFiles = outputFiles;
    %% Save Mancovan parameter file
    icatb_save(fullfile(mancovanInfo.outputDir, [mancovanInfo.prefix, '.mat']), 'mancovanInfo');
    
else
    
    if (~isfield(mancovanInfo, 'outputFiles'))
        error('Please run setup features prior to running mancova');
    end
    
    %% Run mancovan only
    outputFiles = mancovanInfo.outputFiles;
    for nF = 1:length(outputFiles)
        result_files = outputFiles(nF).filesInfo.result_files;
        cF = outputFiles(nF).feature_name;
        fprintf('\n');
        disp(['Running Mancovan on ', cF, ' ...']);
        disp('');
        for nR = 1:length(result_files)
            outFile = fullfile(mancovanInfo.outputDir, result_files{nR});
            if (strcmpi(cF, 'spatial maps'))
                %% Spatial maps
                
                load(outFile, 'comp_est', 'mask_ind', 'sm_params', 'comp_number');
                Stepwise_options = {};
                if (strcmpi(desCriteria, 'mancova'))
                    Stepwise_options = {'reduced', mancovanInfo.modelInteractions.types{:}, 'SVD', 'FIXED', ['FIXED_' num2str(comp_est)]};
                end
                
                fS = whos('-file', outFile);
                fNames = cellstr(char(fS.name));
                chkSM = ~isempty(strmatch('SM', fNames, 'exact'));
                
                if (chkSM)
                    load(outFile, 'SM');
                else
                    disp(['Loading subject spatial maps of component ', num2str(comp_number), ' ...']);
                    SM = icatb_loadComp(sesInfo, comp_number, 'vars_to_load', 'ic', 'subjects', good_sub_inds, 'average_runs', 1, ...
                        'subject_ica_files', subjectICAFiles);
                    % Use only the required indices
                    [dd, iA, iB] = intersect(sesInfo.mask_ind, mask_ind);
                    SM = SM(:, iA);
                    % apply spatial params
                    SM = SM/sm_params.std;
                    SM = SM - sm_params.offset;
                end
                
                
                computeMancova(mancovanInfo, SM, Stepwise_options, step_P, outFile);
                
                
            elseif (strcmpi(cF, 'timecourses spectra'))
                %% Spectra
                
                load(outFile, 'comp_est', 'spectra_tc');
                
                if (strcmpi(feature_params.final.spectra_transform, 'yes'))
                    spectra_tc = log(spectra_tc);
                end
                
                Stepwise_options = {};
                if (strcmpi(desCriteria, 'mancova'))
                    Stepwise_options = {'reduced', mancovanInfo.modelInteractions.types{:}, 'SVD', 'FIXED', ['FIXED_' num2str(comp_est)]};
                end
                
                computeMancova(mancovanInfo, spectra_tc, Stepwise_options, step_P, outFile);
                
            else
                %% Timecourse FNC correlations
                
                load(outFile, 'comp_est', 'fnc_corrs');
                Stepwise_options = {};
                if (strcmpi(desCriteria, 'mancova'))
                    Stepwise_options = {'reduced', mancovanInfo.modelInteractions.types{:}, 'SVD', 'FIXED', ['FIXED_' num2str(comp_est)]};
                end
                % Stepwise_options = {'reduced', mancovanInfo.modelInteractions.types{:}, 'SVD', 'FIXED', ['FIXED_' num2str(comp_est)]};
                
                computeMancova(mancovanInfo, fnc_corrs, Stepwise_options, step_P, outFile);
                
            end
            
            %             disp(['Saving file ', outFile, ' ...']);
            %             icatb_save(outFile, 'MULT', 'UNI', '-append');
            %             if (exist('time', 'var'))
            %                 icatb_save(outFile, 'time', '-append');
            %             end
            %             disp('Done');
            
            
            fprintf('\n');
            clear UNI MULT;
        end
    end
end


if (exist('filesToDelete', 'var') && ~isempty(filesToDelete))
    icatb_cleanupFiles(filesToDelete, outDir);
end

totalTime = toc;

disp('Analysis Complete');

disp(['Total time taken to complete the analysis is ', num2str(totalTime/60), ' minutes']);

diary('off');

function [comp_est, mdl, aic, kic] = order_selection(data)

disp('Estimating dimension ...');

%% Remove mean
data = detrend(data, 0);

%% Arrange data based on correlation threshold
[V1, D1] = icatb_svd(data);

lam = diag(D1);

lam = sort(lam);

lam = lam(end:-1:1);

lam = lam(:)';

N = (size(data, 2));

%% Make eigen spectrum adjustment
tol = max(size(lam)) * eps(max(lam));

if (lam(end) < tol)
    lam(end) = [];
end

%% Correction on the ill-conditioned results (when tdim is large, some
% least significant eigenvalues become small negative numbers)
lam(real(lam) <= tol) = tol;

p = length(lam);
aic = zeros(1, p - 1);
kic = zeros(1, p - 1);
mdl = zeros(1, p - 1);
for k = 1:p-1
    LH = log(prod(lam(k+1:end).^(1/(p-k)) )/mean(lam(k+1:end)));
    mlh = 0.5*N*(p-k)*LH;
    df = 1 + 0.5*k*(2*p-k+1);
    aic(k) =  -2*mlh + 2*df;
    kic(k) =  -2*mlh + 3*df;
    mdl(k) =  -mlh + 0.5*df*log(N);
end

% Find the first local minimum of each ITC
itc = zeros(3, length(mdl));
itc(1,:) = aic;
itc(2,:) = kic;
itc(3,:) = mdl;

%% Use only mdl
dlap = squeeze(itc(end, 2:end)-itc(end, 1:end-1));
a = find(dlap > 0);
if isempty(a)
    comp_est = length(squeeze(itc(end, :)));
else
    comp_est = a(1);
end

disp(['Estimated components is found to be ', num2str(comp_est)]);


if ((comp_est < 2) || all(diff(mdl) < 0))
    comp_est = min([6, min(size(data))]);
    warning('Mancovan:DimensionalityEstimation', 'Estimated components is 1 or MDL function is monotonically decreasing. Using %d instead', comp_est);
end

disp('Done');


function mancovanInfo = mancovan_err_chk(mancovanInfo)
%% Do error check
%

if (~isfield(mancovanInfo.userInput, 'features'))
    error('Please run setup features first');
end

if (isempty(mancovanInfo.userInput.features))
    error('Select features in setup analysis');
end


if (strcmpi(mancovanInfo.designCriteria, 'mancova'))
    % Check covariates
    for nC = 1:length(mancovanInfo.userInput.cov)
        name = mancovanInfo.userInput.cov(nC).name;
        val = mancovanInfo.userInput.cov(nC).value;
        if (isempty(name) && isempty(val))
            error('Please specify name and value for each covariate in setup analysis');
        end
    end
end

% Check components
for nC = 1:length(mancovanInfo.userInput.comp)
    name = mancovanInfo.userInput.comp(nC).name;
    val = mancovanInfo.userInput.comp(nC).value;
    if (isempty(name) && isempty(val))
        error('Please specify name and value for each component network name in setup analysis');
    end
end

if (strcmpi(mancovanInfo.designCriteria, 'mancova'))
    if (~isfield(mancovanInfo.userInput, 'modelInteractions'))
        mancovanInfo = icatb_mancovan_interactions(mancovanInfo);
    end
end

drawnow;


function [MULT, UNI] = run_model(mancovanInfo, data, Stepwise_options, step_P)
%% Run model
%

if (strcmpi(mancovanInfo.designCriteria, 'mancova'))
    
    X = mancovanInfo.X;
    start_terms = mancovanInfo.regressors;
    terms = mancovanInfo.terms;
    term_names = start_terms;
    
    %% Create labels for the columns
    cov_names = cellstr(char(mancovanInfo.cov.name));
    cov_labels = [mancovanInfo.cov.labels];
    
    check = 0;
    while ~check
        try
            %% mulivariate TEST
            [ T_m, p_m, stats_m ] = mStepwise(data, X, terms, step_P, Stepwise_options);
            check = 1;
            
        catch
            err = lasterror;
            if (strcmpi(err.identifier, 'MATLAB:betainc:XOutOfRange'))
                num = str2num(strrep(Stepwise_options{end}, 'FIXED_', ''));
                fprintf('\n');
                warning(['mStepwise did not work with ', num2str(num), ' components']);
                disp(['Rerunning mStepwise with ', num2str(num - 1), ' components']);
                fprintf('\n');
                Stepwise_options{end} = ['FIXED_', num2str(num - 1)];
            else
                rethrow(err);
            end
        end
    end
    
    
    %% find the significant terms in the model
    [sig_terms, I, J]= mUnique(stats_m.Terms);
    
    for ii = 1:length(sig_terms)
        temp_sig{ii} = num2str(sig_terms{ii});
    end
    
    for ii = 1:length(terms)
        temp_term{ii} = num2str(terms{ii});
    end
    
    X_reduced_ind = [];
    for ii = 1:length(temp_sig)
        x_ind = find(strcmp(temp_sig{ii}, temp_term));
        X_reduced_ind = [X_reduced_ind x_ind];
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_reduced_ind = sort(X_reduced_ind);
    X_reduced = X(:, X_reduced_ind);
    
    red_term_names = term_names(X_reduced_ind(2:end)-1); % a little awkward since no constant term in the names
    test_names = red_term_names(I(2:end)-1);% a little awkward since no constant term in the names
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if length(test_names) > 0
        %% univariate TEST
        %[ t, p, stats ] = mT(Y, X, terms, term, options)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ii = 2:length(sig_terms)  % no test for constant
            fprintf('Working on term %d of %d\n', ii-1, length(sig_terms)-1)
            [ t_u{ii-1}, p_u{ii-1}, stats_u{ii-1}] = mT(data, X_reduced, stats_m.Terms, sig_terms{ii}, {'verbose'});
            stats_u{ii-1} = get_contrast_label(stats_u{ii-1}, test_names{ii-1}, cov_names, cov_labels);
        end
    else
        t_u = [];
        stats_u = [];
        p_u = [];
    end
    
    
    MULT.X = X;
    MULT.start_terms = unique(start_terms);
    MULT.p_cutoff = step_P;
    MULT.final_terms = test_names;
    MULT.t = T_m;
    MULT.p = p_m;
    MULT.stats = stats_m;
    
    UNI.tests = test_names;
    UNI.t = t_u;
    UNI.p = p_u;
    UNI.stats = stats_u;
    
else
    
    ttestNames = {mancovanInfo.designCriteria};
    if (strcmpi(mancovanInfo.designCriteria, 'one sample t-test'))
        [ t_u{1}, p_u{1}, stats_u{1}] = mT(data, mancovanInfo.X, [], 0, {'verbose'});
        stats_u{1} = get_contrast_label(stats_u{1}, ttestNames, ttestNames, mancovanInfo.ttestOpts.t.name);
    elseif (strcmpi(mancovanInfo.designCriteria, 'two sample t-test'))
        [t_u{1}, p_u{1}, stats_u{1}] = mT(data, mancovanInfo.X, [], 1, {'verbose'});
        stats_u{1} = get_contrast_label(stats_u{1}, ttestNames, ttestNames, mancovanInfo.ttestOpts.t.name);
    else
        N1 = length(mancovanInfo.ttestOpts.t.val{1});
        [t_u{1}, p_u{1}, stats_u{1}] = mT(data(1:N1, :) - data(N1 + 1:end, :), mancovanInfo.X, [], 0, {'verbose'});
        stats_u{1} = get_contrast_label(stats_u{1}, ttestNames, ttestNames, mancovanInfo.ttestOpts.t.name);
    end
    
    UNI.tests = {mancovanInfo.designCriteria};
    UNI.t = t_u;
    UNI.p = p_u;
    UNI.stats = stats_u;
    MULT = [];
    
end

function s = get_contrast_label(s, tname, allnames, alllabels)

if length(s.Term) == 1 %main effect
    varIND = find(strcmp(tname, allnames));
    labels = alllabels{varIND};
    if (~iscell(labels))
        labels = {labels};
    end
    for jj = 1:size(s.Levels,1)
        if s.Levels(jj,2) == 0 && length(labels) == 1
            s.Contrast{jj} = labels{1};
        else
            s.Contrast{jj} = ['(' labels{s.Levels(jj,1)+1} ') - (' labels{s.Levels(jj,2)+1} ')'];
            
        end
    end
    
else %interaction
    term1_end = strfind(tname, '_X_')-1;
    term2_start = term1_end + 4;
    clear varIND
    varIND(1) = find(strcmp(tname(1:term1_end), allnames));
    varIND(2) = find(strcmp(tname(term2_start:end), allnames));
    labels_1 = alllabels{varIND(1)};
    labels_2 = alllabels{varIND(2)};
    
    if (~iscell(labels_1))
        labels_1 = {labels_1};
    end
    
    if (~iscell(labels_2))
        labels_2 = {labels_2};
    end
    
    
    for jj = 1:size(s.Levels,1)
        if length(labels_1)*length(labels_2) == 1 % both are continuous variables
            s.Contrast{jj} = ['(' labels_1{1} ') X (' labels_2{1} ')'];
        elseif length(labels_1) == 1 || length(labels_2) == 1 %one continuous, one categorical
            
            if length(labels_1) > length(labels_2)
                temp = labels_1;
                labels_1 = labels_2;
                labels_2 = temp;
            end
            
            s.Contrast{jj} = ['(' labels_1{1} ') X [(' labels_2{s.Levels(jj,1)+1} ') - (' labels_2{s.Levels(jj,2)+1} ')]'];
        elseif length(labels_1)*length(labels_2) == 4 %categorical, two by two
            s.Contrast{jj} = ['[(' labels_1{s.Levels(jj,1)+1} ') - (' labels_1{s.Levels(jj,2)+1} ')] X [(' ...
                labels_2{s.Levels(jj,1)+1} ') - (' labels_2{s.Levels(jj,2)+1} ')]'];
        elseif length(labels_1)*length(labels_2) == 6 %categorical, two by three
            if length(labels_1) > length(labels_2)
                temp = labels_1;
                labels_1 = labels_2;
                labels_2 = temp;
            end
            s.Contrast{jj} = ['[(' labels_1{2} ') - (' labels_1{1} ')] X [(' labels_2{s.Levels(jj,1)+1} ') - (' labels_2{s.Levels(jj,2)+1} ')]'];
            
        else
            s.Contrast{jj} = tname;
            %             combinations = [];
            %             for j = 1 : length(labels_1)
            %                 for k = 1 :length(labels_2)
            %                     combinations(end + 1, :) = [ j k ];
            %                 end
            %             end
        end
    end
    
end


function  SM_new = convertDataTo2D(SM, datasetNo, sesInfo)

lastDim = size(SM, 3);

SM_new = zeros(length(datasetNo), lastDim);
for nDataSet = 1:length(datasetNo)
    subNum = ceil(datasetNo(nDataSet)/sesInfo.numOfSess);
    sesNum = mod(datasetNo(nDataSet) - 1, sesInfo.numOfSess) + 1;
    SM_new(nDataSet, :) = SM(subNum, sesNum, :);
end


function computeMancova(mancovanInfo, SM, Stepwise_options, step_P, resultsFile)


if (isfield(mancovanInfo, 'time'))
    
    disp('Computing mancova for timepoint 1 - timepoint 2 ...');
    [MULT, UNI] = run_model(mancovanInfo, SM(mancovanInfo.time.subjects{1}, :) - SM(mancovanInfo.time.subjects{2}, :), Stepwise_options, step_P);
    disp(['Saving file ', resultsFile, ' ...']);
    icatb_save(resultsFile, 'MULT', 'UNI', '-append');
    clear UNI MULT;
    fprintf('\n');
    
    % Time 1
    disp('Computing mancova for timepoint 1 ...');
    [time.MULT{1}, time.UNI{1}] = run_model(mancovanInfo, SM(mancovanInfo.time.subjects{1}, :), Stepwise_options, step_P);
    fprintf('\n');
    
    % Time 2
    disp('Computing mancova for timepoint 2 ...');
    [time.MULT{2}, time.UNI{2}] = run_model(mancovanInfo, SM(mancovanInfo.time.subjects{2}, :), Stepwise_options, step_P);
    icatb_save(resultsFile, 'time', '-append');
    clear UNI MULT;
    fprintf('\n');
    
else
    
    [MULT, UNI] = run_model(mancovanInfo, SM, Stepwise_options, step_P);
    disp(['Saving file ', resultsFile, ' ...']);
    icatb_save(resultsFile, 'MULT', 'UNI', '-append');
    clear UNI MULT;
    fprintf('\n');
    
end
