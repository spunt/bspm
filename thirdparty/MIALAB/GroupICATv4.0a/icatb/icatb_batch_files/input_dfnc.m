%% Use full path for directories and files wherever needed. After entering parameters, use command icatb_dfnc_batch(input_file);

%% Output directory to place results 
outputDir = '/na1/batch_results';

%% ICA parameter file 
ica_param_file = '/na1/rest_hcp_ica_parameter_info.mat';


%% Cell array of dimensions number of network names by 2. Don't duplicate components in different
% network names
comp_network_names = {'BG', 21;                    % Basal ganglia 21st component
                      'AUD', 17;                   % Auditory 17th component
                      'SM', [7 23 24 29 38 56];    % Sensorimotor comps
                      'VIS', [39 46 48 59 64 67];  % Visual comps
                      'DMN', [25 50 53 68];        % DMN comps
                      'ATTN', [34 52 55 60 71 72]; % ATTN Comps
                      'FRONT', [20 42 47 49]};     % Frontal comps


%% TR of the experiment
TR = 2;


%% dFNC params

% 1. tc_detrend - Detrend number used to remove the trends in timecourses.
% Options are 0, 1, 2 and 3.
% 2. tc_despike - Remove spikes from the timecourses. Options are 'yes' and
% 'no'.
% 3. tc_filter - High frequency cutoff.

% 4. a. tc_covariates.filesList - Include character array or cell array of
% covariates to be regressed out from the timecourses. Cell array is of
% length number of subjects * sessions by 1. The order of file names will be
% first subject sessions followed by second subject sessions and so on.
%    b.  tc_covariates.file_numbers - Enter scan numbers to include. Leave
%    empty if you want to select all.
%


% 5. Regularisation method - Options are 'none' and 'L1'. 
% 6. wsize - Window size (scans) 
% 7. window_alpha - Gaussian Window alpha value.
% 8. num_repetitions - No. of repetitions (L1 regularisation).

dfnc_params.tc_detrend = 3;
dfnc_params.tc_despike = 'yes';
dfnc_params.tc_filter = 0.15;

dfnc_params.tc_covariates.filesList = [];
dfnc_params.tc_covariates.file_numbers = [];

dfnc_params.method = 'none';
dfnc_params.wsize = 30;
dfnc_params.window_alpha = 3;
dfnc_params.num_repetitions = 10;





