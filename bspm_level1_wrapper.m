%       images:   functional images
%
%       general_info:   a structure with the following fields
%           analysis:   path for the analysis
%           TR:   repetition time (in seconds)
%           hpf:   high-pass filter cutoff to use (in seconds)
%           autocorrelation:   0=None, 1=AR(1), 2=Weighted Least Squares (WLS)
%           nuisance_file:   txt file with nuisance regressors (leave empty for none)
%           brainmask:   brainmask to use (leave empty for none)
%
%       runs:   a structure with the following fields
%           conditions:
%               name:   string naming the condition
%               onsets:   onsets
%               durations:   durations
%               parameters:   for building parametric modulators (leave empty for none)
%                   name:   string naming the paramter
%                   values:   parameter values (assumed to be orthogonalized)
%
%       contrasts:  a structure with the following fields
%           type:   'T' or 'F'
%           name:   string naming the contrast
%           weights:    vector of contrast weights
%
% -----------------------------------------------------------

% --------------------------------- Copyright (C) 2014 ---------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

studydir = '/Users/bobspunt/Documents/fmri/lois';
subid = 'RA0584*';
runid = '*TOM*';
imid = 'su*nii';
behavid = 'tom*mat';
nuisanceid = 'bad*txt';

subdir = files([studydir filesep subid]);
images = files([studydir filesep subid filesep 'raw' filesep runid filesep imid]);
behav = files([studydir filesep subid filesep 'behav' filesep behavid]);
nuisance = files([studydir filesep subid filesep 'raw' filesep runid filesep nuisanceid]);
load(behav{1});

% TOM - SEEKER COLUMN KEY
% 1 - trial #
% 2 - condition (1=Belief, 2=Photo)
% 3 - intended trial onset
% 4 - intended question onset
% 5 - stimulus index
% 6 - actual story onset
% 7 - actual question onset
% 8 - actual response
% 9 - RT to question onset
% 10 - actual block duration

general_info.analysis = [subdir{1} filesep 'analysis' filesep 'tom_test_smooth'];
general_info.TR = 2.5;
general_info.hpf = 128;
general_info.autocorrelation = 2;
general_info.nuisance_file = nuisance{1};
general_info.brainmask = '';

names = {'Belief' 'Photo'};

for r = 1
    for c = 1:2
        runs(r).conditions(c).name = names{c};
        runs(r).conditions(c).onsets = Seeker(Seeker(:,2)==c,6);
        runs(r).conditions(c).durations = Seeker(Seeker(:,2)==c,10);
    end
end

contrasts(1).type = 'T';
contrasts(1).name = 'Belief-Photo';
contrasts(1).weights = [1 -1];

% RUN
bspm_level1(images, general_info, runs, contrasts)

        
  
 
 
 
 
