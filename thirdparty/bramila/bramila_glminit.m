% Submission script for triton cluster
% See wiki for detailed instructions
clear all
close all
addpath(genpath('/triton/becs/scratch/braindata/shared/toolboxes/bramila/bramila'));
%% Parameters will be stored in single struct to simplify function call
% The data should be structured in a way /dataroot/subject/session/inputfile
dataroot = '/triton/becs/scratch/braindata/shared/2pnhyperEMO'; % where subject folders are located
outdirname = '1stlvl'; % Where would the output go
cfg.TR = 1.7; % TR, get from protocol
cfg.AR = 'none'; % or 'AR(1)' for autoregression correction
cfg.sesses = {'run1','run2','run3','run4','run5'}; % sessions, or runs (see preprocessing guidelines)
cfg.hpf = 512; % high-pass filter, you can play with cutoffcalc in fsl to find a good estimate
cfg.incmoves = 0; % include or not include motion params?
cfg.inputfile = 'epi_preprocessed.nii'; % how the input epi file will be called
cfg.motionfile = 'epi_MCF.nii.par'; % how the input motion parameters file is called
cfg.cnames = {'story'}; % names of conditions or tasks that you model
cfg.regressor_unit = 'secs'; % 'scans' or 'secs', be VERY CAREFUL ABOUT IT
%% Parametric modulations
cfg.pmod = 0; % set to 0 if no parametric modulations, set to 1 if you have some
cfg.pmodname = {'valence','arousal'}; % leave empty if no pmods
% Subject folders
subjects = {
'Emo_listening_4';
};
% Example of how the input files would look like, examine it in case something is wrong
fprintf('Example of input file: %s/%s/%s/%s',dataroot,subjects{1},cfg.sesses{1},cfg.inputfile);
disp('If not fine, check variables: dataroot, subjects, cfg.sesses and cfg.inputfile');
%% Predefine contrasts
% If you want to have your own contrasts, define them here
% Contrasts should be of size(cname), don't include motion and constant
% By default, if not predefined contrasts, you will get main effect (1 0) and one vs other (1 -1; -1 1) contrasts
% cfg.predefined_contrast = []; % if no predefined contrasts
% Example of how to define contrasts
cfg.predefined_contrast{1} = [1 0 0];
cfg.predefined_contrast_name{1} = 'story';
cfg.predefined_contrast{1} = [0 1 0];
cfg.predefined_contrast_name{1} = 'arousal';
cfg.predefined_contrast{1} = [0 0 1];
cfg.predefined_contrast_name{1} = 'valence';
%% Main
% Which function to run on cluster
runfile = 'bramila_glm1st'; % how the function that will be executed is called
for s = 1:length(subjects)
    subject = subjects{s};
    %% Build regressor, varies from person to person, but logic is following
    % cfg.onset should be cell array of size(sesses).
    % each cell contains onsets for each condition in separate cell
    % i.e. if you have 3 tasks or conditions, each cell of cfg.onset contains
    % three cells where onsets of those conditions are specified.
    % Same logic applies to cfg.durations
    % Run this bit to see it as an example
    load(sprintf('%s/%s/Ratings/regressor.mat',dataroot,subject));
    % 1 - onset
    % 2 - duration
    % 3 - boxcar model for speech
    % 4 - arousal
    % 5 - valence
    % 6 - story number (gist)    
    for r = 1:length(cfg.sesses) % for each run        
        for c = 1:length(cfg.cnames)
            cfg.onsets{r}{c} = regressor{r}(:,1);
            cfg.durations{r}{c} = regressor{r}(:,2);
            if cfg.pmod == 1
                cfg.pmodparam{r}{c}{1} = regressor{r}(:,4);
                cfg.pmodparam{r}{c}{2} = regressor{r}(:,5);
            end
        end
    end
    % save the cfgfile
    save(sprintf('%s/%s/%s.mat',dataroot,subject,outdirname),'cfg');
    % job filename
    filename = fullfile(dataroot,subject,[outdirname 'job']);
    % log filename
    logfile = fullfile(dataroot,subject,[outdirname 'logfile']);
    %load the modules
    dlmwrite(filename, '#!/bin/sh', '');
    % Handle the long jobs
    dlmwrite(filename, '#SBATCH -p batch','-append','delimiter','');
    dlmwrite(filename, '#SBATCH -t 04:00:00','-append','delimiter','');
    dlmwrite(filename, '#SBATCH --qos=normal','-append','delimiter','');
    dlmwrite(filename, ['#SBATCH -o "' logfile '"'],'-append','delimiter','');
    dlmwrite(filename, '#SBATCH --mem-per-cpu=2000','-append','delimiter','');
    dlmwrite(filename, 'cd /triton/becs/scratch/braindata/shared/toolboxes/bramila/bramila','-append','delimiter','');
    %command
    dlmwrite(filename, sprintf('matlab -nosplash -nodisplay -nojvm -nodesktop -r "%s(''%s'',''%s'',''%s'')"',runfile,subject,dataroot,outdirname),'-append','delimiter','');
    unix(['sbatch ' filename]);
end


