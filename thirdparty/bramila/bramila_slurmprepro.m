function bramila_slurmprepro(cfg)
% BRAMILA_SLURMPREPRO - Creates cfg, job file and log file, and submits job
% file to slurm on triton cluster.
%   - Usage:
%   bramila_slurmprepro(cfg)
%   - Input:
%   cfg is a struct containing all the parameters from
%   BRAMILA_PREPRO_DEMO.m, among which, obligatory are:
%       cfg.subj = path to folder containing epi.nii and bet.nii
%       cfg.bramilapath = path to BRAMILA toolbox
%	Last edit: DS 2014-02-25
%%
% job filename
if ~exist(cfg.outputfolder, 'dir')
	mkdir(cfg.outputfolder)
end
filename = fullfile(cfg.outputfolder,'preprocessing_job');
% log filename
logfile = fullfile(cfg.outputfolder,'preprocessing_logfile');
% save cfg separately so that job has access to it
save(fullfile(cfg.outputfolder,'cfg.mat'),'cfg','-v7.3')
%load the modules
dlmwrite(filename, '#!/bin/sh', '');
dlmwrite(filename, '#SBATCH -p batch','-append','delimiter','');
dlmwrite(filename, '#SBATCH -t 04:00:00','-append','delimiter','');
dlmwrite(filename, '#SBATCH --qos=normal','-append','delimiter','');
dlmwrite(filename, ['#SBATCH -o "' logfile '"'],'-append','delimiter','');
% Adjust memory requirement
memreq_line = check_mem_per_core(sprintf('%s/epi.nii',cfg.inputfolder));
dlmwrite(filename, memreq_line,'-append','delimiter','');
dlmwrite(filename, 'module load matlab','-append','delimiter','');
dlmwrite(filename, 'module load fsl','-append','delimiter','');
dlmwrite(filename, 'source $FSLDIR/etc/fslconf/fsl.sh','-append','delimiter','');

dlmwrite(filename, ['cd ' cfg.bramilapath],'-append','delimiter','');

%command
dlmwrite(filename,sprintf('matlab -nosplash -nodisplay -nodesktop -r "bramila_preprocessor(''%s/cfg.mat'');exit;"',cfg.outputfolder),'-append','delimiter','');
unix(['sbatch ' filename]);
