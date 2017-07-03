
function bramila_run_slurm_scripts(blockID,fullpath)
% BRAMILA_RUN_SLURM_SCRIPTS - simple function that, for each block of your parallelized code, generates the triton job file and submits it to the queue
%   - Usage:
%   	bramila_run_slurm_scripts(blockID,fullpath)
%			blockID = used as unique identifier for the specific block
%			fullpath = a path with output folder where the job files, log files and results will be stored
%	- Notes:
%		There is lots of space for improvement by adding more parameters (e.g. subfolders for storage, common file storage etc etc)

filename = fullfile([fullpath '/' num2str(blockID) '.job']);
logfile=fullfile([fullpath '/' num2str(blockID) '.log']);
dlmwrite(filename, '#!/bin/sh', '');
dlmwrite(filename, '#SBATCH -p batch','-append','delimiter','');
dlmwrite(filename, '#SBATCH -t 04:00:00','-append','delimiter','');
dlmwrite(filename, '#SBATCH --qos=normal','-append','delimiter','');
dlmwrite(filename, ['#SBATCH -o "' logfile '"'],'-append','delimiter','');
dlmwrite(filename, '#SBATCH --mem-per-cpu=30000','-append','delimiter','');
dlmwrite(filename, 'module load matlab','-append','delimiter','');
dlmwrite(filename, ['cd ' fullpath],'-append','delimiter','');
dlmwrite(filename,sprintf('matlab -nosplash -nodisplay -nodesktop -r "bramila_runtheblock(%d);exit"',blockID),'-append','delimiter','');
% and now we run it
unix(['sbatch ' filename]);
