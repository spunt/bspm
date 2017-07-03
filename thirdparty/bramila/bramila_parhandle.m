function cfg = bramila_parhandle(varargin)
% BRAMILA_PARHANDLE - Checks for parallelization possibilities in current
% preprocessing configuration and reslices the cfg subjectwise.
% It uses the function testGrid from ISC toolbox.
%   - Usage:
%   cfg = bramila_parhandle(cfg,subjects)
%   cfg = bramila_parhandle(cfg,subjects,subjects_out)
%   - Input:
%   cfg is a struct with following parameters
%       subjects = cell array with paths to folders containing epi.nii and
%       bet.nii
%       cfg = other preprocessing settings
%   - Output:
%       cfg = cell array containing resliced cfg, where cfg{1}.inputfolder is
%       filled with individual subject folder
%	Last edit: EG 2014-07-14

cfg = varargin{1};
subjects = varargin{2};

% fix the backward compatibility regarding subjects_out=subjects
same_in_out=0;
if(nargin==3)
	subjects_out = varargin{3};
	%check that subj in and out have same length AND force them to be different
	NI=length(subjects);
	NO=length(subjects_out);
	if(NI~=NO)
		error('There is a mismatch between the number of input/output folders');
	end
	for n=1:NI
		if(strcmp(subjects{n},subjects_out{n}))
			same_in_out=1;
		end
	end
else	
	subjects_out = subjects;
	same_in_out=1;
end

if(same_in_out==1)
	error('Input and output folder cannot be the same anymore. Please specify a different output folder.')
end

% You need the ISC toolbox to run the below
local=0;
if(isfield(cfg,'local'))
	local=cfg.local;
end

if(local == 0)
	gridtype = testGrid; % used to decide, whether to use matlabpool (dione, for example) or slurm (triton)
else
	gridtype = 0;
end

if(gridtype == 0)
	[aaa bbb]=system('uname -n');
	disp(['We are on machine ' bbb(1:end-1) ', preprocessing will run locally']);	
end

Nsub=length(subjects);

% adding an exception for local machine jouni. If a job is submitted from jouni, it will not run on the grid.
%if(strcmp(bbb(1:5),'jouni'))
%	disp('We are on jouni, preprocessing will run locally');
%	gridtype=0;
%end


if Nsub > 1
    if(strcmp(gridtype,'slurm'))
        disp('You have slurm, running things on cluster!')
    else
        poolobj = gcp('nocreate'); % Check if pool is open
        if isempty(poolobj)
            disp('You have no slurm, but you can use matlabpool')
            parpool
        else
            disp('matlabpool already open')
        end
    end
end

% Build individual cfg for each subject to run them in loop
tempcfg = cfg; 
tempcfg.gridtype = gridtype;
cfg = [];
for i = 1:Nsub
    cfg{i} = tempcfg;
    cfg{i}.inputfolder = subjects{i};
	cfg{i}.outputfolder = subjects_out{i};
end
