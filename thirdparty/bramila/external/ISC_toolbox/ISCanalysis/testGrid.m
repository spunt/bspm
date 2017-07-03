function grid_type = testGrid
% tests if Slurm, or SGE grid is available
% returns 'sge or 'slurm' in string variable
% if no grid system is detected returns empty string (grid disabled)
% if windows system is detected returns empty string (grid disabled)
%
% 26.6.2013
% Juha Pajula, Tampere University of Technology 
% juha.pajula@tut.fi 

% by default grid computing is disabled. 
grid_type='';

%check the OS type (grids works here only on linux):
OSsys = computer;
if(strcmp(OSsys,'PCWIN') || strcmp(OSsys,'PCWIN64'))
    disp('Windows system detected, disabling grid computing over SGE/Slurm')
    return;
else
    [status1,~]=unix('command -v qsub'); %for SGE status == 0, if qsub is found
    [status2,~]=unix('command -v srun'); %for Slurm

    if ~status1
        disp('SGE system detected');
        grid_type = 'sge';
    end
    if ~status2
        disp('Slurm system detected');
        grid_type = 'slurm';
    end
end    
%NOTE: This is possible to expand over other grid engines if needed:
% [statusX,~]=unix('command -v execute_command'); %for system X
% if ~statusX; grid_type = 'systemX';end
%
% possible extension might be for example Condor: if it is enabled the
% parsering must be added to gridParser function too, It uses the output
% from here to determine how the grid command is parsed (now there is sge 
% and slurm supported).
end