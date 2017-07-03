function [status,outpt]=gridParser(command, Params, InPuts, grid_type,note)
%Function which parses the SGE/Slurm commandline script for ISCToolbox
%Usage: [status,outpt]=gridParser(command, Params, InPuts, grid_type)
% status = 1 if succeed, 0 if failure
% outpt = output message from grid engine
%
% command = string containing the name of usd command/function
% Params = Params struct of ISCtoolbox
% InPuts = cell of all variables of the 'command' (exept Params)
% grid_type = grid_type from testGrid function ('sge' or 'slurm')
%
% Example: 
% [~,outpt]=gridParser('PearsonFilon',Params,{nrSession,1:freqComps},grid_type);
%
% Created 26.6.2013
% Juha Pajula, Tampere University of Technology 
% juha.pajula@tut.fi 
     
status = false;
outpt = 'No grid engine commands submitted';

home_folder = pwd;
%define the script destination
script_folder=[Params.PublicParams.dataDestination,'scripts'];
%move to script folder (the original path is in home_folder variable)
%this is done to avoid windows/unix slash problem
cd(script_folder);

%find the toolboxpath
toolboxpath=which('ISCanalysis'); toolboxpath=toolboxpath(1:end-25);
%find the Params definition path
Params_path = [Params.PublicParams.dataDestination, Params.PublicParams.dataDescription];

%to ensure that ISCtoolbox is on the matlab path also in the gridengine
%nodes:
%load the Params
init_com = ['addpath(genpath(''' toolboxpath ''')); load(''' Params_path ''');']; 

%create the matlab command for ISCtoolbox phase (must end "exit"):
ISC_com = [command,'(Params'];
file_part = [];
for k = 1:length(InPuts)
    ISC_com = [ISC_com, ', ', '[' num2str(InPuts{k}) ']'];
    file_part = [file_part, '_' num2str(sum(InPuts{k})) ];
end
ISC_com = [ISC_com ,'); exit'];

%define scriptname for qsub/sbatch:
%file_name = ['ISCgrid_',command,'_',num2str(cputime),'.sh'];
file_name = [note,'_',command,file_part,'.sh'];
time_com = ['/usr/bin/time -f "\n%E elapsed,\n%U user(s),\n%S system(s),\n%M memory (kb)\n" -o ' file_name '.time.txt -a'];

%if the grid is SGE
if(strcmp(grid_type,'sge'))
    %write the scriptfile and run it, XYZ needs fo qsub:
    dlmwrite(file_name,[time_com ' matlab -nosplash -nodisplay -nojvm -nodesktop -r "' init_com ISC_com '"'], '')
    dlmwrite(file_name,'exit','-append','delimiter','');
    [status,outpt]=unix(['qsub -cwd -V -l arch=lx24-amd64,h_rt=120:0:0,mf=8G ' file_name ]); %status = 1 if succesful submission
end

%if the grid is Slurm
if(strcmp(grid_type,'slurm'))
    %load the modules
    dlmwrite(file_name, '#!/bin/sh', '');
    dlmwrite(file_name, 'module load matlab','-append','delimiter','');
    dlmwrite(file_name,[time_com ' matlab -nosplash -nodisplay -nojvm -nodesktop -r "' init_com ISC_com '"'],'-append','delimiter','');
    dlmwrite(file_name,'exit','-append','delimiter','');
    [status,outpt]=unix(['sbatch -J "' command '" --partition=sgn --mem=8192 --time=2-0 --error="' note '.e%j" --output="' note '.o%j" ' file_name ]); %status = 1 if succesful submission
end

pause(2) %wait 2 sec before next submission
cd(home_folder);
end
