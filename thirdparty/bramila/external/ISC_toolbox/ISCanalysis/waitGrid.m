function waitGrid(grid_type, pausetime, log_path)
%wait_results=waitGrid(grid_type, pausetime)
% grid_type is from testGrig function
% pausetime defines in seconds how long functon waits before it tests again
% if processes have finished.
% when when all processes are ready function returns false
%
% Created 05.08.2013 
% Juha Pajula, Tampere University of Technology 
% juha.pajula@tut.fi 

wait_results=true;
[~,us_name] = unix('echo $LOGNAME');

disp('Waiting processes to end: ')
while(wait_results)
    pause(pausetime);
    
    %if the grid is SGE
    if(strcmp(grid_type,'sge'))
        [~,outpt]=unix(['qstat -u ', us_name]); %status = 1 if succesful submission
        if(isempty(outpt))
            %wait_results=false;
            return %if ended return
        end
    end

    %if the grid is Slurm
    if(strcmp(grid_type,'slurm'))
        [~,outpt]=unix(['squeue -h -u ', us_name]); %status = 1 if succesful submission
        if(isempty(outpt))
            %wait_results=false;
            return %if ended return
        end
    end
   error_files=dir([log_path, '/*.e*']); %searching the error logs
   %error_files = sum(cell2mat({error_files(:).bytes}));
   %if error files have some content, trigger an error 
   if(sum(cell2mat({error_files(:).bytes}))~=0)
        error('Error occured in grid computing. See the *.e* files in /scripts folder for details')   
   end
     
    %print dots to show advance in process
    fprintf('\b',1); disp('.')
    
end
