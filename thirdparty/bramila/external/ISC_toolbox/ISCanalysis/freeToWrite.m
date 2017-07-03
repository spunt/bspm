function status = freeToWrite(mode,lock_path,pref)
% function writes and destroys a lock file depending on the mode
% Usage: status = freeToWrite(lock_path,mode,pref)
% 
% lock_path = full path where the lock is written/will be written
% pref = string containing the filename without ending ".lock"
%
% mode = 'check' -> function stays polling if the file exists, when it does
% not anymore exists it writes a new and returns true
%
% mode = 'release' -> function checks that the lock is there and then
% deletes the lock file
%
%Example: if(freeToWrite('check',full_path,'process1')
%            disp('the lock is set')
%            save('Params')
%            [~]=freeToWrite('release',full_path,'process1')
%         end
%
% Juha Pajula, 2013, juha.pajula@tut.fi

switch mode
    case 'check'
        while(size(dir([lock_path,'*.lock']),1)~=0)
            pause(ceil(rand(1)*100)/10)
            disp('Waiting to get the lock.')
        end
        dlmwrite([lock_path,pref, '.lock'], ' ')
        status = true;
        disp(['Set Lock: ', pref, '.lock']);
        return
        
    case 'release'
        if(size(dir([lock_path,pref,'.lock']),1)==1)
            delete([lock_path,pref,'.lock'])
            status=true;
            disp(['Released Lock: ', pref, '.lock']);
            return
        else
            error('Lock not found!')
        end
    otherwise
        error('Some error within lock system.')
end

