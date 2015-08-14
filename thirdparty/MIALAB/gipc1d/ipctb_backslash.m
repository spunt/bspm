function [sPath]=ipctb_backslash(sPath)
    if isunix         
        sSlash = '/'; %UNIX
    else
        sSlash = '\'; %PC
    end
    nPathLen = size(sPath,2);
    if char(sPath(nPathLen)) ~= sSlash;        
        if isunix 
        	sPath = [sPath '/']; 
        else 
            % pc or mac
        	sPath = [sPath '\']; 
        end 
        
        
    end
