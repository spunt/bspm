function bDirCreated = ipctb_fCreateDir(sDir, bDelExisting)

[retSucc retMess retId] = mkdir(sDir);
if strcmp(retId, 'MATLAB:MKDIR:DirectoryExists')
    sRet = questdlg(['OK to erase previous files in Folder ' [sDir] ' ?'],'FNC','Yes','No','Yes');
    if strcmp(sRet(1), 'Y')
        delete ([ipctb_backslash(sDir) '*']);
        bDirCreated = true;
    else
        bDirCreated = false;
    end
else
    bDirCreated = true;    
end