function hObject = waveidtb_loadParameters(hObject)

global GUI_GLOBAL_DATA;
waveidtb_defaults;

if isempty(GUI_GLOBAL_DATA.ParameterMATfileName)
    
%     waveidtb_defaults;
    if isunix
        [fileName , direc] = uigetfile('waveidtb_ParameterInfo_*.MAT','/homes/','Select the WaveIDT Parameter MAT file. ');
    else
        [fileName , direc] = uigetfile('waveidtb_ParameterInfo_*.MAT','C:\','Select the WaveIDT Parameter MAT file. ');
    end;

    GUI_GLOBAL_DATA.ParameterMATfileName = fullfile(direc,fileName);
    load(GUI_GLOBAL_DATA.ParameterMATfileName);
    clc;
    fprintf('Denoising Paremeters were loaded from memory.\n');
    disp(['Filename loaded (with path): ',GUI_GLOBAL_DATA.ParameterMATfileName,'.\n']);
    
else
    
    load(GUI_GLOBAL_DATA.ParameterMATfileName);
    fprintf('Denoising Paremeters were loaded from memory.\n');
    fprintf(['Filename loaded (with path): ',GUI_GLOBAL_DATA.ParameterMATfileName,'\n']);
    
end;

guidata(hObject ,  GUI_GLOBAL_DATA);
fprintf('Continue to "Run Denoising".\n');

end
