function hObject = waveidtb_getOutputDir(hObject)

global GUI_GLOBAL_DATA;
waveidtb_defaults;
clc;
fprintf('Setup NEW wavelet denoising analysis.\n');
if isunix
    [ GUI_GLOBAL_DATA.outDir ] = uigetdir('/homes/','Select directory to store the Denoising Parameters (MAT file)');
else
    [ GUI_GLOBAL_DATA.outDir ] = uigetdir('C:\','Select directory to store the Denoising Parameters (MAT file)');
end;

% If the user chooses a directory, set this flag otherwise leave it to
% default;
if GUI_GLOBAL_DATA.outDir 
    [ GUI_GLOBAL_DATA.outDirFlag ] = 'UserDefined';
end;
guidata(hObject,GUI_GLOBAL_DATA);
end