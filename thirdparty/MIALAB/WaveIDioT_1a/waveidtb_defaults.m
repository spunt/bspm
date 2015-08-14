function waveidtb_defaults
% All the WaveIDT toolbox defaults are stored in this file

global GUI_GLOBAL_DATA;

GUI_GLOBAL_DATA.outDirFlag = 'Default';

% Defualt Output Directory for saving the results
if isunix
    GUI_GLOBAL_DATA.outDir = '/homes/';
else
    GUI_GLOBAL_DATA.outDir = pwd;
end;

% Denoising parameters
GUI_GLOBAL_DATA.WaveType = 'sym2';
GUI_GLOBAL_DATA.NumLev = 2;
GUI_GLOBAL_DATA.dataSelectOption = 2;
GUI_GLOBAL_DATA.origDataPrefix = 'w'; % SPM stores the normalized data with this prefix by default
GUI_GLOBAL_DATA.fileFormat = 'nii';
GUI_GLOBAL_DATA.ParameterMATfileName = '';
GUI_GLOBAL_DATA.WriteMATfileName = '';
GUI_GLOBAL_DATA.timePoints = 1:300;
