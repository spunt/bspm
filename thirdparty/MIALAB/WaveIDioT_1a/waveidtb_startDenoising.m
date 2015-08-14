function waveidtb_startDenoising(hObject)
%% Wavelet based 3D denoising of real fMRI data for group data
% Created on: July 26 2010
% Revised on: April 29 2011
% Author: Siddharth Khullar

% Description: Wavelet_denoising_main_RealData.m to handle group data
% History: Version 1
%          Data Structure added on 07/28/2010.

global GUI_GLOBAL_DATA;

ParameterInfo.waveType = GUI_GLOBAL_DATA.WaveType; % Type of Wavelet Filters used.
ParameterInfo.levels = GUI_GLOBAL_DATA.NumLev; % Number of levels of denoising
ParameterInfo.OutPath = GUI_GLOBAL_DATA.outDir; % Output path where denoised data is stored
ParameterInfo.fileNames = GUI_GLOBAL_DATA.fileNames;
ParameterInfo.fileFormat = GUI_GLOBAL_DATA.fileFormat;

%% Mask generation (Using the first subject's data)
[mask ,x , y , z, w ] = waveidtb_generateMask(ParameterInfo);

%% Loop to run through subjects
denoisedMap = zeros(x,y,z);
fprintf('A total of %2.1i images were selected.\n' , length(ParameterInfo.fileNames));
volNo = 1;
for n = 1:length(ParameterInfo.fileNames)

    fprintf('Denoising image no. %1.1i' , n);
    % Read the image
    img = spm_read_vols( spm_vol( ParameterInfo.fileNames{n} ) );
    denoisedMap(:,:,:,volNo) = waveidtb_waveletDenoise(img , mask , ParameterInfo);
    clear img;
    % Write after w volumes have been denoised
    if volNo == w
        fprintf('Writing the denoised data for this subject.');
        ParameterInfo.filesWritten = ParameterInfo.fileNames((n-w+1):n , 1);
        waveidtb_writeData(denoisedMap, ParameterInfo , n )
        fprintf('Done for this subject/session (%2.1i time points ).\n',volNo);
        volNo = 0; % Reset Volume no.
        
    end;
    
    volNo = volNo + 1;
end;
