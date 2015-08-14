    function waveidtb_writeData(denoisedMap, ParameterInfo , n )
    
    global GUI_GLOBAL_DATA;
    
    V = spm_vol( ParameterInfo.filesWritten{n} );
    Vmat = V.mat;
    WritingParamters = struct;
    WritingParamters.Vmat = Vmat;
    WritingParameters.fileNames = ParameterInfo.filesWritten;
    WritingParameters.waveType = ParameterInfo.waveType; % Type of Wavelet Filters used.
    WritingParameters.levels = ParameterInfo.levels; % Number of levels of denoising
    WritingParameters.DenoisingInfoDirectory = ParameterInfo.OutPath;
    [direc,temp] = fileparts(ParameterInfo.fileNames{n});
    fileName = fullfile(direc,['WD',temp,'.nii']);
    fprintf('\nWriting the denoised volumes to a NIFTI file.\n');
    waveidtb_create4DNiftifile(fileName, denoisedMap, Vmat);
    fprintf('NIFTI file saved as ');disp(fileName);
    fprintf('Saving Write parameters as MAT File.\n');
    MATfileName = fullfile(direc,['WaveIDT_',temp,'_',datestr(now,'ddmmmyyyy_HHMMSS')]);
    save( MATfileName ,'WritingParameters');
    fprintf('MAT file saved as : ');disp(MATfileName);
    GUI_GLOBAL_DATA.WriteMATfileName = MATfileName;