function handles = validateDataAndParams(handles)

failedTests = zeros(1,7);

if ~( strcmp(handles.Pub.dataDestination(end),'/') || ...
        strcmp(handles.Pub.dataDestination(end),'\') )
    handles.validFlag = false;
    disp(' ')
    disp('Destination directory must end with "slash".')
    disp('Validation failed.')
    failedTests(1) = 1;
    return
end    

handles.validFlag = true;

disp(' ')
disp('Searching and validating fMRI data...')

handles = validateParams(handles,'subj');
if ~handles.validFlag
    disp(' ')
    disp('Validation failed.')
    failedTests(2) = 1;
    return
end

if ( strcmp(handles.Pub.fileFormatSubj,'nii') || handles.Pub.useTemplate == 1 )
    disp(' ')
    disp('Searching nifti-package...')
    if ~( ( exist('load_nii_hdr.m','file') == 2 ) && ...
            ( exist('load_nii_img.m','file') == 2 ) && ...
            ( exist('load_nii.m','file') == 2 ) )
        handles.validFlag = false;
        disp(' ')
        disp('Package to read Nifti-files to Matlab not found.')
        disp('Package can be loaded from: http://www.rotman-baycrest.on.ca/~jimmy/NIFTI/')
        disp('Validation failed.')
        failedTests(4) = 1;
        return
    end

end
if failedTests(4) == 0
    disp('Ok!')
end

disp(' ')
disp('Set private parameters...')
try
    [handles.Priv,handles.Pub] = setPrivParams(handles.Pub);
    %handles.Priv
catch
    handles.validFlag = false;
    disp(' ')
    disp(lasterr)
    disp('Validation failed.')
    failedTests(3) = 1;
    return
end
if failedTests(3) == 0
    disp('Ok!')
end

if handles.checkboxFreq == 1
    disp(' ')
    disp('Checking frequency parameters...')
    handles = validateParams(handles,'freqBands');
    if ~handles.validFlag
        disp(' ')
        disp('Validation failed.')
        failedTests(5) = 1;
        return
    end
end


if handles.Pub.useTemplate == 1
    disp(' ')
    disp('Searching standard templates...')
    
    handles = validateParams(handles,'template');
    if ~handles.validFlag
        disp(' ')
        disp('Validation failed.')
        failedTests(6) = 1;
        return
    end
end

disp(' ')
disp('Checking TR...')
handles = validateParams(handles,'TR');
if ~handles.validFlag
    disp(' ')
    disp('Validation failed.')
    failedTests(7) = 1;
    return
end

if sum(failedTests) == 0
    disp(' ')
    disp('Validation successful!!')
    disp(' ')
    Params.PublicParams = handles.Pub;
    Params.PrivateParams = handles.Priv;
    save([handles.Pub.dataDestination handles.Pub.dataDescription],'Params')
    disp(['Parameter file ''' handles.Pub.dataDescription '.mat'' saved to directory ' handles.Pub.dataDestination '.'])
    handles.Params = Params;
    handles.ParamsValid = Params;
else
    disp(' ')
    disp('Failed tests:')
    if failedTests(1)
        disp('Destination directory')
    end
    if failedTests(2)
        disp('Subject directories')
    end
    if failedTests(3)
        disp('Setting Private parameters')
    end
    if failedTests(4)
        disp('Nifti-file package')
    end
    if failedTests(5)
        disp('Frequency band settings')
    end
    if failedTests(6)
        disp('Template settings')
    end
    
    
end
