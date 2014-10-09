function [] = bspm_convert_dcm(subs, format, force)
% BSPM_CONVERT_DCM
%
% USAGE: bspm_convert_dcm(subs, format, force)
%
%   subs =  paths to subject folders containing dicom directory
%   format  =  'img' or 'nii' (default = 'nii')
%   force = option to force conversion (ignores default directory check)
%

% ------------------------- Copyright (C) 2014 -------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014


% check arguments
if nargin<3, force = 0; end
if nargin<1, disp('USAGE: bspm_convert_dcm(subs, format, force)'); return
elseif nargin<2, disp('Format not specified. Using the default; 3D nifti'); format = 'nii'; end
if ~iscell(subs), subs = cellstr(subs); end

% loop over subject folders
nsubs = length(subs);
for sub = 1:nsubs
    
    subDIR = subs{sub};
    [p, subnam, e] = fileparts(subDIR);
    tmp = files([subDIR filesep '*'], 'dironly', 1);
    
    if ~force
        if length(tmp)~=1
            fprintf('\n%s looks strange... skipping\n', subnam);
            continue
        end
        [p n e] = fileparts(char(tmp));
        if ~strcmp(n,'dicom')
            movefile(tmp{1},[subDIR filesep 'dicom']);
        end
    else
        idx = cellstrfind(tmp,'dicom');
        if isempty(idx), fprintf('\nCouldn''t find dicom directory for %s... skipping\n', subnam); continue; end
    end
    % make directories
    % -----------------
    dirnames = {'raw' 'analysis' 'notes' 'behav'};
    for i = 1:length(dirnames);
        mkdir([subDIR filesep dirnames{i}]);
    end

    dicomDIR = [subDIR filesep 'dicom'];
    dothese = files([dicomDIR filesep '*']);

    for s = 1:length(dothese)

        currentDIR = dothese{s};
        dicomfiles = files([currentDIR filesep '*.dcm']);
        if isempty(dicomfiles)
            dicomfiles = files([currentDIR filesep '*']);
        end
        current_hdr = spm_dicom_headers(dicomfiles{1});

        % get some info from header file
        studyStartTime=current_hdr{1}.PerformedProcedureStepID(end-5:end-2);
        sequenceTime=round((current_hdr{1}.AcquisitionTime-current_hdr{1}.StudyTime)/60);
        sequenceTime=num2str(sequenceTime);
        sequenceName=regexprep(current_hdr{1}.ProtocolName,' ','_');
        sequenceType=regexprep(current_hdr{1}.ScanningSequence,'\','_');
        sequenceType=regexprep(sequenceType,' ','_');

        % now...
        [p n e] = fileparts(currentDIR);
        outputDIR = [subDIR filesep 'raw' filesep sequenceType '_' n];
        mkdir(outputDIR);
        hdr = current_hdr;
        save([outputDIR filesep 'dicom_header.mat'], 'hdr');

        % this is what the dicom convert job looks like:
        matlabbatch{1}.spm.util.dicom.data = cellstr(dicomfiles);
        matlabbatch{1}.spm.util.dicom.root = 'flat';
        matlabbatch{1}.spm.util.dicom.outdir{1} = outputDIR;
        matlabbatch{1}.spm.util.dicom.convopts.format = format;
        matlabbatch{1}.spm.util.dicom.convopts.icedims = 0;

        % run
        spm('defaults','fmri'); spm_jobman('initcfg');  
        spm_jobman('run',matlabbatch);
        clear matlabbatch

    end;

end;
 
 
 
 
