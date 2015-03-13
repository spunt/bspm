function matlabbatch = bspm_convert_dcm(subdirs)
% BSPM_CONVERT_DCM
%
% USAGE: matlabbatch = bspm_convert_dcm(subdirs)
%
%   subdirs = paths to subject folders containing dicom directory
%

% ------------------------- Copyright (C) 2014 -------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<1, disp('USAGE: bspm_convert_dcm(subdirs)'); return; end
if ischar(subdirs), subdirs = cellstr(subdirs); end
nsubs = length(subdirs);
count = 0; 
for sub = 1:nsubs

    subDIR          = subdirs{sub};
    [p, subnam, e]  = fileparts(subDIR);
    tmp = files([subDIR filesep '*'], 'dironly', 1);
    [p, n, e] = fileparts(char(tmp));
    if ~strcmp(n,'dicom'), movefile(tmp{1},[subDIR filesep 'dicom']); end
    dirnames = {'raw' 'analysis' 'notes' 'behav'};
    for i = 1:length(dirnames), mkdir([subDIR filesep dirnames{i}]); end
    dicomDIR    = fullfile(subDIR, 'dicom');
    dothese     = files(fullfile(dicomDIR, '*'));
    for s = 1:length(dothese)

        currentDIR = dothese{s};
        dicomfiles = files([currentDIR filesep '*.dcm']);
        if isempty(dicomfiles), dicomfiles = files([currentDIR filesep '*']); end
        dcmref  = dicomfiles{floor(length(dicomfiles)/2)}; 
        dcminfo = bspm_get_dicom_info(dcmref, 0); 
        outputDIR = fullfile(subDIR, 'raw', sprintf('%s_%s_%02d', dcminfo.sequenceinfo.type(1:2), dcminfo.sequenceinfo.name, dcminfo.sequenceinfo.order));
        mkdir(outputDIR);
        copyfile(dcmref, fullfile(outputDIR, 'dicom_reference.dcm'));
        if strcmpi(dcminfo.sequenceinfo.type(1:2), 'EP')
            save(fullfile(outputDIR, 'dicom_info.mat'), 'dcminfo');
        end
        
        % | Job
        count = count + 1; 
        matlabbatch{count}.spm.util.import.dicom.data = cellstr(dicomfiles);
        matlabbatch{count}.spm.util.import.dicom.root = 'flat';
        matlabbatch{count}.spm.util.import.dicom.outdir{1} = outputDIR;
        matlabbatch{count}.spm.util.import.dicom.protfilter = '.*';
        matlabbatch{count}.spm.util.import.dicom.convopts.format = 'nii';
        matlabbatch{count}.spm.util.import.dicom.convopts.icedims = 0;

    end;

end;
% run
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end
end
