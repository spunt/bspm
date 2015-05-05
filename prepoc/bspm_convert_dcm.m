function bspm_convert_dcm(subdirs, varargin)
% BSPM_CONVERT_DCM
%
% USAGE: bspm_convert_dcm(subdirs, varargin)
%
%   subdirs = paths to subject folders containing dicom directory
%

% ------------------------- Copyright (C) 2014 -------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
def = { ...
    'dcmdirpat',        '*', ...
	'do3dto4d',         0,	...
	'omitfirstN',		0,	...
    'outputdir',        'raw', ...
    'adddirs',          {'analysis' 'behav'} ...
	};
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
if ischar(subdirs), subdirs = cellstr(subdirs); end
nsubs = length(subdirs);
for sub = 1:nsubs

    subDIR          = subdirs{sub};
    [p, subnam, e]  = fileparts(subDIR);
    fprintf('\n| Working on: %s', subnam);
    alldcm          = files(fullfile(subDIR, '**', '*dcm'));
    if isempty(alldcm)
        fprintf('\n!!! COULD NOT FIND .DCM FILES FOR THIS SUBJECT, MOVING ON !!!'); 
        continue
    end
    [dothese, ia, ic] = unique(cellfun(@fileparts, alldcm, 'unif', false));
    fprintf('\n| Running DICM2NII on %d dicom folders\n', length(dothese));
    disp(dothese);
    rawDIR = fullfile(subDIR, 'raw'); mkdir(rawDIR); 
    fprintf('| Output directory: %s\n', rawDIR);
    for s = 1:length(dothese)
        
        dicomfiles  = alldcm(ic==s); 
        dcmref      = dicomfiles{floor(length(dicomfiles)/2)}; 
        dcminfo     = bspm_get_dicom_info(dcmref, 0); 
        outputDIR   = fullfile(rawDIR, sprintf('%s_%s_%02d', dcminfo.sequenceinfo.type(1:2), dcminfo.sequenceinfo.name, dcminfo.sequenceinfo.order));
        mkdir(outputDIR);
        dicm2nii(dothese{s}, outputDIR, 4);
        if strcmp(dcminfo.sequenceinfo.type, 'EP')
           nii = files(fullfile(outputDIR, '*nii')); 
           if omitfirstN, delete(nii{1:omitfirstN}); nii(1:omitfirstN) = []; end
           if do3dto4d
               bnii_3dto4d(nii, 'compress', 1, 'delete3d', 1, 'delimiter', '_'); 
           end
        end 
    end
    
    % | Change name of parent dicom folder to "dicom"
    if ~strcmpi(parentpath(dothese), 'dicom'), 
        movefile(parentpath(dothese), fullfile(subDIR, 'dicom')); 
    end
    
    % | Make Additional Directories
    for i = 1:length(adddirs), mkdir([subDIR filesep adddirs{i}]); end

end;
end
