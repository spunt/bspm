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
    'funcprefix',       'fadolphs', ...
    'anatprefix',       'sadolphs'  ...
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
    if isempty(alldcm), fprintf('\n!!! NO .DCM FILES FOR THIS SUBJECT, MOVING ON...'); continue; end
    [dothese, ia, ic] = unique(cellfun(@fileparts, alldcm, 'unif', false));
    fprintf('\n| Running DICM2NII on %d dicom folders\n', length(dothese));
    disp(dothese);
    rawDIR = fullfile(subDIR, 'raw'); mkdir(rawDIR); 
    fprintf('| Output directory: %s\n', rawDIR);
    
    seqid = cell(length(dothese), 2);
    
    for s = 1:length(dothese)

        dicomfiles  = alldcm(ic==s); 
        dcmref      = dicomfiles{ceil(length(dicomfiles)/2)}; 
        dcminfo     = bspm_get_dicom_info(dcmref, 0); 
        outputDIR   = fullfile(rawDIR, sprintf('%s_%s_%02d', dcminfo.sequenceinfo.type(1:2), dcminfo.sequenceinfo.name, dcminfo.sequenceinfo.order));
        mkdir(outputDIR);
        dicm2nii(dothese{s}, outputDIR, 4);
        try
            dcminfo.dcmHeaders = load(fullfile(outputDIR, 'dcmHeaders.mat'));
            delete(fullfile(outputDIR, 'dcmHeaders.mat'));
        catch
        end
        tmp = dcminfo.sequenceinfo; 
        seqid{s,1} = sprintf('%s%s%s%d', tmp.name, tmp.type, tmp.pulsename, tmp.timestamp);
        seqid{s,2} = outputDIR; 
        save(fullfile(outputDIR, 'dicominfo.mat'), 'dcminfo');
        nii         = files(fullfile(outputDIR, '*nii'));
        
        if strcmp(tmp.type, 'EP')
           if omitfirstN, delete(nii{1:omitfirstN}); nii(1:omitfirstN) = []; end
           if ~isempty(funcprefix)
                [pth, fn, fe] = cellfileparts(nii);
                outfn = strcat(pth, filesep, sprintf('%s_%03d_', funcprefix, dcminfo.sequenceinfo.order), lower(fn), fe);
                cellfun(@movefile, nii, outfn)
           end
           if do3dto4d, bnii_3dto4d(nii, 'compress', 1, 'delete3d', 1, 'delimiter', '_'); end
        elseif all(strncmp(tmp.type, 'GR', 2), ~isempty(anatprefix))
            [pth, fn, fe] = cellfileparts(nii);
            
            outfn = strcat(pth, filesep, sprintf('%s_%03d_', anatprefix, dcminfo.sequenceinfo.order), lower(fn), fe);
            cellfun(@movefile, nii, outfn)
        end 
        if s~=1
            if all([strcmp(tmp.type, 'GR') strcmp(seqid{s,1}, seqid{s-1,1})])
                nii         = files(fullfile(outputDIR, '*nii'));
                movefile(char(nii), seqid{s-1,2});
                rmdir(outputDIR, 's');
            end
            
        end
        
        
    end
    
    % | Change name of parent dicom folder to "dicom"
    try
        movefile(parentpath(dothese), fullfile(subDIR, 'dicom')); 
    catch
    end
    
    % | Make Additional Directories
    for i = 1:length(adddirs), mkdir([subDIR filesep adddirs{i}]); end

end;
end
