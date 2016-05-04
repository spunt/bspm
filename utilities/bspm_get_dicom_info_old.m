function dicominfo = bspm_get_dicom_info(in,disptag)
% BSPM_GET_DICOM_INFO
%
% USAGE: dicominfo = bspm_get_dicom_info(in,disptag)
%
%   ARGUMENTS
%       in = dicom file
%       disptag = 1 (default) will display (requires f(n) strucdisp)
%
%   OUTPUT EXAMPLE
%       dicominfo.parameterinfo.TR = 2500;
%       dicominfo.parameterinfo.voxelsize = 3;
%       dicominfo.parameterinfo.matrixsize = 64;
%       dicominfo.parameterinfo.echotime = 30;
%       dicominfo.parameterinfo.flipangle = 80;
%       dicominfo.parameterinfo.bandwidth = 2604;
%       dicominfo.sequenceinfo.name = TOM;
%       dicominfo.sequenceinfo.type = EP;
%       dicominfo.sequenceinfo.pulsename = *epfid2d1_64;
%       dicominfo.sequenceinfo.timestamp = MR20130328094606;
%       dicominfo.sequenceinfo.order = 8;
%       dicominfo.subjectinfo.subjectid = AM_032813;
%       dicominfo.subjectinfo.age = 039Y;
%       dicominfo.subjectinfo.sex = F;
%       dicominfo.sliceinfo.spacing = 3;
%       dicominfo.sliceinfo.orientation = Tra>Cor(-21.2);
%       dicominfo.sliceinfo.acquisitiontimes[1] = 1252.5;
%       dicominfo.sliceinfo.order
%       dicominfo.sliceinfo.number = 46;
%

% ----------------------- Copyright (C) 2014 -----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, mfile_showhelp; return; end
if nargin < 2, disptag = 1; end
if iscell(in), in = char(in); end

% check file
try 
    load(in);
catch
    hdr = spm_dicom_headers(in);
end
allfield = structfields(hdr{1});

% scanner information
if isfield(hdr{1},'.ManufacturersModelName'), dicominfo.scannerinfo.model = hdr{1}.ManufacturersModelName; end
if cellstrfind(allfield,'32Ch'), dicominfo.scannerinfo.coil = 32; else dicominfo.scannerinfo.coil = 12; end
if isfield(hdr{1}, 'InstitutionName'), dicominfo.scannerinfo.facility = hdr{1}.InstitutionName; end
if isfield(hdr{1}, 'InstitutionAddress'), dicominfo.scannerinfo.location = hdr{1}.InstitutionAddress; end
    
% general information
dicominfo.parameterinfo.TR          = hdr{1}.RepetitionTime;
dicominfo.parameterinfo.voxelsize   = [hdr{1}.PixelSpacing' hdr{1}.SliceThickness];
dicominfo.parameterinfo.matrixsize  = hdr{1}.AcquisitionMatrix;
dicominfo.parameterinfo.echotime    = hdr{1}.EchoTime;
dicominfo.parameterinfo.flipangle   = hdr{1}.FlipAngle;
dicominfo.parameterinfo.bandwidth   = hdr{1}.PixelBandwidth;

% dicominfo.sequenceinfo
dicominfo.sequenceinfo.name = cleanupname(hdr{1}.ProtocolName); 
dicominfo.sequenceinfo.type = cleanupname(hdr{1}.ScanningSequence);
dicominfo.sequenceinfo.acquisitiontype = cleanupname(hdr{1}.MRAcquisitionType); 
dicominfo.sequenceinfo.pulsename = cleanupname(hdr{1}.SequenceName);
dicominfo.sequenceinfo.timestamp = cleanupname(hdr{1}.PerformedProcedureStepID); 
dicominfo.sequenceinfo.order = hdr{1}.SeriesNumber;
idx = cellstrfind(allfield,'TotalScanTimeSec');
tmp = allfield{idx};
idx = strfind(tmp,'TotalScanTimeSec');
tmp = tmp(idx:idx+100);
idx = strfind(tmp,'=');
tmp = tmp(idx(1)+1:idx(1)+5);
dicominfo.sequenceinfo.duration_secs = str2double(tmp);
str = 'sPat.lAccelFactPE                        = ';
idx = cellstrfind(allfield,str);
tmp = allfield{idx};
idx = strfind(tmp,str);
idx = idx+length(str);
dicominfo.sequenceinfo.ipatfactor = str2num(tmp(idx:idx+1));
str = 'dReadoutFOV';
idx = cellstrfind(allfield,str);
tmp = allfield{idx};
idx = strfind(tmp,str);
tmp = tmp(idx:idx+50);
idx = strfind(tmp,'=');
tmp = strtrim(tmp(idx(1)+1:idx(1)+5));
dicominfo.sequenceinfo.FOVread = str2num(tmp);

% dicominfo.subjectinfo
if isfield(hdr{1},'PatientName')
    dicominfo.subjectinfo.subjectid = strtrim(hdr{1}.PatientName);
    dicominfo.subjectinfo.age = strtrim(hdr{1}.PatientAge);
    dicominfo.subjectinfo.sex = strtrim(hdr{1}.PatientSex);
else
    dicominfo.subjectinfo.subjectid = strtrim(hdr{1}.PatientsName);
    dicominfo.subjectinfo.age = strtrim(hdr{1}.PatientsAge);
    dicominfo.subjectinfo.sex = strtrim(hdr{1}.PatientsSex);
end

% dicominfo.sliceinfo
if strcmp(hdr{1}.MRAcquisitionType, '2D')
    

    
    dicominfo.sliceinfo.effechospacing = get_echo_spacing(hdr); 
    
    if isfield(hdr{1},'Slicethickness'), dicominfo.sliceinfo.thickness = hdr{1}.SliceThickness; end
    if isfield(hdr{1},'SpacingBetweenSlices'), dicominfo.sliceinfo.spacing = hdr{1}.SpacingBetweenSlices; end
    if isfield(hdr{1},'Private_0051_100e'), dicominfo.sliceinfo.orientation = strtrim(hdr{1}.Private_0051_100e); end
    if isfield(hdr{1},'Private_0019_1029')
        slicetimes = hdr{1}.Private_0019_1029';
        dicominfo.sliceinfo.acquisitiontimes = slicetimes;
        slicetimes(:,2) = 1:length(slicetimes);
        slicetimes = sortrows(slicetimes,1);
        dicominfo.sliceinfo.order = slicetimes(:,2);
        dicominfo.sliceinfo.number = length(dicominfo.sliceinfo.acquisitiontimes);
    end
    
    
end

% display
if disptag
    tmp = which('strucdisp.m'); 
    if ~isempty(tmp), strucdisp(dicominfo); end
end

end

% | SUBFUNCTIONS

function cname = cleanupname(name); 
    cname = regexprep(strtrim(name), ' ', '_');
    cname = regexprep(cname, '\', '');
end

% - everything below adapted from dicm2nii.m
% - see: http://www.mathworks.com/matlabcentral/fileexchange/42997
function dwell = get_echo_spacing(hdr)
   


    s = hdr{1}.CSAImageHeaderInfo;
    csaidx = structfind(s, 'name', 'BandwidthPerPixelPhaseEncode');
    hz = s(csaidx).item(1).val; 

    
    hz = csa_header(s, 'BandwidthPerPixelPhaseEncode');
    dwell = 1000 ./ hz / dim(iPhase); % in ms
    if isempty(dwell) % true for syngo MR 2004A
        % ppf = [1 2 4 8 16] represent [4 5 6 7 8] 8ths PartialFourier
        % ppf = asc_header(s, 'sKSpace.ucPhasePartialFourier');
        lns = asc_header(s, 'sKSpace.lPhaseEncodingLines');
        dur = csa_header(s, 'SliceMeasurementDuration');
        dwell = dur ./ lns; % ./ (log2(ppf)+4) * 8;
    end
    if isempty(dwell) % next is not accurate, so as last resort
        dur = csa_header(s, 'RealDwellTime') * 1e-6; % ns to ms
        dwell = dur * asc_header(s, 'sKSpace.lBaseResolution');
    end
    if isempty(dwell)
        dwell = double(tryGetField(s, 'EffectiveEchoSpacing')) / 1000; % GE
    end
    % http://www.spinozacentre.nl/wiki/index.php/NeuroWiki:Current_developments
    if isempty(dwell) % Philips
        wfs = tryGetField(s, 'WaterFatShift');
        epiFactor = tryGetField(s, 'EPIFactor');
        dwell = wfs ./ (434.215 * (double(epiFactor)+1)) * 1000;
    end
    if ~isempty(dwell)
        if s.isDTI
            readout = dwell * dim(iPhase) / 1000; % in sec
            descrip = sprintf('readout=%.3g;%s', readout, descrip);
            s.ReadoutSeconds = readout;
        else
            descrip = sprintf('dwell=%.3g;%s', dwell, descrip);
            s.EffectiveEPIEchoSpacing = dwell;
        end
    end
end



 
 
 
 
