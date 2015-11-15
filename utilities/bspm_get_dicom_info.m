function varargout = bspm_get_dicom_info(in, disptag)
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

if nargin < 1, error('USAGE: bspm_get_dicom_info(in, disptag)'); end
if nargin < 2, disptag = 1; end
if iscell(in), in = char(in); end

% | Load file
[p, n, e] = fileparts(in);
if strcmpi(e, '.mat')
    h = load(in);
    hdr = getnestedfield(h, 'hdr');
    if iscell(hdr), hdr = hdr{1}; end
    hdr = orderfields(hdr); 
else
    hdr = orderfields(dicm_hdr(in)); 
end
allfn   = structfields(hdr);

% scanner information
dicominfo.scannerinfo = struct( ...
    'manufacturer', hdr.('Manufacturer'), ...
    'model', hdr.('ManufacturerModelName'), ...
    'facility', hdr.('InstitutionName'), ...
    'location', hdr.('InstitutionAddress'), ...
    'fieldstrength', hdr.('MagneticFieldStrength'));
if cellstrfind(allfn,'32Ch')
    dicominfo.scannerinfo.coil = 32; 
elseif cellstrfind(allfn,'12Ch')
    dicominfo.scannerinfo.coil = 12; 
end

% general information
dicominfo.parameterinfo = struct( ...
    'TR', hdr.('RepetitionTime'), ...
    'voxelsize', [hdr.('PixelSpacing')' hdr.('SliceThickness')], ...
    'echotime', hdr.('EchoTime'), ...
    'echonumber', [], ...
    'flipangle', hdr.('FlipAngle'), ...
    'bandwidth', hdr.('PixelBandwidth'));
try dicominfo.parameterinfo.echonumber = hdr.('EchoNumber'); catch; end

% dicominfo.subjectinfo
dicominfo.subjectinfo = struct( ...
    'subjectid', hdr.('PatientName'), ...
    'age', hdr.('PatientAge'), ...
    'sex', hdr.('PatientSex'), ...
    'weight', hdr.('PatientWeight'));

% dicominfo.sequenceinfo
dicominfo.sequenceinfo = struct( ...
    'name', hdr.('ProtocolName'), ...
    'type', hdr.('ScanningSequence'), ...
    'imagetype', hdr.('MRAcquisitionType'), ...
    'pulsename', hdr.('SequenceName'), ...
    'timestamp', hdr.('PerformedProcedureStepID'), ...
    'order', hdr.('SeriesNumber'), ...
    'acquisitiontime', hdr.('AcquisitionTime'), ...
    'seriestime', hdr.('SeriesTime'), ...
    'studytime', hdr.('StudyTime'));

idx = cellstrfind(allfn,'TotalScanTimeSec');
tmp = allfn{idx};
idx = strfind(tmp,'TotalScanTimeSec');
tmp = tmp(idx:idx+100);
idx = strfind(tmp,'=');
tmp = tmp(idx(1)+1:idx(1)+5);
dicominfo.sequenceinfo.duration_secs = str2double(tmp);
str = 'sPat.lAccelFactPE                        = ';
idx = cellstrfind(allfn,str);
tmp = allfn{idx};
idx = strfind(tmp,str);
idx = idx+length(str);
dicominfo.sequenceinfo.ipatfactor = str2num(tmp(idx:idx+1));
str = 'dReadoutFOV';
idx = cellstrfind(allfn,str);
tmp = allfn{idx};
idx = strfind(tmp,str);
tmp = tmp(idx:idx+50);
idx = strfind(tmp,'=');
tmp = strtrim(tmp(idx(1)+1:idx(1)+5));
dicominfo.sequenceinfo.FOVread = str2num(tmp);

% dicominfo.sliceinfo
if strcmp(hdr.ScanningSequence, 'EP')
    
    % | Volume Dimensions
    try
        dicominfo.parameterinfo.dim = double([hdr.('AcquisitionMatrix')([1 4])' length([hdr.('MosaicRefAcqTimes')])]);

    
    % | Slice Timing Information
    dicominfo.sliceinfo = struct( ...
    'thickness', hdr.('SliceThickness'), ...
    'spacebetweenslices', hdr.('SpacingBetweenSlices'), ...
    'orientation', hdr.('ImageOrientation'));
    slicetimes = double(hdr.('MosaicRefAcqTimes'));
    dicominfo.sliceinfo.acquisitiontimes = slicetimes;
    slicetimes(:,2) = 1:length(slicetimes);
    slicetimes = sortrows(slicetimes,1);
    dicominfo.sliceinfo.order = slicetimes(:,2);
    dicominfo.sliceinfo.number = length(dicominfo.sliceinfo.acquisitiontimes);
    [es, uwdir] = get_echo_spacing(hdr, dicominfo);
    
    % | UNWARP INFORMATION
    dicominfo.unwarpinfo.effechospacing     = es;
    dicominfo.unwarpinfo.readouttime        = es*hdr.('AcquisitionMatrix')(1);
    dicominfo.unwarpinfo.unwarpdirection    = uwdir; 
    
    catch

    end

end

% display
if disptag
    tmp = which('strucdisp.m'); 
    if ~isempty(tmp), strucdisp(dicominfo); end
end
if nargout > 0, varargout{1} = dicominfo; end
    

end
% | SUBFUNCTIONS (adapted from code from DICM2NII)
function [es, uwdir] = get_echo_spacing(s, d)
    
    dim = d.parameterinfo.dim; 
    [ixyz, R, pixdim, xyz_unit] = xform_mat(s, dim);
    R(1:2,:) = -R(1:2,:); % dicom LPS to nifti RAS, xform matrix before reorient
    [phPos, fps_bits] = phaseDirection(s);
    [~, perm] = sort(ixyz); % may permute 3 dimensions in this order
    if (strcmp(tryGetField(s, 'MRAcquisitionType'), '3D'))
        dim(3)>1 && (~isequal(perm, 1:3)) % skip if already standard view
        R(:, 1:3) = R(:, perm); % xform matrix after perm
        fps_bits = fps_bits(perm);
        ixyz = ixyz(perm); % 1:3 after re-orient
        dim = dim(perm);
        pixdim = pixdim(perm);
        if isfield(s, 'bvec'), s.bvec = s.bvec(:,perm); end
    end
    iSL = find(fps_bits==16);
    iPhase = find(fps_bits==4); % axis index for phase_dim in re-oriented img
    % Flip image to make first axis negative and other two positive
    ind4 = ixyz + [0 4 8]; % index in 4xN matrix
    flip = R(ind4)<0; % flip an axis if true
    flip(1) = ~flip(1); % first axis negative: comment this to make all positive
    rotM = diag([1-flip*2 1]); % 1 or -1 on diagnal
    rotM(1:3, 4) = (dim-1) .* flip; % 0 or dim-1
    R = R / rotM; % xform matrix after flip
    if flip(iPhase), phPos = ~phPos; end
    hz = s.('BandwidthPerPixelPhaseEncode');
    es = (1000./hz) / dim(iPhase); % in ms
    if isempty(phPos), pm = '?'; elseif phPos, pm = ''; else pm = '-'; end
    axes = 'xyz';
    uwdir = [pm axes(iPhase)]; 
end
function [phPos, fps_bits] = phaseDirection(s)
    iPhase = 2; % COLUMN
    foo = tryGetField(s, 'InPlanePhaseEncodingDirection', ''); % no for Philips
    if strcmp(foo, 'ROW'), iPhase = 1; end
    ixyz = xform_mat(s);
    iPhase = ixyz(iPhase); % now can be 1/2/3 for RL/AP/IS
    phPos = [];
    if strncmpi(s.Manufacturer, 'SIEMENS', 7)
        phPos = csa_header(s, 'PhaseEncodingDirectionPositive'); % 0 or 1
    elseif strncmpi(s.Manufacturer, 'GE', 2)
        fld = 'ProtocolDataBlock';
        if isfield(s, fld) && isfield(s.(fld), 'VIEWORDER')
            phPos = s.(fld).VIEWORDER == 1; % 1 == bottom_up
        end
    elseif strncmpi(s.Manufacturer, 'Philips', 7) % no InPlanePhaseEncodingDirection
        fld = 'MRStackPreparationDirection';
        if isfield(s, 'Stack') && isfield(s.Stack.Item_1, fld)
            d = s.Stack.Item_1.(fld);
            iPhase = strfind('LRAPSIFH', d(1));
            iPhase = ceil(iPhase/2); 
            if iPhase>3, iPhase = 3; end % 1/2/3 for RL/AP/IS
            if any(d(1) == 'LPHS'), phPos = false;
            elseif any(d(1) == 'RAFI'), phPos = true;
            end
        end
    end
    fps_bits = [1 4 16]; % 4 for phase_dim
    if iPhase == ixyz(1), fps_bits = [4 1 16]; end
end
function val = csa_header(s, key)
try val = s.CSAImageHeaderInfo.(key); 
catch, val = [];
end
 end
function val = tryGetField(s, field, dftVal)
if isfield(s, field), val = s.(field); 
elseif nargin>2, val = dftVal;
else val = [];
end
end
function [ixyz, R, pixdim, xyz_unit] = xform_mat(s, dim)
R = reshape(tryGetField(s, 'ImageOrientationPatient', [1 0 0 0 1 0]), 3, 2);
R(:,3) = cross(R(:,1), R(:,2)); % right handed, but sign may be opposite
foo = abs(R);
[~, ixyz] = max(foo); % orientation info: perm of 1:3
if ixyz(2) == ixyz(1), foo(ixyz(2),2) = 0; [~, ixyz(2)] = max(foo(:,2)); end
if any(ixyz(3) == ixyz(1:2)), ixyz(3) = setdiff(1:3, ixyz(1:2)); end
if nargout<2, return; end
iSL = ixyz(3); % 1/2/3 for Sag/Cor/Tra slice
cosSL = R(iSL, 3);

thk = tryGetField(s, 'SpacingBetweenSlices');
if isempty(thk), thk = tryGetField(s, 'SliceThickness'); end
pixdim = tryGetField(s, 'PixelSpacing');
if isempty(thk) || isempty(pixdim), xyz_unit = 0; else xyz_unit = 2; end % mm
if isempty(thk), thk = 1; end
if isempty(pixdim), pixdim = [1 1]; end
pixdim = [pixdim(:); thk];
R = R * diag(pixdim); % apply vox size
% Next is almost dicom xform matrix, except mosaic trans and unsure slice_dir
R = [R tryGetField(s, 'ImagePositionPatient', [0 0 0]'); 0 0 0 1];

% rest are former: R = verify_slice_dir(R, s, dim, iSL)
if dim(3)<2, return; end % don't care direction for single slice
pos = []; % SliceLocation for last or center slice we try to retrieve

if s.Columns > dim(1) % Siemens mosaic: use Columns since no transpose to img
    R(:,4) = R * [ceil(sqrt(dim(3))-1)*dim(1:2)/2 0 1]'; % real slice location
    vec = csa_header(s, 'SliceNormalVector'); % mosaic has this
    if ~isempty(vec) % exist for all tested data
        if sign(vec(iSL)) ~= sign(cosSL), R(:,3) = -R(:,3); end
        return;
    end
    
    % only reach here if SliceNormalVector not exists
    if isfield(s, 'CSASeriesHeaderInfo')
        i = 1; % use 2nd slice unless revNumbering & dim(3)==2
        if dim(3)==2 && csa_header(s, 'ProtocolSliceNumber'), i = 0; end
        % sSliceArray.asSlice[0].sPosition.dSag/Cor/Tra
        ori = ['Sag'; 'Cor'; 'Tra']; % 1/2/3
        key = sprintf('sSliceArray.asSlice[%g].sPosition.d%s', i, ori(iSL,:));
        pos = asc_header(s, key);
        if ~isempty(pos)
            pos = [R(iSL, 1:2) pos] * [-dim(1:2)/2 1]'; % 2nd slice location
        end
    end
    
elseif isfield(s, 'LastFile') && isfield(s.LastFile, 'ImagePositionPatient')
    % s.LastFile works for most GE, Philips and all non-mosaic Siemens data
    R(1:3, 3) = (s.LastFile.ImagePositionPatient - R(1:3,4)) / (dim(3)-1);
    pixdim(3) = abs(R(iSL,3) / cosSL); % override slice thickness of dcm hdr
    return;
end

% May be useful for Philips dicom: use volume centre info
if isempty(pos) && isfield(s, 'Stack')
    ori = ['RL'; 'AP'; 'FH']; % x y z
    pos = tryGetField(s.Stack.Item_1, ['MRStackOffcentre' ori(iSL,:)]);
    if ~isempty(pos)
        pos = [R(iSL, 1:2) pos] * [-dim(1:2)/2 1]'; % mid-slice location
    end
end

% GE: LastScanLoc is always different from the slice in 1st file
if isempty(pos) && isfield(s, 'LastScanLoc')
    pos = s.LastScanLoc;
    if iSL<3, pos = -pos; end % LastScanLoc uses RAS convention!
end

if ~isempty(pos) % have SliceLocation for last/center slice
    flip = sign(pos-R(iSL,4)) ~= sign(R(iSL,3)); % same direction?
else % we do some guess work and warn user
    errorLog(['Please check whether slices are flipped: ' ProtocolName(s)]);
    pos1 = R(iSL, 3:4) * [dim(3)-1 1]'; % last SliceLocation based on current R
    pos2 = R(iSL, 3:4) * [1-dim(3) 1]'; % opposite slice direction
    % if pos1 is larger than the other dir, and is way outside head
    flip = all(abs(pos1) > [abs(pos2) 150]); % arbituary 150 mm
end
if flip, R(:,3) = -R(:,3); end % change to opposite direction
end

 

