function [im, rec_data] = read_parrec(filein,varargin)
%
%[im rec_data] = readpar(parfile);
%[im rec_data] = readpar(parfile,phase);
%[im rec_data] = readpar(parfile,phase,field);
%read scan parameters from exported PAR files (Philips research tools)
%
% output :
% - par   : scan parameters and image scaling factors -
%           float_data = (int_data * RS + RI) / (RS*SS)
%
%

phase = -1;
field = 1;
if (nargin>1)
    phase=varargin{1};
end
if (nargin>2)
    field=varargin{2};
end



par_data = readpar(filein);
if par_data.phases ==1
    phase=1;
end

index=strfind(filein,'.PAR');
rec_data_file=[filein(1:index-1) '.REC'];

% Initialise output

info = par_data;
rec_data = readrec(rec_data_file,info);
data = (rec_data(:,:,:,:,1,1,field));


if (field==3)
    % set the max and min values!
    maxval = info.phaseEncZ + info.phaseEncX + info.phaseEncY;
    
    data = data*2*maxval/(2^(12))-maxval;
    
else
    data=data*512/max(data(:));
end







AFRtoLPS = [0 0 1; 1 0 0; 0 1 0];
% ORIENTATION -------------------------------------------------

magicmatrix = [-1 0 0; 0 0 1; 0 -1 0];
TRA = [0 1 0; -1 0 0; 0 0 -1];
SAG = [-1 0 0; 0 0 1; 0 -1 0];
COR = [0 1 0; 0 0 1; -1 0 0];

switch (str2num(par_data.sliceorient))
    case 1
        Torientation = TRA;
    case 2
        Torientation = SAG;
    case 3
        Torientation = COR;
    otherwise
        Torientation = eye(3);
end


ap = par_data.angAP * pi / 180.0;
fh = par_data.angFH * pi / 180.0;
rl = par_data.angRL * pi / 180.0;

Tap = [1 0 0; 0 cos(ap) -sin(ap); 0 sin(ap) cos(ap)];
Tfh = [cos(fh) 0 sin(fh); 0 1 0; -sin(fh) 0 cos(fh)];
Trl = [cos(rl) -sin(rl) 0; sin(rl) cos(rl) 0; 0 0 1];

orientation_matrix = AFRtoLPS * Trl * Tap * Tfh * magicmatrix' * Torientation'; % c++ version
%orientation_matrix =  Trl * Tap * Tfh * Torientation'; % c++ version
%orientation_matrix =  Trl * Tap * Tfh * magicmatrix' * Torientation';
%orientation_matrix = inv(Torientation * magicmatrix * inv (Trl * Tap * Tfh)); %matlab version
%orientation_matrix = Trl * Tap * Tfh *magicmatrix' * Torientation'  ; %matlab version
im.orientation = orientation_matrix;
%orientation_matrix_str = num2str(reshape(orientation_matrix,1,9));


s = size(data)';

% ORIGIN ----------------------------------------------------

midoffset =    info.vox(1:3)' .* (s(1:3)- ones(3,1)) ./2.0;
midoffset = orientation_matrix * midoffset;

%origin = [par_data1.offRL par_data1.offAP par_data1.offFH ]';
origin = [par_data.offAP par_data.offFH par_data.offRL ]'; % like this in par files


origin = AFRtoLPS * origin;

origin = origin - midoffset;

%-------------------------------

if phase>-1
    im.data = squeeze(data(:,:,:,phase));
    spacing = [info.vox]';
else
    for i = 1:info.phases
        timings(i) = mean(info.slicetime(info.ph==i));
    end
    temporal_spacing = mean(timings(2:end)-timings(1:end-1))/1000;
    spacing = [info.vox temporal_spacing]';
    origin = [origin ; mean(info.delay)];
end
M = eye(numel(s));
M(1:3,1:3)=orientation_matrix;


im = ImageType(s,origin,spacing,M);
im.data = data;



end




