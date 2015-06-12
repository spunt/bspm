function LOCMAX = bspm_localmaxima(FNAME, varargin)
% BSPM_LOCALMAXIMA Get local maxima from T-image
%
%  USAGE: LOCMAX = bspm_localmaxima(FNAME, varargin)
%
% FNAME
%
%   Dis    [20]     - min distance among clusters {mm} [8]
%   Num    [3]      - number of maxima per cluster [3]
%   pval   [.001]   - uncorrected P value
%   k      [5]      - extent threshold {voxels}
%   STAT   ['T']    - distribution {Z, T, X or F} 
%   SIGN   ['+']    - sign of effect {+, -, or +/-}
%
% OUTPUT
%   Col 1: t-value
%   Col 2: cluster extent
%   Col 3: x-coordinates
%   Col 4: y-coordinates
%   Col 5: z-coordiantes
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-05-15
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
def = { ... 
	'Dis',		20,     ...
	'Num',		3,      ...
	'STAT',		'T',	...
	'k',		5,      ...
	'pval',		.001,	...
    'SIGN',     '+'     ...
	};
vals = setargs(def, varargin);
if nargin < 1, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end

% | Load Data
if iscell(FNAME), FNAME = char(FNAME); end
oh              = spm_vol(FNAME);
od              = spm_read_vols(oh);
od(isnan(od))   = 0;

% | PARAMS
% S       = sum(od(:)~=0);    % - search Volume {voxels}
DIM     = oh.dim';          % - image dimensions {voxels}
[X,Y,Z] = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
XYZ = [X(:)';Y(:)';Z(:)'];  % - location of voxels {voxel coords}
M       = oh.mat;         	% - voxels --> mm matrix
% VOX = abs(diag(M(:,1:3)));  % - voxel dimensions {mm}
RCP     = XYZ; RCP(4,:)=1;   
XYZmm   = M(1:3,:)*RCP;     % - location of voxels {mm}

% | DF
tmp     = oh.descrip;
idx{1}  = regexp(tmp,'[','ONCE');
idx{2}  = regexp(tmp,']','ONCE');
df      = str2num(tmp(idx{1}+1:idx{2}-1));
u       = spm_invTcdf(1-pval, df); 

% | raw data to XYZ
switch SIGN
    case '+'
        LOCMAX    = getmaxima(od, u, k, XYZ, M, XYZmm, Dis, Num);
    case '-'
        LOCMAX    = getmaxima(od*-1, u, k, XYZ, M, XYZmm, Dis, Num);
    otherwise
        POS     = getmaxima(od, u, k, XYZ, M, XYZmm, Dis, Num);
        NEG     = getmaxima(od*-1, u, k, XYZ, M, XYZmm, Dis, Num);
        NEG(:,1) = NEG(:,1)*-1; 
        LOCMAX    = [POS; NEG]; 
end
end
% - SUBFUNC
function PEAK = getmaxima(img, u, k, XYZ, M, XYZmm, Dis, Num)
supra = (img(:) > u)';
if sum(supra)
    tmp             = spm_clusters(XYZ(:, supra));
    clbin           = repmat(1:max(tmp), length(tmp), 1)==repmat(tmp', 1, max(tmp));
    C(supra)        = sum(repmat(sum(clbin), size(tmp, 2), 1) .* clbin, 2)';
    I(1,supra)      = tmp;
end
C(C < k)            = 0; 
idx                 = find(C > 0);
Z                   = img(idx);
XYZ                 = XYZ(:,idx);
XYZmm               = XYZmm(:,idx);
[N,Z,XYZ,A,L]       = spm_max(Z,XYZ);
XYZmm               = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];

%-Cluster and local maxima p-values & statistics
%----------------------------------------------------------------------
npeak = 0; 
while numel(find(isfinite(Z)))

    %-Find largest remaining local maximum
    %------------------------------------------------------------------
    [U,i]   = max(Z);            %-largest maxima
    j       = find(A == A(i));   %-maxima in cluster
    npeak   = npeak + 1;         %-number of peaks
    extent  = N(i); 
    PEAK(npeak,:) = [U extent XYZmm(:,i)']; %-assign peak
    
    %-Print Num secondary maxima (> Dis mm apart)
    %------------------------------------------------------------------
    [l,q] = sort(-Z(j));                              % sort on Z value
    D     = i;
    for i = 1:length(q)
        d = j(q(i));
        if min(sqrt(sum((XYZmm(:,D)-repmat(XYZmm(:,d),1,size(D,2))).^2)))>Dis
            if length(D) < Num
                D          = [D d];
                npeak   = npeak + 1;         %-number of peaks
                PEAK(npeak,:) = [Z(d) extent XYZmm(:,d)']; %-assign peak
            end
        end
    end
    Z(j) = NaN;     % Set local maxima to NaN
end
end




