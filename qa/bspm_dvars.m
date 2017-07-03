function [dvars, framewise] = bspm_dvars(epi, varargin)
% BSPM_DVARS
%
% USAGE: TS = bspm_dvars(epi, cutoff, maskfn, nosave)
%       
%

% ------------------------ Copyright (C) 2014 ------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

%% CHECK INPUTS
def = { 'dvars_thresh',     2.5,  ...
        'framewise_thresh', 0.5,  ...
        'include_rp',       1,    ...
        'makeplot',         0,    ...
        'prefix',           'badscan', ...
        'maskfile',         []};
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end

% | PATH TO BRAMILA TOOLS
cfg.bramilapath = fullfile(getenv('HOME'), 'Github', 'bspm', 'thirdparty', 'bramila'); % bramila toolbox path
% cfg.StdTemplate=bspm_template; % 2mm MNI template from FSL
% cfg.TR = 2.5; % TR from scanning protocol, used in bramila
bramiladir = fullfile(getenv('HOME'), 'Github', 'bspm', 'thirdparty', 'bramila');
addgenpath(bramiladir, 'excl',  '_archive_'); 

% | get data and apply implicit masking
if ischar(epi), epi = cellstr(epi); end
fprintf('\n | - Loading Data');
maskthresh  = 0.8; 
[v,h]       = bnii_read(epi, 'reshape', 1);
[nvox nvol] = size(v);
if ~isempty(maskfile)
%     v(v < repmat(mean(v)*maskthresh, nvox, 1)) = NaN;
    m = bspm_reslice(maskfile, strcat(epi, ',1'), 0, 1);
    v(m(:)==0,:) = 0; 
end
cfg.vol     = reshape(v, h.dim(2:5));

% | create output filename
outfile     = sprintf('%s_dvars%dframewise%d_%s.txt', prefix, dvars_thresh*100, framewise_thresh*100, strtrim(datestr(now,'mmm_DD_YYYY')));
epidir      = fileparts(epi{1});

% | motion parameters
rp                  = load(char(files(fullfile(epidir, 'rp_*txt'))));
rp(:,4:6)           = rp(:,4:6)*57.3;
framewise           = zeros(nvol,1);
framewise(2:end)    = max(abs(diff(rp)), [], 2);

% | use BRAMILA tools to compute dvars and framewise displacement
fprintf('\n | - Identifying Scans with DVARS > %2.2f and framewise displacement > %2.2f\n', dvars_thresh, framewise_thresh);
cfg.plot            = makeplot;
dvars               = bramila_dvars(cfg);

    
end
% =========================================================================
% SUBFUNCTIONS
% =========================================================================
function zout = oneoutzscore(in)
% leave-one-out zscore
if size(in,1)==1, in=in'; end
n = length(in); 
in = repmat(in, 1, n); 
theoneout = in(logical(eye(n)))';  
theleftin = reshape(in(logical(~eye(n))),n-1,n);
zout = (theoneout-mean(theleftin))./std(theleftin);
zout = zout';
end