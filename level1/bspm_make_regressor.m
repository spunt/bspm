function [X, X0] = bspm_make_regressor(nscan, TR, onsets, durations, varargin)
% BSPM_MAKE_REGRESSOR Build a single timeseries regressor convolved with HRF 
%
%   USAGE: [X, X0] = bspm_make_regressor(nscan, TR, onsets, durations, varargin)
%
%   NECESSARY ARGUMENTS
%       nscan       - N scans (timepoints)
%       TR          - Repetition time (secs)
%       onsets      - Event onsets (secs)
%       durations   - Event durations (secs; enter 0 to model as Events)
%      
%   OPTIONAL ARGUMENTS (Name-Value Pairs):
%       TRbin       - N Bins Per Scan (for upsampling)
%       TRons       - Bin to Use as Onset
%       TempDeriv   - Tag to include Temporal Derivative
%       PMods       - Parameters: rows are events, columsn are parameters
%
%   OUTPUTS
%       X           - convolved X matrix
%       X0          - un-convolved X matrix
%

% --------------------------------- Copyright (C) 2014 ---------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
def = { 'TRbin',            16, ...
        'TRons',            1, ...
        'TempDeriv',        0, ...
        'PMods',            []};
bspm_setdefaults(def, varargin); 
if nargin<4, error('USAGE: [X, X0] = bspm_make_regressor(nscan, TR, onsets, durations, varargin)'); end
if size(onsets,1)==1 & size(onsets,2)~=1, onsets = onsets'; end
if size(durations,1)==1 & size(durations,2)~=1, durations = durations'; end
% if mod(TRbin, TR), TRbin = ceil(TRbin/TR)*TR; end
osrate  = TR/(TR/TRbin);
xlength = nscan*osrate;
onso    = fix((onsets/TR)*osrate)+TRons;
duro    = fix((durations/TR)*osrate);
if duro==0, duro = 1; end
if length(duro)==1, duro = repmat(duro,length(onso),1); end
hrf     = spmhrf(TR/TRbin, TRbin, TempDeriv);
xo      = zeros(xlength,1);
for i = 1:length(onso)
    xo(onso(i):onso(i)+duro(i)-1) = 1;
end
xoc     = conv(xo, hrf(:,1));
xoc     = xoc(1:xlength);
if TempDeriv
    xoc2    = conv(xo, hrf(:,2)); 
    xoc     = [xoc xoc2(1:xlength)];
end
X       = xoc(1:osrate:xlength, :);
X0      = xo(1:osrate:xlength);
% | PARAMETRIC MODULATION
if ~isempty(PMods)
    PMods = scalepms(PMods);
    po = repmat(xo,1,size(PMods,2));
    poc = []; 
    for p = 1:size(PMods,2)
        cpm = PMods(:,p);
        cpm = cpm - mean(cpm);
        for i = 1:length(onso)
            po(onso(i):onso(i)+duro(i)-1,p) = 1*cpm(i);
        end
        cpoc = conv(po(:,p),hrf(:,1));
        cpoc = cpoc(1:xlength); 
        if TempDeriv
            p2 = conv(po(:,p),hrf(:,2));
            cpoc = [cpoc p2(1:xlength)];
        end
        poc = [poc cpoc]; 
    end
    X   = [X poc(1:osrate:xlength,:)];
    X0  = [X0 po(1:osrate:xlength,:)];
end
end
% | SUBFUNCTIONS
function out    = scalepms(in)
if nargin<1, error('out = scalepms(in)'); end
nrow = size(in, 1);
out = (in - repmat(min(in), nrow, 1)) ./ (repmat(max(in), nrow, 1) - repmat(min(in), nrow, 1)); 
end
function bf     = spmhrf(RT, TRbin, TempDeriv)
    if nargin<3, TempDeriv = 0; end
    p   = [6 16 1 1 6 0 32];
    dt  = RT/TRbin;
    u   = (0:ceil(p(7)/dt)) - p(6)/dt;
    bf = spm_Gpdf(u,p(1)/p(3),dt/p(3)) - spm_Gpdf(u,p(2)/p(4),dt/p(4))/p(5);
    bf = bf((0:floor(p(7)/RT))*TRbin + 1);
    bf = bf'/sum(bf);
    if TempDeriv
        dp      = 1;
        p(6)    = p(6) + dp;
        u       = (0:ceil(p(7)/dt)) - p(6)/dt;
        D       = spm_Gpdf(u,p(1)/p(3),dt/p(3)) - spm_Gpdf(u,p(2)/p(4),dt/p(4))/p(5);
        D       = D((0:floor(p(7)/RT))*TRbin + 1);
        D       = D'/sum(D);
        D       = (bf(:,1) - D)/dp;
        bf      = [bf D(:)];
    end
end