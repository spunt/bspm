function [X, X0] = bspm_make_regressor(nscan, TR, ons, dur, varargin)
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
        'TRons',            8, ...
        'Derivatives',      0, ...
        'PMods',            []};
vals = setargs(def, varargin); 
if nargin<4, mfile_showhelp; disp(vals); return; end
if size(ons,1)==1 & size(ons,2)~=1, ons = ons'; end
if size(dur,1)==1 & size(dur,2)~=1, dur = dur'; end
if size(PMods, 1)==1, PMods = PMods'; end

% | Time Units
dt      = TR/TRbin;

% | Basis Functions
optDerivative = {'hrf',...
                'hrf (with time derivative)',...
                'hrf (with time and dispersion derivatives)'};
bf = struct( ...
        'UNITS',        'secs', ...
        'T',            TRbin,  ...
        'T0',           TRons,  ...
        'name',         optDerivative{Derivatives + 1},  ...
        'Volterra',     1,      ...
        'dt',           dt);
bf  = spm_get_bf(bf); 

% | Peri-stimulus times {seconds}
pst   = [0:(nscan-1)]*TRbin*dt - min(ons);
for j = 1:length(ons)
    w      = [0:(nscan-1)]*TRbin*dt - ons(j);
    v      = find(w >= 0);
    pst(v) = w(v);
end
u     = ons.^0;

% | Parametric Modulators
if ~isempty(PMods)
    PMods   = scalepms(PMods);
    PMods   = PMods - repmat(mean(PMods), size(PMods, 1), 1); 
    u       = [u PMods]; 
end

% | And scale so sum(u*dt) = number of events, if event-related
if ~any(dur), u     = u/dt; end

% | Create stimulus functions (32 bin offset)
ton       = round(ons/dt) + 33;               % onsets
tof       = round(dur/dt) + ton + 1;          % offset
X0        = sparse((nscan*TRbin + 128),size(u,2));
ton       = max(ton,1);
tof       = max(tof,1);
for j = 1:length(ton)
    if size(X0,1) > ton(j)
        X0(ton(j),:) = X0(ton(j),:) + u(j,:);
    end
    if size(X0,1) > tof(j)
        X0(tof(j),:) = X0(tof(j),:) - u(j,:);
    end
end
X0        = cumsum(X0);                         % integrate
X0        = X0(1:(nscan*TRbin + 32),:);         % stimulus

% | Convolve stimulus functions with basis functions
X     = [];
for k = 1:size(X0, 2)
    for p = 1:size(bf.bf, 2)
        x = X0(:,k);
        d = 1:length(x);
        x = conv(full(x), bf.bf(:,p));
        x = x(d);
        X = [X x];
    end
end
  
% | Resample regressors at acquisition times (32 bin offset)
X0 = full(X0); 
if ~isempty(X) 
    X = X((0:(nscan - 1))*TRbin + TRons + 32,:);
end


end



% 
% osrate  = TR/(TR/TRbin);
% xlength = nscan*osrate;
% onso    = fix((onsets/TR)*osrate)+TRons;
% duro    = fix((durations/TR)*osrate);
% if duro==0, duro = 1; end
% if length(duro)==1, duro = repmat(duro,length(onso),1); end
% hrf     = spmhrf(TR/TRbin, TRbin, TempDeriv);
% xo      = zeros(xlength,1);
% for i = 1:length(onso)
%     xo(onso(i):onso(i)+duro(i)-1) = 1;
% end
% xoc     = conv(xo, hrf(:,1));
% xoc     = xoc(1:xlength);
% if TempDeriv
%     xoc2    = conv(xo, hrf(:,2)); 
%     xoc     = [xoc xoc2(1:xlength)];
% end
% X       = xoc(1:osrate:xlength, :);
% X0      = xo(1:osrate:xlength);
% % | PARAMETRIC MODULATION
% if ~isempty(PMods)
%     PMods = scalepms(PMods);
%     po = repmat(xo,1,size(PMods,2));
%     poc = []; 
%     for p = 1:size(PMods,2)
%         cpm = PMods(:,p);
%         cpm = cpm - mean(cpm);
%         for i = 1:length(onso)
%             po(onso(i):onso(i)+duro(i)-1,p) = 1*cpm(i);
%         end
%         cpoc = conv(po(:,p),hrf(:,1));
%         cpoc = cpoc(1:xlength); 
%         if TempDeriv
%             p2 = conv(po(:,p),hrf(:,2));
%             cpoc = [cpoc p2(1:xlength)];
%         end
%         poc = [poc cpoc]; 
%     end
%     X   = [X poc(1:osrate:xlength,:)];
%     X0  = [X0 po(1:osrate:xlength,:)];
% end

% | SUBFUNCTIONS
function out    = scalepms(in)
if nargin<1, error('out = scalepms(in)'); end
nrow = size(in, 1);
out = (in - repmat(min(in), nrow, 1)) ./ (repmat(max(in), nrow, 1) - repmat(min(in), nrow, 1)); 
end