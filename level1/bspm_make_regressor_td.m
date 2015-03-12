function X = bspm_make_regressor_td(nvols, TR, TRbin, ons, dur, pm)
% BSPM_MAKE_REGRESSOR
%
%   USAGE: X = bspm_make_regressor(nvols, TR, TRbin, ons, dur, pm)
%   
%   ARGUMENTS
%       nvols = # of volumes
%       TR = TR (in secs)
%       TRbin = # of time bins per scan (for oversampling)
%       ons = onsets (in secs)
%       dur = durations (in secs) - if 0, model all as event
%       pm = parameters of interest (each column is a different parameter)
%
%   OUTPUT
%       X = X matrix
%

% -------------------------- Copyright (C) 2014 --------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin==5, pm = []; end
if nargin<5, display('USAGE: X = bspm_make_regressor(nvols, TR, TRbin, ons, dur, pm)'); end
if size(ons,1)==1 && size(ons,2)~=1, ons = ons'; end
if size(dur,1)==1 && size(dur,2)~=1, dur = dur'; end
if mod(TRbin, TR), TRbin = ceil(TRbin/TR)*TR; end
osrate = TR/(TR/TRbin);
xlength = nvols*osrate;
onso = fix((ons/TR)*osrate)+1;
duro = fix((dur/TR)*osrate)+1;
if duro==0, duro = 1; end
if length(duro)==1
    duro = repmat(duro,length(onso),1);
end
hrf = spm_hrf(TR/TRbin);
hrf_td = spm_hrf_td(TR/TRbin);
xo = zeros(xlength,1);
for i = 1:length(onso)
    xo(onso(i):onso(i)+duro(i)-1) = 1;
end
xoc = conv(xo, hrf);
xoc_td = conv(xo, hrf_td);
xoc = xoc(1:xlength);
xoc_td = xoc_td(1:xlength);
X1 = xoc(1:osrate:xlength);
X2 = xoc_td(1:osrate:xlength);
if ~isempty(pm)
    pm = bob_scalematrix(pm);
    po = repmat(xo,1,size(pm,2));
    for p = 1:size(pm,2)
        cpm = pm(:,p);
        cpm = cpm - mean(cpm);
        for i = 1:length(onso)
            po(onso(i):onso(i)+duro(i)-1,p) = 1*cpm(i);
        end
        poc(:,p) = conv(po(:,p),hrf);
        poc_td(:,p) = conv(po(:,p),hrf_td);
    end
    poc = poc(1:xlength,:);
    poc_td = poc_td(1:xlength,:);
    X1(:,2:2+size(pm,2)-1) = poc(1:osrate:xlength,:);
    X2(:,2:2+size(pm,2)-1) = poc_td(1:osrate:xlength,:);
end
X = zeros(size(X1,1),2*size(X1,2));
X(:,1:2:end) = X1;
X(:,2:2:end) = X2;
 
 
 
 
