function [X, X0] = bspm_make_regressor(nvols, TR, TRbin, ons, dur, pm, TRons, TDtag)
% BSPM_MAKE_REGRESSOR
%
%   USAGE: X = bspm_make_regressor(nvols, TR, TRbin, ons, dur, pm, TRons, TDtag)
%   
%   ARGUMENTS
%       nvols = # of volumes
%       TR = TR (in secs)
%       TRbin = # of time bins per scan (for oversampling)
%       ons = onsets (in secs)
%       dur = durations (in secs) - if 0, model all as event
%       pm = parameters of interest (each column is a different parameter)
%       TRons = bin to use as onset
%       TDtag = include temporal derivative
%
%   OUTPUT
%       X       = convolved X matrix
%       X0      = un-convolved X matrix
%

% --------------------------------- Copyright (C) 2014 ---------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<8, TDtag = 0; end
if nargin<7, TRons = 1; end
if nargin<6, pm = []; end
if nargin<5, error('USAGE: X = bspm_make_regressor(nvols, TR, TRbin, ons, dur, pm, TRons, TDtag)'); end
if size(ons,1)==1 && size(ons,2)~=1, ons = ons'; end
if size(dur,1)==1 && size(dur,2)~=1, dur = dur'; end
if mod(TRbin, TR), TRbin = ceil(TRbin/TR)*TR; end
osrate = TR/(TR/TRbin);
xlength = nvols*osrate;
onso = fix((ons/TR)*osrate)+TRons;
duro = fix((dur/TR)*osrate);
if duro==0, duro = 1; end
if length(duro)==1, duro = repmat(duro,length(onso),1); end
hrf = spm_hrf(TR/TRbin);
xo = zeros(xlength,1);
for i = 1:length(onso)
    xo(onso(i):onso(i)+duro(i)-1) = 1;
end
xoc = conv(xo, hrf);
xoc = xoc(1:xlength);
X = xoc(1:osrate:xlength);
X0 = xo(1:osrate:xlength);
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
    end
    poc = poc(1:xlength,:);
    X(:,2:2+size(pm,2)-1) = poc(1:osrate:xlength,:);
    X0(:,2:2+size(pm,2)-1) = po(1:osrate:xlength,:);
end

if TDtag
    X1 = X;
    X01 = X0;
    clear xo xoc X X0 po cpm po poc 
    hrf = spm_hrf_td(TR/TRbin);
    xo = zeros(xlength,1);
    for i = 1:length(onso)
        xo(onso(i):onso(i)+duro(i)-1) = 1;
    end
    xoc = conv(xo, hrf);
    xoc = xoc(1:xlength);
    X = xoc(1:osrate:xlength);
    X0 = xo(1:osrate:xlength);
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
        end
        poc = poc(1:xlength,:);
        X(:,2:2+size(pm,2)-1) = poc(1:osrate:xlength,:);
        X0(:,2:2+size(pm,2)-1) = po(1:osrate:xlength,:);
    end
    X2 = X;
    X02 = X0;
    X = zeros(size([X1 X2]));
    X0 = zeros(size([X01 X02]));
    X(:,1:2:end) = X1;
    X(:,2:2:end) = X2;
    X0(:,1:2:end) = X01;
    X0(:,2:2:end) = X02;
    
end

 
 
 
 
