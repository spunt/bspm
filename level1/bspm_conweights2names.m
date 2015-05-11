function name = bspm_conweights2names(weights, covnames)
% BSPM_CONWEIGHTS2NAMES Make constrast names from weights
%
%  USAGE: name = bspm_conweights2names(weights, covnames)
%
%  INPUTS
%	weights:    1xN matrix of weights
%   covnames:   1xN cell array of covariate names
% __________________________________________________________________________ 
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-02-24
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 2, disp('USAGE: names = bspm_conweights2names(weights, covnames)'); return; end
if size(weights, 2)~=length(covnames), error('WEIGHTS and COVNAMES must be the same length'); end 
posidx = []; posidx = find(weights>0);
negidx = []; negidx = find(weights<0);
if isempty(negidx)
    name = strcat(covnames{posidx});
elseif isempty(posidx)
    name = ['INV_' strcat(covnames{negidx})];
else
    name = [strcat(covnames{posidx}) '_-_' strcat(covnames{negidx})];
end