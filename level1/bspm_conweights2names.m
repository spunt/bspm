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
if nargin < 2, mfile_showhelp; return; end
if size(weights, 2)~=length(covnames), error('WEIGHTS and COVNAMES must be the same length'); end
ncon = size(weights, 1);
name = cell(ncon, 1); 
for c = 1:ncon

    posidx = []; posidx = find(weights(c,:)>0);
    negidx = []; negidx = find(weights(c,:)<0);

    if isempty(negidx)
        name{c} = strcat(covnames{posidx}, '_-_Baseline');
    elseif isempty(posidx)
        name{c} = ['Baseline_-_' strcat(covnames{negidx})];
    else
        name{c} = [strcat(covnames{posidx}) '_-_' strcat(covnames{negidx})];
    end

end
if ncon==1, name = char(name); end