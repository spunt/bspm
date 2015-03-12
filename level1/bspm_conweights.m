function weights = bspm_conweights(ncond)
% BSPM_CONWEIGHTS
%
% Generate a matrix of non-redundant pairwise comparisons
% across conditions. Each row has a unique comparison.
%
% USAGE: weights = bspm_contrast_weights(ncond)
%
% ARGUMENTS
%   ncond = number of conditions
%

% ------------------ Copyright (C) 2014 ------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1, disp('USAGE: weights = bspm_conweights(ncond)'); return; end

[p,q] = meshgrid(1:ncond);
pairs = [p(:) q(:)];
pairs = unique(pairs,'rows');
useless = find(sum(pairs==fliplr(pairs),2));
for i = 1:ncond-1;
    tmp = find(pairs(:,2)==i);
    useless = [useless; tmp(i+1:end)];
end
pairs(useless,:) = [];
basevec = zeros(1,ncond);
for i = 1:size(pairs,1)
    tmp = basevec;
    tmp(pairs(i,:)) = [1 -1];
    weights(i,:) = tmp;
end
 
 
 
 
