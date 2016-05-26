function weights = bspm_conweights(design)
% BSPM_CONWEIGHTS
%
% Generate a matrix of non-redundant pairwise comparisons
% across conditions. Each row has a unique comparison.
%
% USAGE: weights = bspm_contrast_weights(design)
%
% ARGUMENTS
%   design, e.g., [2 3], [2 2]
%

% ------------------ Copyright (C) 2014 ------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<1, mfile_showhelp; return; end
ncell = prod(design);
nfact = length(design);
if nfact > 1
    design = sort(design);
    [p,q] = meshgrid(1:design(1), 1:design(2));
    pairs = [p(:) q(:)];
    [pairs, k, b] = unique(sort(pairs, 2),'rows');
    pairs(pairs(:,1)==pairs(:,2),:) = []; 
    wc  = zeros(size(pairs, 1), design(2));
    w1   = wc; 
    for ii = 1:size(pairs,1)
        w1(ii, pairs(ii,:)) = [1 -1];
    end
    w = [[w1 w1]; [w1 wc]; [wc w1]; [w1 w1*-1]];
    for ii = 1:design(2)
        w2(ii, find(q==ii)) = [1 -1];
    end
    w = [w; w2; sum(w2)];
    weights = w./repmat(sum(w==1,2), 1, ncell); 
else
    [p,q] = meshgrid(1:design);
    pairs = [p(:) q(:)];
    [pairs, k, b] = unique(sort(pairs, 2),'rows');
    pairs(pairs(:,1)==pairs(:,2),:) = []; 
    weights  = zeros(size(pairs, 1), design);
     for ii = 1:size(pairs,1)
        weights(ii, pairs(ii,:)) = [1 -1];
    end
end
end
 
 
 
