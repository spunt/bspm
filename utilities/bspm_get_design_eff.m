function eff = bspm_get_design_eff(X,L,W)
% BSPM_GET_DESIGN_EFF
% 
% USAGE: eff = bspm_get_design_eff(X,L,W)
%
% ARGUMENTS
%   X = design (assumes no constant term present)
%   L = contrasts of interest (rows are contrast vectors)
%   W = weights for each contrast (default = equal weighting)
%

% -------------------- Copyright (C) 2014 --------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<2, mfile_showhelp; return; end
if nargin<3, W = ones(1,size(L,1)); end
if size(L,2)~=size(X,2), error('X and L must have same number of columns!'); return; end
if size(W,2)~=size(L,2), error('L and W must have same number of columsn!'); return; end

% mods
L(:,end+1) = 0;
L = L';
ncontrasts = size(L,2);

% compute
X(:,end+1) = 1;
for c = 1:ncontrasts
    eff(c) = 1/trace(L(:,c)'*pinv(X'*X)*L(:,c));
end
eff = eff*W';

% display
display(eff)
 
 
 
 
