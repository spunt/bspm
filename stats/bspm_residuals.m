function res = bspm_residuals(Y,X)
% BOB_RESIDUALS
%
%   USAGE: res = bspm_residuals(Y,X)
% 
%   ARGUMENTS
%       Y:  n-by-m matrix of observations, where different observations
%            appear in different columns
%       X:  n-by-m design matrix (do not include constant)
%
% ------------------------------------------------------------
if nargin<2, display('USAGE: res = bspm_residuals(Y,X)'); return; end
X = [X ones(size(X(:,1)))]; % add constant
res = Y - X*(X\Y); % compute residuals
res = res + repmat(nanmean(Y), size(X,1), 1); % add mean to residuals

 
 
 
 
