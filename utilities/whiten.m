function [Z,W,U,D] = whiten(X)

%[Z,W,U,D] = whiten(X);
%Whiten whitens ("spheres") X, returning decorrelated transform of X, Z, the sphering matrix, W, the matrix of eigenvactors 
%of covariance of X, U, and diagonal matrix, D. X is a matrix with the column index, j, corresponding
% to the jth measurment. Whiten returns as many principle components as Rank(X).

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin==0, mfile_showhelp; return; end

sdim = min(size(X));
% sdim = rank(X);  

n = size(X,1);

%Centering
m = mean(X);                   
Xcnt = X-m(ones(1,size(X,1)),:);



C = (1/(n-1))*Xcnt'*Xcnt;       %Covariance matrix


%Principal components & whitening

[U,D] = svds(C,sdim);

W = D^(-.5)*U';

Z = (W*Xcnt')';