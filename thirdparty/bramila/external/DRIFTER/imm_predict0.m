function [M,P,W] = imm_predict0(M,P,W,T,A,Q,B,U)
% IMM_PREDICT0 - Perform IMM prediction step
%
% Syntax:
%   [M,P,W] = IMM_PREDICT0(M,P,W,T,A,Q,B,U)
%
% Author:
%   Simo S�rkk�, 2003
%
% In:
%   M - DxN mean state estimates of previous step
%   P - DxDxN state covariances from previous step
%   W - 1xN weights (posterior probabilities of modes)
%   T - Mode transition probabilities T(i,j) = p(i|j)
%   A - DxD transition matrix of discrete model or DxDxN
%       matrix containing separate matrices for each mode
%       (optional, default identity)
%   Q - DxD process noise covariance of discrete model or DxDxN
%       matrix containing separate covariances for each mode
%       (optional, default zero matrix)
%   B - DxS input effect matrix or DxSxN matrix containing separate
%       input effect matrices for each mode (optional, default identity)
%   U - Sx1 constant input or SxN matrix containing separate
%       inputs for each mode                (optional, default empty)
%
% Out:
%   X  - Predicted state means as DxN matrix
%   P  - Predicted state covariance as DxDxN matrix
%   W  - Predicted weights as 1xN matrix
%   
% Description:
%   Perform Interacting Multiple Model (IMM)
%   prediction step. The model is
%
%     p(i|j) = T_{ij}
%     x_i[k] = A_i*x_i[k-1] + B*u[k-1] + q_i,  q_i ~ N(0,Q_i)
%     i = 1,...,N
%
% See also:
%   IMM_UPDATE, KF_PREDICT, LTI_DISC, EKF_PREDICT, EKF_UPDATE
%
% History:
%   18.02.2003  The first official version.

% Copyright (C) 2003 Simo Särkkä
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

%%
  %
  % Check arguments
  %
  if nargin < 5
    A = [];
  end
  if nargin < 6
    Q = [];
  end
  if nargin < 7
    B = [];
  end
  if nargin < 8
    U = [];
  end
  
  %
  % Apply defaults
  %
  if isempty(A)
    A = eye(size(M,1));
  end
  if isempty(Q)
    Q = zeros(size(M,1));
  end
  if isempty(B) && ~isempty(U)
    B = eye(size(M,1),size(U,1));
  end
  if size(A,3)==1
    A = repmat(A,[1 1 size(M,2)]);
  end
  if size(Q,3)==1
    Q = repmat(Q,[1 1 size(M,2)]);
  end

  %
  % Calculate mixing probabilities
  %
  mu_ij = T.*repmat(W',1,size(W,2));
  CJ = sum(mu_ij,1);
  % CJ(find(CJ==0)) = eps;
  CJ(CJ==0) = eps;   % Equivalent to find() above
  mu_ij = mu_ij ./ repmat(CJ,size(mu_ij,1),1);
  
  %
  % Mixing
  %
  M0 = zeros(size(M));
  P0 = zeros(size(P));

  for j=1:size(M,2)
    for i=1:size(M,2)
        if mu_ij(i,j) > eps %!!!!
            M0(:,j) = M0(:,j) + M(:,i)*mu_ij(i,j);
        end
    end
    for i=1:size(M,2)
        if mu_ij(i,j) > eps %!!!!
            C = P(:,:,i) + (M(:,i) - M0(:,j))*(M(:,i) - M0(:,j))';
            P0(:,:,j) = P0(:,:,j) + mu_ij(i,j) * C;
        end
    end
  end

  %
  % Mode matched prediction
  %
  for i=1:size(M,2)
    if isempty(U)
      M(:,i) = A(:,:,i) * M0(:,i);
    else
      M(:,i) = A(:,:,i) * M0(:,i) + B(:,:,i) * U(:,i);
    end
    P(:,:,i) = A(:,:,i) * P0(:,:,i) * A(:,:,i)' + Q(:,:,i);
  end

  %
  % Weights are actually the
  % normalization constants
  %
  W = CJ ./ sum(CJ);


