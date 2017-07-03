function [M,P,W,K,IM,IS,LH] = imm_update0(M,P,W,y,H,R)
% IMM_UPDATE0 - IMM update step
%
% Syntax:
%   [M,P,W,K,IM,IS,LH] = IMM_UPDATE0(M,P,W,Y,H,R)
%
% In:
%   M - Dx1xN mean state estimates
%   P - DxDxN state covariances
%   W - 1xN vector of mode weights
%   Y - Ex1 measurement vector.
%   H - ExD measurement matrix or ExDxN matrix
%       of mode conditional matrices
%   R - ExE measurement noise covariance or ExExN
%       matrix of mode conditional covariances
%
% Out:
%   X  - Updated state means
%   P  - Updated state covariances
%   W  - Updated mode weights
%   K  - Computed Kalman gains
%   IM - Means of predictive distributions of Y.
%   IS - Covariances or predictive means of Y.
%   LH - Predictive probability (likelihood) of measurement.
%   
% Description:
%   Perform Interactive Multiple Model (IMM) update step.
%   The IMM measurement model is:
%
%     y[k] = H_i*x_i[k] + r_i,  r_i ~ N(0,R_i)
%     i = 1,...,N
%
%   Predictive measurement distribution is defined as
%
%     p(y[k] |�y[1:k-1]) = sum w_i N(y[k] | IM_i[k], IS_i[k])
%
% See also:
%   IMM_PREDICT, KF_UPDATE
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


  %
  % Check which arguments are there
  %
  if nargin < 6
    error('Too few arguments');
  end

  %
  % Duplicate singleton measurement parameters
  %
  if size(H,3)==1
    H = repmat(H,[1 1 size(M,2)]);
  end
  if size(R,3)==1
    R = repmat(R,[1 1 size(M,2)]);
  end

  %
  % Mode conditional update steps
  %
  IM = zeros(size(y,1),size(M,2));
  IS = zeros(size(y,1),size(y,1),size(M,2));
  K  = zeros(size(M,1),size(y,1),size(M,2));
  LH = zeros(1,size(M,2));

  for i=1:size(M,2)
    IM(:,i)   = H(:,:,i) * M(:,i);
    IS(:,:,i) = R(:,:,i) + H(:,:,i) * P(:,:,i) * H(:,:,i)';
    K(:,:,i)  = P(:,:,i) * H(:,:,i)' / IS(:,:,i);
    M(:,i)    = M(:,i) + K(:,:,i) * (y - IM(:,i));
    P(:,:,i)  = P(:,:,i) - K(:,:,i) * IS(:,:,i) * K(:,:,i)';
    %LH(1,i)   = gauss_pdf(y,IM(:,i),IS(:,:,i));
    DX = y-IM(:,i);  
    E = 0.5*DX'*(IS(:,:,i)\DX);
%    E = E + 0.5 * size(M,1) * log(2*pi) + 0.5 * log(det(IS(:,:,i)));
    E = E + 0.5 * log(det(IS(:,:,i)));
    LH(1,i) = exp(-E);
  end

  %
  % Mode probability update
  %
  W  = W .* LH;
  LH = sum(W);
  W  = W ./ LH;
