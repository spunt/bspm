function [M,P,D] = rts_smooth(M,P,A,Q)
% RTS_SMOOTH - Rauch-Tung-Striebel smoother
%
% Syntax:
%   [M,P,S] = RTS_SMOOTH(M,P,A,Q)
%
% In:
%   M - NxK matrix of K mean estimates from Kalman filter
%   P - NxNxK matrix of K state covariances from Kalman Filter
%   A - NxN state transition matrix or NxNxK matrix of K state
%       transition matrices for each step.
%   Q - NxN process noise covariance matrix or NxNxK matrix
%       of K state process noise covariance matrices for each step.
%
% Out:
%   M - Smoothed state mean sequence
%   P - Smoothed state covariance sequence
%   D - Smoother gain sequence
%   
% Description:
%   Rauch-Tung-Striebel smoother algorithm. Calculate "smoothed"
%   sequence from given Kalman filter output sequence
%   by conditioning all steps to all measurements.
%
% Example:
%   m = m0;
%   P = P0;
%   MM = zeros(size(m,1),size(Y,2));
%   PP = zeros(size(m,1),size(m,1),size(Y,2));
%   for k=1:size(Y,2)
%     [m,P] = kf_predict(m,P,A,Q);
%     [m,P] = kf_update(m,P,Y(:,k),H,R);
%     MM(:,k) = m;
%     PP(:,:,k) = P;
%   end
%   [SM,SP] = rts_smooth(MM,PP,A,Q);
%
% See also:
%   KF_PREDICT, KF_UPDATE

% Copyright (C) 2003-2006 Simo Särkkä
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
  % Check which arguments are there
  %
  if nargin < 4
    error('Too few arguments');
  end

  %
  % Extend A and Q if they are NxN matrices
  %
  if size(A,3)==1
    A = repmat(A,[1 1 size(M,2)]);
  end
  if size(Q,3)==1
    Q = repmat(Q,[1 1 size(M,2)]);
  end

  %
  % Run the smoother
  %
  D = zeros(size(M,1),size(M,1),size(M,2));
  for k=(size(M,2)-1):-1:1
    P_pred   = A(:,:,k) * P(:,:,k) * A(:,:,k)' + Q(:,:,k);
    D(:,:,k) = P(:,:,k) * A(:,:,k)' / P_pred;
    M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - A(:,:,k) * M(:,k));
    P(:,:,k) = P(:,:,k) + D(:,:,k) * (P(:,:,k+1) - P_pred) * D(:,:,k)';
  end

