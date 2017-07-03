function [m,C] = imm_estimate0(M,P,W)
% IMM_ESTIMATE0 - Calculate mean and covariance of IMM mixture
%
% Syntax:
%   [m,C] = IMM_ESTIMATE0(M,P,W)
%
% Author:
%   Simo S�rkk�, 2003
%
% In:
%   M - Dx1xN mean state estimates
%   P - DxDxN state covariances
%   W - 1xN vector of mode weights
%
% Out:
%   m - State estimate mean
%   C - State estimate covariance
%   
% Description:
%   Calculate mean and covariance of IMM mixture.
%   This is just a helper routine to ease the interpretation
%   of mixture Gaussian output of IMM filter.
%
% See also:
%   IMM_PREDICT, IMM_UPDATE
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
  m = zeros(size(M,1),1);
  for i=1:size(M,2)
    m = m + W(i) * M(:,i);
  end
  C = zeros(size(M,1),size(M,1));
  for i=1:size(M,2)
    C = C + W(i) * (P(:,:,i) + (M(:,i) - m) * (M(:,i) - m)');
  end
