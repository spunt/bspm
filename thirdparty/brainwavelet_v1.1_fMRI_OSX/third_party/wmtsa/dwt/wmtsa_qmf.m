function [b] = wmtsa_qmf(a, inverse)
% wmtsa_qmf -- Calculate quadrature mirror filter (QMF).
%
%****f* wmtsa.dwt/wmtsa_qmf
%
% NAME
%   wmtsa_qmf -- Calculate quadrature mirror filter (QMF).
%
% USAGE
%   [b] = function_name(a, [inverse])
%
% INPUTS
%   * a           -- filter coefficients (vector).
%   * inverse     -- (optional) flag for calculating inverse QMF (Boolean).
%                    Default: inverse = 0 (FALSE).
%
% OUTPUTS
%    b            - QMF coefficients (vector).
%
% SIDE EFFECTS
%
%
% DESCRIPTION
%    wmtsa_qmf calculates the quadrature mirror filter (QMF) of
%    for the specified filter coefficients.  If a is a vector,
%    the QMF of the vector is calculated.  If a is a matrix or higher
%    order array, the QMF is calculated along the first dimension.
%
%   The inverse flag, if set, calculates the inverse QMF.  inverse
%   is a Boolean values specified as (1/0, y/n, T/F or true/false).
%
% EXAMPLE
%    % h is the QMF of g.
%    g = [0.7071067811865475 0.7071067811865475];
%    h = wmtsa_qmf(g);
%
%    % g is the inverse QMF of h.
%    h = [0.7071067811865475 -0.7071067811865475];
%    g = wmtsa_qmf(h, 1);
%
% ALGORITHM
%      g_l = (-1)^(l+1) * h_L-1-l
%      h_l = (-1)^l * g_L-1-l
%    See pages 75 of WMTSA for additional details.
%
% REFERENCES
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%   Time Series Analysis. Cambridge: Cambridge University Press.
%
% SEE ALSO
%   yn
%
% TOOLBOX
%   wmtsa/dwt
%
% CATEGORY
%   Filters:  Utilities
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2005-02-02
%
% COPYRIGHT
%   (c) 2005  Charles R. Cornish
%
% CREDITS
%
%
% REVISION
%   $Revision: 612 $
%
%***

%   $Id: wmtsa_qmf.m 612 2005-10-28 21:42:24Z ccornish $

usage_str = ['Usage:  [b] = ', mfilename, ...
             '(a, [inverse])'];


error(nargerr(mfilename, nargin, [1:2], nargout, [0:1], 1, usage_str, 'struct'));

if (~exist('inverse', 'var') || isempty(inverse))
  inverse = 0;
end

L = length(a);

if (wmtsa_isvector(a))
  b = flipvec(a);
else
  b = flipdim(a, 1)
end

if (yn(inverse))
  first = 1;
else
  first = 2;
end

b(first:2:end) = -b(first:2:end);
  
return
