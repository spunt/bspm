function [tf] = wmtsa_isscalar(x)
% wmtsa_isscalar -- Determine if item is a point.
%
%****f* wmtsa.utils/wmtsa_isscalar
%
% NAME
%   wmtsa_isscalar -- Determine if item is a scalar.
%
% USAGE
%   [tf] = isscalar(x)
%
% INPUTS
%   * x          -- item to check (object).
%
% OUTPUTS
%   * tf          -- flag indicating whether item is a point (numeric Boolean).
%
% DESCRIPTION
%   In MATLAB, points, vectors and matrices all have a dimensionality of two,
%   (i.e. ndims(x) = 2).  A point is the degenerate case where an array, has
%   size of one in all dimensions, i.e. the array is singleton in all dimensions.
%   Function checks whether item is a point by checking that its length
%   in all dimensions is one.
%
% ERRORS
%   WMTSA:InvalidNumArguments
%
% TOOLBOX
%   wmtsa/utils
%
% CATEGORY
%   WMTSA Utilities
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2005-Feb-03
%
% COPYRIGHT
%
%
% CREDITS
%
%
% REVISION
%   $Revision: 612 $
%
%***

%   $Id: wmtsa_isscalar.m 612 2005-10-28 21:42:24Z ccornish $

usage_str = ['Usage:  [tf] = ', mfilename, '(x)'];

%%  Check input arguments and set defaults.
error(nargerr(mfilename, nargin, [1:2], nargout, [0:2], 1, usage_str, 'struct'));

sz = size(x);

if ( (sum(sz) / length(sz)) == 1)
  tf = 1;
else
  tf = 0;
end

return
