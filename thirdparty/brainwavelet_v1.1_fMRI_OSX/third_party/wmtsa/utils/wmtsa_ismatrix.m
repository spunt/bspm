function [tf, sz] = wmtsa_ismatrix(x, type)
% wmtsa_ismatrix -- Determine if item is a matrix.
%
%****f* wmtsa.utils/wmtsa_ismatrix
%
% NAME
%   wmtsa_ismatrix -- Determine if item is a matrix.
%
% USAGE
%   [tf] = wmtsa_ismatrix(x, [type])
%
% INPUTS
%   * x           -- item to check (object).
%   * type        -- (optional) type of matrix (character string).
%                    Valid Values: 'point', 'truevector', 'truematrix'
%                    Default: No type specified.
%
% OUTPUTS
%   * tf         -- flag indicating whether item is a vector (Boolean).
%   * sz         -- size of x (vector of length 2).
%
% DESCRIPTION
%   By definition a matrix is a two dimensional array, which degenerates
%   to a vector or a point if it has, respectively, one or two singleton 
%   dimensions.  In MATLAB, points, vectors and matrices all have a 
%   dimensionality of two, (i.e. ndims(x) = 2).  A vector is the case of
%   where one dimension has length equal to one. A point is the degenerate 
%   case where both dimensions have lengths equal to one.
%   The function checks whether item
%   * item has 2 dimensions
%   If type is specified, it checks whether x is:
%   * a 'point' having two singleton dimensions.
%   * a 'truevector' having one singleton dimensions.
%   * a 'truematrix' having no singleton dimensions.
%
%   If the output argument 'sz' is specified, the size of the dimensions
%   of the item are returned.
%
% USAGE
%   tf = wmtsa_ismatrix(x)
%
%   tf = wmtsa_ismatrix(x, type)
%
%   [tf, sz] = wmtsa_ismatrix(x)
%
% ERRORS
%   WMTSA:InvalidNumArguments, WMTSA:InvalidArgumentValue
%
%
% EXAMPLE
%   x = ones(10,10);
%     % A matrix
%   tf = wmtsa_ismatrix(x)
%     % Result: tf =  1
%   tf = wmtsa_ismatrix(x, 'point')
%     % Result: tf =  0
%   tf  = wmtsa_ismatrix(x, 'truevector')
%     % Result: tf =  0
%   [tf, sz]  = wmtsa_ismatrix(x)
%     % Result: tf =  1, sz = [10 10]
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
%   2004-Apr-26
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

%   $Id: wmtsa_ismatrix.m 612 2005-10-28 21:42:24Z ccornish $


usage_str = ['Usage:  [tf, sz] = ', mfilename, ...
             '(x, [type])'];

%%  Check input arguments and set defaults.
error(nargerr(mfilename, nargin, [1:2], nargout, [0:2], 1, usage_str, 'struct'));


if(exist('type', 'var'))
  switch type
   case 'point'
    tf = wmtsa_isscalar(x);
   case 'truevector'
    tf = wmtsa_isvector(x, 'truevector');
   case 'truematrix'
    tf =  (ndims(x) == 2 && (size(x,1) ~= 1) && (size(x,2) ~= 1));
   otherwise
    error('WMTSA:InvalidArgumentValue', 'Invalid value for type');
  end  % switch
else
  tf = (ndims(x) == 2);
end % if


if (nargout > 1)
  sz = size(x);
end

return
