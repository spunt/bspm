function [y] = flipvec(x)
% flipvec -- Flip a vector.
%
%****f* wmtsa.utils/flipvec.m
%
% NAME
%   flipvec -- Flip a vector.
%
% SYNOPSIS
%   [y] = flipvec(x)
%
% INPUTS
%   * x         -- vector of values.
%
% OUTPUTS
%   * y         -- vector of values in reversed order.
%
% DESCRIPTION
%   Function checks whether item is a vector (row or column), i.e.
%   has ndim = 2 and a singleton dimension, and then
%   flips the vector (i.e. reverse order).
%
% ERRORS
%    WMTSA:InvalidNumArguments
%    WMTSA:NotAVector
%
% SEE ALSO
%   isvector
%
% TOOLBOX
%   wmtsa/utils
%
% CATEGORY
%
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2004-Apr-27
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

%   $Id: flipvec.m 612 2005-10-28 21:42:24Z ccornish $

usage_str = ['Usage:  [y] = ', mfilename, '(x)'];

%% Check arguments
error(nargerr(mfilename, nargin, [1:1], nargout, [0:1], 1, usage_str, 'struct'));
                      

[is_a_vector, nsdim] = wmtsa_isvector(x);
if (~is_a_vector)
  error('WMTSA:notAVector', ...
        encode_errmsg('WMTSA:notAVector', ...
                      wmtsa_err_table, 'x'));;
end

if (isempty(nsdim))
  % A point vector
  y = x;
else
  % A true vector
  y = flipdim(x, nsdim);
end
  

return
