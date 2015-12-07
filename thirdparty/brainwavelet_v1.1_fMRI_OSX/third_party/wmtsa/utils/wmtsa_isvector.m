function [tf, nsdim] = wmtsa_isvector(x, type)
% wmtsa_isvector -- Determine if item is a vector.
%
%****f* wmtsa.utils/wmtsa_isvector
%
% NAME
%   wmtsa_isvector -- Determine if item is a vector.
%
% SYNOPSIS
%   [tf, nsdim] = wmtsa_isvector(x, [type])
%
% INPUTS
%   * x          -- item to check (object).
%   * type       -- (optional) type of vector (character string).
%
% OUTPUTS
%   * tf         -- flag indicating whether item is a vector (Boolean).
%   * nsdim      -- the non-singleton dimension of the vector (integer).
%
% DESCRIPTION
%   Function checks if the item is a vector by determining whether it has:
%   * two dimensions
%   * at least one singleton dimension (length of dimension = 1).
% 
%   The optional input argument 'type' specifies whether to check for a 
%   particular type of vector.  Valid values for type include:
%   * 'row'            -- row vector with a singleton dimension of 1
%   * 'col','column'   -- column vector with a singleton dimension of 2
%   * 'nonsingleton', 'truevector'-- a vector having one non-singleton dimension, i.e.
%                         either a row or column vector.
%   * 'point'          -- the degenerate case where both dimensions are singletons.
%   There is no default value for 'type'.  If 'type' is not specified, any vector 
%   (row, column, point) returns a Boolean true.
%                         
%   If the output argument 'nsdim' is specified, the ordinal value of the
%   non-singleton dimension of the vector is returned.  
%   Valid values for nsdim are:
%   * 1       -- first dimension is non-singleton, i.e. a row vector
%   * 2       -- second dimension is non-singleton, i.e. a column vector
%   * <empty> -- both dimensions are singleton, i.e. a 'point' vector.
%
% USAGE
%   tf = wmtsa_isvector(x)
%
%   tf = wmtsa_isvector(x, type)
%
%   [tf, nsdim] = wmtsa_isvector(x)
%
% ERRORS
%   WMTSA:InvalidNumArguments, WMTSA:InvalidArgumentValue 
%
% EXAMPLE
%   x = [1:10];
%     % A row vector
%   tf = wmtsa_isvector(x)
%     % Result: tf =  1
%   tf = wmtsa_isvector(x, 'row')
%     % Result: tf =  1
%   tf  = wmtsa_isvector(x, 'col')
%     % Result: tf =  0
%   y = x';
%     % y is a column vector.
%   tf  = wmtsa_isvector(y, 'row')
%     % Result: tf =  0
%   tf  = wmtsa_isvector(y, 'col')
%     % Result: tf =  1
%   [tf, nsdim] = wmtsa_isvector(y)
%     % Result: tf =  1, nsdim = 1
%   [tf, nsdim] = wmtsa_isvector(x)
%     % Result: tf =  1, nsdim = 2
%
% NOTES
%   1.  Starting with version 7, MATLAB features a isvector function.
%       wmtsa_isvector is compatiable with MATLAB version but supplies 
%       additional functionality.
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
%   (c) 2004, 2005 Charles R. Cornish
%
% CREDITS
%
%
% REVISION
%   $Revision: 612 $
%
%***

%   $Id: wmtsa_isvector.m 612 2005-10-28 21:42:24Z ccornish $

usage_str = ['Usage:  [tf, nsdim] = ', mfilename, ...
             '(x, [type])'];

%% Check arguments
error(nargerr(mfilename, nargin, [1:2], nargout, [0:2], 1, usage_str, 'struct'));

tf = 0;
nsd = [];

tf = ((ndims(x) == 2) && ...
      ((size(x,1) == 1) || (size(x,2) == 1)));

if (tf && ( (nargout > 1) || exist('type', 'var') ) )
  if (size(x,1) ~= 1)
    nsd = 1;
  elseif (size(x,2) ~= 1)
    nsd = 2;
  else
    nsd = [];
  end

  if(exist('type', 'var'))
    switch type
   
     case 'row'
        tf = (nsd == 2);
   
     case {'col', 'column'}
        tf = (nsd == 1);
     
     case {'nonsingleton', 'truevector'}
      tf = ~isempty(nsd);
   
     case 'point'
      tf = isempty(nsd);
    
     otherwise
      error(['WMTSA:', mfilename, ':invalidArgumentValue'], ...
            encode_errmsg('WMTSA:invalidArgumentValue', ...
                          wmtsa_err_table, 'type',  num2str(type)));
    end
  end
end


if (nargout > 1)
  nsdim = nsd;
end


return
