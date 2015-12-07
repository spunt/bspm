function [tf, msg] = verify_datatype(var, datatypes, var_name)
% verify_datatype -- Verify the datatype(s) of a variable.
%
%****f* wmtsa.utils/verify_datatype
%
% NAME
%   verify_datatype -- Verify the datatype(s) of a variable.
%
% SYNOPSIS
%   [tf, msg] = verify_datatype(var, datatypes, var_name)
%
% INPUTS
%   * var        -- variable to verify (object)
%   * datatypes  -- expected data type(s) of the arg.
%                   See verify_datatypes function for possibles datatypes 
%                   to check (string or cell array of strings).
%   * var_name   -- (optional) alternative name of variable to use in 
%                   error message (string).
%
% OUTPUTS
%   * tf         -- flag indicating whether object as specified datatype(s).
%                   (Boolean).
%   * msg     -- diagnostic error message (string).
%
% SIDE EFFECTS
%   Function call requires a mininum of two input arguments; otherwise error.
%
% DESCRIPTION
%   verify_datatypes checks whether the specified object has the specified
%   data type(s). If the datatypes argument is a cell array of strings, then the
%   function checks whether the variable datatype matches each and all of the 
%   datatypes.  The function returns a logical true if the object's datatypes
%   match those specified by the 'datatypes' argument; otherwise a logical false
%   is returned.
%
%   Possible datatypes to check include:
%   * 'posint'              -- All are positive integers --> integer value(s) > 0.
%   * 'int0'                -- All are positive integers plus zero --> integer value(s) >= 0.
%   * 'int','integer'       -- All are integers --> any integer value(s).
%   * 'real'                -- All are real numbers.
%   * 'num','numeric'       -- All are numeric.
%   * 'struct','structure'  -- Is a structure.
%   * 'char','character','string' - Is a character string.
%   * 'scalar               -- Is a point (size of all dimensions = 1).
%   * 'vec','vector'        -- Is a vector (i.e. MxN, with M and/or N = 1).
%   * 'matrix'              -- Is a matrix (i.e. dim = 2).
%   * 'nonsingleton', 'truevector' -- Is a vector (i.e. MxN with M *or* N = 1).
%   * 'row','rowvector'     -- Is a row vector (i.e. M x 1).
%   * 'col','columnvector'  -- Is a column vector (i.e. 1 x N).
%   * 'finite'              -- All are finite.
%   * 'nonsparse'           -- Is a non-sparse matrix.
%
% EXAMPLE
%
%
% ERRORS
%  WMTSA:InvalidNumArguments, WMTSA:InvalidArgumentType, WMTSA:InvalidArgumentValue
%
% NOTES
%
%
% SEE ALSO
%   argterr
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
%    2005-02-17
%
% COPYRIGHT
%   (c) 2005 Charles R. Cornish
%
% MATLAB VERSION
%   7.0
%
% CREDITS
%
%
% REVISION
%   $Revision: 612 $
%
%***

%   $Id: verify_datatype.m 612 2005-10-28 21:42:24Z ccornish $

usage_str = ['Usage:  [tf, msg] = ', mfilename, ...
             '(var, datatypes, [var_name])'];

%% Check arguments
error(nargerr(mfilename, nargin, [2:3], nargout, [0:2], 1, usage_str, 'struct'));


if (iscellstr(datatypes))
  % OK - do nothing
elseif (ischar(datatypes))
  % Convert to character cell array.
  datatypes = {datatypes};
else
  error('WMTSA:verify_datatype:invalidArgumentType', ...
        ['Argument ''datatypes'' must be a charcter string or cell ' ...
         'array of strings.']);
end

if (~exist('var_name', 'var') || isempty(var_name))
  var_name = inputname(1);
end

%% Initialize output arguments.
tf = 1;
msg = [];

%% Check datatypes against specified list.
for (ii = 1:length(datatypes))  
  type = lower(datatypes{ii});
  switch type
   
   %%%
   %%% Integers
   %%%
   
   case {'int', 'integer'}
    % Any integer
    if (~isnumeric(var) | any(var ~= fix(var)))
      msg = ['Expected ''', var_name, ''' to be an integer.'];
      break;
    end
   
   case 'int0'
    % Positive integers plus zero (integers >= 0)
    if (~isnumeric(var) | any(var < 0) | any(var ~= fix(var)))
      msg = ['Expected ''', var_name, ''' to be an integer with value(s) >= 0.'];
      break;
    end
    
   case 'posint'   
    % Positive integers (integers > 0)
    if (~isnumeric(var) | any(var < 1) | any(var ~= fix(var)))
      msg = ['Expected ''', var_name, ''' to be an integer with value(s) > 0.'];
      break;
    end
    
   %%%
   %%% Numbers
   %%%
  
    case {'real'}
    % Any real number
    if (~isreal(var))
      msg = ['Expected ''', var_name, ''' to be numeric.'];
      break;
    end
   
   case {'num', 'numeric'}
    % Any number
    if (~isnumeric(var))
      msg = ['Expected ''', var_name, ''' to be numeric.'];
      break;
    end
  
   %%%
   %%% Characters/Strings
   %%%
   
   case {'char', 'character', 'string'}
    % A string or character array
    if (~ischar(var))
      msg = ['Expected ''', var_name, ''' to be a string or character array.'];
      break;
    end
  
   case {'scalar'}
    % A scalar
    if (~wmtsa_isscalar(var))
      msg = ['Expected ''', var_name, ''' to be a scalar.'];
      break;
    end

   %%%
   %%% Vectors
   %%%
   
   case {'col', 'columnvector'}
    % A col vector of any datatype
    if (~wmtsa_isvector(var, 'col'))
      msg = ['Expected ''', var_name, ''' to be a column vector.'];
      break;
    end
   
   case {'row', 'rowvector'}
    % A row vector of any datatype
    if (~wmtsa_isvector(var, 'row'))
      msg = ['Expected ''', var_name, ''' to be a row vector.'];
      break;
    end
    
   case {'truevector'}
    % A true vector (i.e. not a point)
    if (~wmtsa_isvector(var, 'vector'))
      msg = ['Expected ''', var_name, ''' to be a true vector.'];
      break;
    end
   
   case {'vec', 'vector'}
    % A vector of any data type
    if (~wmtsa_isvector(var))
      msg = ['Expected ''', var_name, ''' to be a vector.'];
      break;
    end
  
   %%%
   %%% Matrices/Arrays
   %%%

   case {'mat', 'matrix'}
    % A matrix of any data type
    if (~wmtsa_ismatrix(var))
      msg = ['Expected ''', var_name, ''' to be a matrix.'];
      break;
    end
    
   case {'nonsparse'}
    % A nonsparse matrix
    if (~issparse(var))
      msg = ['Expected ''', var_name, ''' to be non-sparse.'];
      break;
    end
   
   case {'finite'}
    % All elements are finite
    if (~isfinite(var))
      msg = ['Expected ''', var_name, ''' to be finite.'];
      break;
    end

   %%%
   %%% Structs
   %%%
   
   case {'struct', 'structure'}
    % A structure
    if (~isstruct(var))
      msg = ['Expected ''', var_name, ''' to be a structure.'];
      break;
    end
    
    
   otherwise
    error('WMTSA:verify_data_type:unknownDataTypeValue', ...
          ['Unrecognzied value for datatype (', type, ').']);
    
  end

end

if (~isempty(msg))
  tf = 0;
end

  
return

