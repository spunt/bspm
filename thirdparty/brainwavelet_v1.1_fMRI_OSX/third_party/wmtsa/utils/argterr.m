function [err] = argterr(func, arg, datatypes, arg_size, mode, msg, return_type)
% argterr - Check the data type and size of a function argument.
%
%****f* wmtsa.utils/argterr
%
% NAME
%   argterr - Check the data type and size of a function argument.
%
% SYNOPSIS
%   [errmsg] = argterr(func, arg, type, [arg_size], [mode])
%   [errmsg] = argterr(func, arg, type, [arg_size], [mode], [msg], 'string')
%   [errstruct] = argterr(func, arg, type, [arg_size], [mode], [msg], 'struct')
%
% INPUTS
%   * func          -- checked function (string or function handle).
%   * arg           -- the argument to check (object)
%   * datatypes     -- expected data type(s) of the arg.
%                      See verify_datatypes function for possibles datatypes to check.
%                      (string or string cell array).
%   * arg_size      -- (optional) expected size of array  (integer vector).
%   * mode          -- (optional) output display mode (charcter string)
%   * msg           -- (optional) message string to be displayed (string).
%   * return_type   -- (optional) type of output argument (string).
%                       'string' = return error message string (default).
%                       'struct' = return error message struct.
%
% OUTPUTS
%   * errmsg        -- error message (string).
%   * errstruct     -- error struct with fields:  message, identifier.
%
% SIDE EFFECTS
%   Function call requires a minimum of three input arguments; otherwise error.
%
% DESCRIPTION
%   argterr checks the data type(s) and optionally the size of an argument to a 
%   function call.  If arg does not have the specified data type(s), an error
%   message or error struct is returned.
%
%   arg may be a single
%
%   Possible values for 'datatypes' to check include:
%   * 'posint'              -- All are positive integers --> integer value(s) > 0.
%   * 'int0'                -- All are positive integers plus zero --> integer value(s) >= 0.
%   * 'int','integer'       -- All are integers --> any integer value(s).
%   * 'num','numeric'       -- All are numeric.
%   * 'struct','structure'  -- Is a structure.
%   * 'char','character','string' - Is a character string.
%   * 'scalar               -- Is a point (size of all dimensions = 1).
%   * 'vec','vector'        -- Is a vector (i.e. MxN, with M and/or N = 1).
%   * 'nonsingleton','truevector'  -- Is a vector (i.e. MxN with M *or* N = 1).
%   * 'row','rowvector'     -- Is a row vector (i.e. M x 1).
%   * 'col','columnvector'  -- Is a column vector (i.e. 1 x N).
%   * 'finite'              -- All are finite.
%   * 'nonsparse'           -- Is a non-sparse matrix.
% 
%   The input argument 'mode' controls the display of diagnostic information
%   to the consolue and has possible values:
%   * 0 -- silent
%   * 1 -- verbose
%   The default value is 1 (silent).
% 
% EXAMPLE
%   arg = [0 2];
%     % arg is an integer.
%   err = argterr('myfunction', arg, 'posint')
%     % Result: err = 1
%   err = argterr('myfunction', arg, 'int', [1 2], '')
%     % Result: err = 0
%   [err, errmsg] = argterr('myfunction', arg, {'int', 'vector'})
%     % Result: err = 0
%
% NOTES
%   argterr is modelled after errargt, which is part of MATLAB wavelet toolkit.
%
% SEE ALSO
%   verify_datatype
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2003-06-22
%
% COPYRIGHT
%   (c) 2003, 2004, 2005 Charles R. Cornish
%
% CREDITS
%   argterr is inspired by the function errargt 
%   by M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi
%   which is part of the MATLAB wavelet toolkit.
%
% MATLAB VERSION
%   7.0
%
% REVISION
%   $Revision: 612 $
%
%***

% $Id: argterr.m 612 2005-10-28 21:42:24Z ccornish $
  
%% Set Defaults
err = [];
valid_mode = [0 1];
valid_return_type = {'string', 'struct'};

usage_str = ['Usage: [err, errmsg] = ', mfilename, ...
               '(func, arg, datatypes, [arg_size], [mode], [msg], [return_type])'];

%% Check arguments
error(nargerr(mfilename, nargin, [3:7], nargout, [0:1], ...
              1, usage_str, 'struct'));

if (ischar(func))
  funcname = func;
elseif (isa(func, 'function_handle'))
  funcname = func2str(func);
else
  error('WMTSA:nargerr:invalidInputArgumentValue', ...
        'func must be a string or function handle');
end

arg_name = inputname(2);

if (iscellstr(datatypes))
  % OK - do nothing
elseif (ischar(datatypes))
  % Convert to character cell array.
  datatypes = {datatypes};
else
  error('WMTSA:invalidArgumentType', ...
        ['Argument ''datatypes'' must be a charcter string or charcter cell ' ...
         'array.']);
end

if (~exist('mode', 'var'))
  mode = 1;
end

if (isempty(find(mode == valid_mode)))
  err_s.message = ['Invalid value for input argument ', 'mode'];
  err_s.identifier = ['WMTSA:nargerr:invalidInputArgumentValue'];
  error(err_s);
end

if (~exist('return_type', 'var') || isempty(return_type))
  return_type = 'string';
end

if (isempty(strmatch(return_type, valid_return_type)))
  err_s.message = ['Invalid value for input argument ', 'return_type'];
  err_s.identifier = ['WMTSA:nargerr:invalidInputArgumentValue'];
  error(err_s);
end

% Check argument datatypes
[tf, errmsg] = verify_datatype(arg, datatypes, arg_name);

err_s = [];

if (~tf)
  err_s.identifier = ['WMTSA:argterr:invalidArgumentDataType'];
  err_s.message = errmsg;
else
  % Check argument dimension sizes
  if (exist('arg_size', 'var') && ~isempty(arg_size))
    actual_arg_size = size(arg);
    if (~isequal(actual_arg_size, arg_size))
      err_s.identifier = ['WMTSA:argterr:incorrectArgumentDimensions'];
      err_s.message = ['Expected and actual sizes of ', arg_name, ' are not equal'];
    end
  end
end

if (~isempty(err_s))

  if (mode)
    display_errmsg(funcname, err_s.message)
  end

  if(exist('msg', 'var'))
    disp(msg);
  end

  if (strmatch(return_type, 'string'))
    err = err_s.message;
  elseif (strmatch(return_type, 'struct'))
    err = err_s;
  end
end


return


function display_errmsg(funcname, errmsg)
  
  blank = ' ';

  if (~isempty(deblank(funcname)))
    msg = ['  ', funcname, ': ', xlate(errmsg)];
  else
    msg = errmsg;
  end
  
  msg_length = length(msg);
  banner_length = msg_length + 2;
  star_str = '*';
  star_str = star_str(:, ones(1, banner_length));
  minus_str = '-';
  minus_str = minus_str(:, ones(1, banner_length));
  
  
  disp(blank);
  disp(star_str);
  disp('ARGUMENT TYPE ERROR ...');
  disp(minus_str);
  disp(msg);
  disp(star_str);
  disp(blank);

return
  
  

