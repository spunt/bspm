function [err] = nargerr(func, ...
                         numargin, exp_numargin, ...
                         numargout, exp_numargout, ...
                         mode, msg, return_type)
% nargerr -- Check number of arguments to a function.
%
%****f* wmtsa.utils/nargerr
%
% NAME
%    nargerr -- Check number of arguments to a function.
%
% SYNOPSIS
%   [errmsg] = nargerr(func, numargin, argin, numargout, argout, [mode], [msg])
%   [errmsg] = nargerr(func, numargin, argin, numargout, argout, [mode], [msg], 'string')
%   [errstruct] = nargerr(func, numargin, argin, numargout, argout, [mode], [msg], 'struct')
%
% INPUTS
%   * func          -- name of checked function (string or function handle).
%   * numargin      -- number of input arguments passed during function call (integer).
%   * exp_numargin  -- expected number of input arguments.
%                      (integer or range of integers).
%   * numargout     -- number of output arguments passed during function call (integer).
%   * exp_numargout -- expected number of output arguments.
%                      (integer or range of integers).
%   * mode          -- (optional) output display mode (integer).
%                        0 = silent
%                        1 = verbose (default).
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
%   Function call requires a minimum of 5 input arguments; otherwise error.
%
% DESCRIPTION
%   nargerr checks whether the number of input and output arguments of a function
%   call each are within the valid range. If the number of input (numargin) or 
%   output (numargout) arguments are not within range of expected values for input 
%   (exp_numargin) or output (exp_numargout) arguments, respectively, the 
%   function returns an error message string or error message struct.
%
%   The 'exp_numarg' arguments may be specified as a(n):
%   * integer
%   * range of integers
%   * string expression for range of integers
%   For range intergers, e.g. [1:3], the value of 0 may be used
%   for the lower bound and Inf for the upper bound to specify
%   no lower or upper limit, respectively.
%   The string expression option also allows specification of no 
%   lower or upper bounds, e.g.
%     '2:'  -- A minimum of two arguments and no upper bound.
%     ':2'  -- A maximum of two arguments.
%
%   'mode' controls whether information is displayed to the command window.
%
%   'msg' is an message to be displayed to the command window, if nonempty
%   and regardless of the value of 'mode'.  This allows a user-defined customized
%   message for display.
%   
%   'return_type' specifies the format of the output argument:
%      'string' -- output is an error message (errmsg).
%      'struct' -- output is a struct (errstruct) with fields message and identifier.
%    The 'struct' option may be used in with the error function to throw an
%    error.  In this example,
%
%       error(nargerr('myfunction', 2, 1, 2, 0, '', '', 'struct'))
%
%    the number of input and output arguments are outside the expected range.
%    The nargerr function returns an error struct, which in turn causes the 
%    error function to throw an error.
%
% ERRORS
%   WMTSA:nargerr:invalidInputArgumentValue
%   WMTSA:nargerr:incorrectStringFormat
%   MATLAB:nargchk:notEnoughInputs     (thrown by nargchk)
%   MATLAB:nargchk:tooManyInputs       (thrown by nargchk)
%   MATLAB:nargoutchk:notEnoughOutputs (thrown by nargchk)
%   MATLAB:nargoutchk:tooManyOutputs   (thrown by nargchk)
%
% EXAMPLE
%   % Default execution.
%   [errmsg] = nargerr(mfilename, nargin, [1:3], nargout, [0:2]);
%
%   % No display to command window.
%   [errmsg] = nargerr(mfilename, nargin, [1:3], nargout, [0:2], 0);
%
%   % Display optional msg to command window.
%   msg = 'Usage: myfuction a b c [d]';
%   [errmsg] = nargerr(mfilename, nargin, [1:3], nargout, [0:2], 0, msg);
%
%   % Use strings to specify range of integers with ':' syntax
%   [err, errmsg] = nargerr(mfilename, nargin, ':3', nargout, '');
%
%   % Throw an error if one found.
%   error(nargerr(mfilename, 2, 1, 2, 0, '', '', 'struct'))
%
% NOTES
%   1. This function requires MATLAB 7, which allows error structures as input
%      to the error function.
%   2. To skip number of argument checking for input or output arguments, specify
%      exp_numargin or exp_numargout, respectively, as empty vectors ([]) 
%      or empty strings ('').
%
% SEE ALSO
%   nargchk, nargoutchk, errargn
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
%   nargerr was inspired by the errargn function
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

% $Id: nargerr.m 612 2005-10-28 21:42:24Z ccornish $

%% Set Defaults
err = [];
valid_mode = [0 1];
valid_return_type = {'string', 'struct'};

usage_str = ['Usage: [err] = ', mfilename, ...
             '(func, numargin, exp_numargin, numargout, exp_numargout, ' ...
             '[mode], [msg], [return_type])'];


%% Check arguments.
err = nargchk(5, 8, nargin, 'struct');
if (~isempty(err))
  error(err.identifier, err.message);
end

  

if (ischar(func))
  funcname = func;
elseif (isa(func, 'function_handle'))
  funcname = func2str(func);
else
  error('WMTSA:nargerr:invalidInputArgumentValue', ...
        'func must be a string or function handle');
end

  

if (~exist('mode', 'var') || isempty(mode))
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

%% Check number of input arguments.
if (isempty(exp_numargin))
  % Do nothing -- skip
else
  err = check_numarg(numargin, exp_numargin, 'input', return_type);
end


%% Check number of output arguments.
if (isempty(err) && ~isempty(exp_numargout))
  err = check_numarg(numargout, exp_numargout, 'output', return_type);
end

%% Display optional fancy error message and additional message.
if (~isempty(err))
  if (mode)
    if (isstruct(err))
      errmsg = err.message;
    else
      errmsg = err;
    end
    display_errmsg(funcname, errmsg);
  end
  if(exist('msg', 'var'))
    disp(msg);
  end
end

return

function err = check_numarg(numarg, exp_numarg, argtype, return_type)
%% Function to check the arguments to func.
  if (isnumeric(exp_numarg))
    low = min(exp_numarg);
    high = max(exp_numarg);
  elseif (ischar(exp_numarg) & ~isletter(exp_numarg))
    delpos = strfind(exp_numarg, ':');
    if (isempty(delpos))
      low = str2num(exp_numarg);
      high = str2num(exp_numarg);
    else
      if (~isempty(delpos))
        low = str2num(exp_numarg(1:delpos-1));
        if (isempty(low))
          low = 0;
        end
        high = str2num(exp_numarg(delpos+1:end));
        if (isempty(high))
          high = Inf;
        end
      end
    end
    else
    err_s.message = 'Cannot parse expected number argument string';
    err_s.identifier = ['WMTSA:nargerr:incorrectStringFormat'];
    error(err_s);
  end

  switch argtype
   case 'input'
    hnargchk = @nargchk;
   case 'output'
    hnargchk = @nargoutchk;
  end
  err = hnargchk(low, high, numarg, return_type);
return


function display_errmsg(funcname, errmsg)
%% Function to display 'fancy' error message.
  
  blank = ' ';

  if (~isempty(funcname))
    msg = [ ' ' funcname ':  ' xlate(errmsg)];
  else
    msg = [ ' ' xlate(errmsg)];
  end
  msg_length = length(msg);

  banner_length = msg_length + 2;
  star_str = '*';
  star_str = star_str(:, ones(1, banner_length));
  minus_str = '-';
  minus_str = minus_str(:, ones(1, banner_length));
  
  disp(blank);
  disp(star_str);
  disp('ERROR ...');
  disp(minus_str);
  disp(msg);
  disp(star_str);
  disp(blank);

return
