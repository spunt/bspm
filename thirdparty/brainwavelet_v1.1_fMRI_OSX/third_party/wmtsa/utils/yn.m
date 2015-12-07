function [tf] = yn(x, return_format)
% yn -- Determine Boolean value.
%
%****f* wmtsa.utils/yn
%
% NAME
%   yn -- Determine Boolean value.
%
% SYNOPSIS
%   [tf] = yn(x, [return_format])
%
% INPUTS
%   * x             -- value to check (Boolean numeric or character string)
%   * return_format -- (optional) format of returned value, see DESCRIPTION
%                      for details (string).
%
% OUTPUTS
%  * tf             -- Boolean result (Boolean numeric or character string).
%
% SIDE EFFECTS
%
%
% DESCRIPTION
%  yn evaluates the value of x and determines whether is a reduces to a
%  Boolean (true/false) value.  It is a handy utlity for converting between
%  formats for Boolean values, i.e. T/F -> true/false -> 1/0 -> Y/N -> yes/no,
%  as well as for evaluating expressions within the callers workspace. 
% 
%  When x is numeric, the function checks the values of the array and
%  returns a Boolean value depending on whether all values are zero (false)
%  or otherwise (true).
%  
%  When x is a character string and has the values:
%  * 'T', 'TRUE',  'YES', 'Y' (case-insensitive), it returns a value of 'true'
%  * 'F', 'FALSE', 'NO', 'N' (case-insensitive),  it returns a value of 'false'
%
%  For other characters strings, the function attempts to evaluate the 
%  expression using the specified variables contained in caller's workspace.  
%  If the expression can  not be evaluated, an error is thrown.
% 
%  If x is neither numeric or character, an error is thrown.
%
%  The 'return_format' input argument controls format of the output argument.
%  Valid values for 'return_format' are case-insenstive and include:
%  * 'numeric'    -- returns 1/0 for Boolean true/false
%  * 'tf'         -- returns T/F for Boolean true/false
%  * 'truefalse'  -- returns TRUE/FALSE for Boolean true/false
%  * 'yn'         -- returns Y/N for Boolean true/false
%  * 'yesno'      -- returns YES/NO for Boolean true/false.
%  The default value of 'return_format' is 'numeric'.
%
% USAGE
%   tf = yn(x)
% 
%   tf = yn(x, return_type)
%
% ERRORS
%   WMTSA:InvalidNumArguments, WMTSA:InvalidStringArgumentValue,
%   WMTSA:InvalidArgumentType, WMTSA:InvalidArgumentValue
%
% EXAMPLE
%   tf = yn(1)
%     % Returns tf = 1
%   tf = yn(1, 'tf')
%     % Returns tf = 'T"
%   tf = yn('F')
%     % Returns tf = 0
%   tf = yn('T', 'yesno')
%     % Returns tf = 'YES'
%   x = 1;
%   y = 0;
%   tf = yn('x == y', 'truefalse')
%     % Returns tf = 'FALSE'
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
%   2005-02-16
%
% COPYRIGHT
%   (c) 2005 Charles R. Cornish
%
% CREDITS
%
%
% REVISION
%   $Revision: 612 $
%
%***

%   $Id: yn.m 612 2005-10-28 21:42:24Z ccornish $

defaults.return_format = 'numeric';  

usage_str = ['Usage:  [tf] = ', mfilename, ...
             '(x, [return_format])'];
  
error(nargerr(mfilename, nargin, [1:2], nargout, [0:1], 1, usage_str, 'struct'));

if (~exist('return_format', 'var') || isempty(return_format))
  return_format = defaults.return_format;
end

if (isnumeric(x))
  if (x)
    tf = 1;
  else
    tf = 0;
  end
elseif (ischar(x))
  switch upper(x)
   case {'T', 'TRUE', 'YES', 'Y'}
    tf = 1;
   case {'F', 'FALSE', 'N', 'NO'}
    tf = 0;
   otherwise
    try
      tf = evalin('caller', x);
    catch
      error('WMTSA:InvalidStringArgumentValue', ...
            ['Cannot evaluate (', x, ') into Boolean.']);
    end
  end
else
  error('WMTSA:InvalidArgumentType', ...
        ['Cannot evaluate x as Boolean.']);
end

switch lower(return_format)
 case 'numeric'
  % Do nothing, return as is
 case 'tf'
  if (tf)
    tf = 'T';
  else
    tf = 'F';
  end
 case 'yn'
  if (tf)
    tf = 'Y';
  else
    tf = 'N';
  end
 case 'truefalse'
  if (tf)
    tf = 'TRUE';
  else
    tf = 'FALSE';
  end
 case 'yesno'
  if (tf)
    tf = 'YES';
  else
    tf = 'NO';
  end
 otherwise
  error('WMTSA:InvalidArgumentValue', ...
        ['Unrecognzied value for ''return_format'' argument', ...
         '(', return_value, ').']);
end

return
  
  
    
    



return
