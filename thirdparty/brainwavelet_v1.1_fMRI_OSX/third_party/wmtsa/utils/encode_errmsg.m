function [errmsg] = encode_errmsg(err_id, err_table_ps, varargin)
% encode_errmsg  -- Encode error message for specified err_id.
%
%****f* wmtsa.utils/encode_errmsg
%
% NAME
%   encode_errmsg  -- Encode error message for specified err_id.
%
% USAGE
%   [errmsg] = encode_errmsg(err_id, err_table_ps, [varargin])
%
% INPUTS
%   * err_id         -- error message id (character string).
%   * err_table_ps   -- error lookup table (character string or struct).
%   * varargin       -- (optional) supplemental values to encode in errmsg.
%
% OUTPUTS
%   * errmsg         -- encoded error message (character string).
%
% SIDE EFFECTS
%
%
% DESCRIPTION
%   encode_errmsg encodes an error message (errmsg) for a specified error message 
%   id (err_id).  The function loads the error message table (err_table) from the
%   specified path, searches for the matching (err_table.err_id) entry and 
%   returns the error message template (err_table.errmsg).  Based on the number 
%   of message arguments (err_table.nargs), the function encodes the errmsg using
%   the variable number of arguments (varargin) passed on the function call.
%
%   The err_table is a structure array with the following fields:
%   * err_id     -- error message id (character string).
%   * err_msg    -- error message template (character string).
%   * nargs      -- number of supplement arguments to use for encoding errmsg.
%
%   The input argument 'err_table_ps' may be either a character string or a struct.
%   If a  character string, err_table_ps is full path to a function or script 
%   containing the error table struct to run and load.  If a struct, then the 
%   struct passed as the value in err_table_ps argument is used.
%
% EXAMPLE
%   % Error table load via function wmtsa_err_table.
%   % Name of required argument is 'transform'.
%   errmsg = encode_errmsg('WMTSA:missingRequiredArgument', ...
%                           wmtsa_err_table, 'transform'));
%
% WARNINGS
%
%
% ERRORS
%
%
% NOTES
%
%
% SEE ALSO
%
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
%   2004-06-25
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

%   $Id: encode_errmsg.m 612 2005-10-28 21:42:24Z ccornish $


%% Check arguments
error(nargerr(mfilename, nargin, '2:', nargout, [0:1], 1, '', 'struct'));

%% Load error lookup table.
if (ischar(err_table_ps))
  run(err_table_ps);
elseif (isstruct(err_table_ps))
  err_table = err_table_ps;
end


%% Get list of err_id's.
err_id_list = {err_table.err_id};

%% And match requested err_id.
err_num = strmatch(err_id, err_id_list, 'exact');
if (isempty(err_num))
  error('WMTSA:encode_errmsg:unknownErrId', ...
        ['No match for err_id (', err_id, ').']);
end

%% Write the error message to string.
nargs = err_table(err_num).nargs;
args = cell([1 nargs]);
args(1:length(varargin)) = varargin;

err_msg_template = err_table(err_num).err_msg;

[errmsg, spf_errmsg] = sprintf(err_msg_template, args{:});

if (~isempty(spf_errmsg))
  error('WMTSA:encode_errmsg:unknownError', ... 
        ['Error occurred while encoding errmsg string:  ', spf_errmsg]);
end  

return

