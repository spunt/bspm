function err_table = wmtsa_err_table

% $Id: wmtsa_err_table.m 612 2005-10-28 21:42:24Z ccornish $
% Error table for WMTSA toolkit.

% Example
%  encode_errmsg('WMTSA:invalidArgumentDataType', wmtsa_err_table, 'wtf'));
%   
  
% Define the error table fields
err_table = struct( 'err_id',    {}, ...
                    'err_msg',   {}, ...
                    'err_descr', {}, ...
                    'nargs',     {});

% Template for a new entry
% err_table(n).err_id    = '';
% err_table(n).err_msg   = '';
% err_table(n).err_descr = '';
% err_table(n).nargs = 0;

err_table(1).err_id    = 'WMTSA:unknownError';
err_table(1).err_msg   = 'An unknown error has occurred.';
err_table(1).err_descr = 'An unknown error has occurred.';
err_table(1).nargs = 0;

err_table(2).err_id    = 'WMTSA:invalidArguments';
err_table(2).err_msg   = '%s';
err_table(2).err_descr = ...
    'Invalid arguments specified in function call.';
err_table(2).nargs = 1;


err_table(3).err_id    = 'WMTSA:invalidNumArguments';
err_table(3).err_msg   = '%s';
err_table(3).err_descr = ...
    'Invalid number of arguments specified in function call.';
err_table(3).nargs = 1;

err_table(4).err_id    = 'WMTSA:invalidArgumentDataType';
err_table(4).err_msg   = 'Argument (%s) has an invalid datatype.';
err_table(4).err_descr = 'Argument has incorrect datatype.';
err_table(4).nargs = 1;
% Example
%  encode_errmsg('WMTSA:invalidArgumentDataType', wmtsa_err_table, 'wtf'));

err_table(5).err_id    = 'WMTSA:invalidArgumentValue';
err_table(5).err_msg   = 'Argument (%s) has an invalid value (%s).';
err_table(5).err_descr = 'Argument has incorrect value.';
err_table(5).nargs = 2;
% Example
% encode_errmsg('WMTSA:invalidArgumentValue', wmtsa_err_table, 'transform', ...
%                        transform));


err_table(6).err_id    = 'WMTSA:invalidArgumentNDims';
err_table(6).err_msg   =  '';
err_table(6).err_descr = 'Argument has incorrect number of dimensions.';
err_table(6).nargs = 0;

err_table(7).err_id    = 'WMTSA:invalidArgumentSize';
err_table(7).err_msg   = '%s %s';
err_table(7).err_descr = 'Argument has incorrect size.';
err_table(7).nargs = 2;

err_table(8).err_id    = 'WMTSA:notAVector';
err_table(8).err_msg   = 'Item %s is not a vector';
err_table(8).err_descr = 'Item is not a vector.';
err_table(8).nargs = 1;

err_table(9).err_id    = 'WMTSA:missingRequiredArgument';
err_table(9).err_msg   = 'Required argument (%s) is not specified.';
err_table(9).err_descr = 'Required argument is not specified.';
err_table(9).nargs = 1;

err_table(10).err_id    = 'WMTSA:invalidNLevelsValue';
err_table(10).err_msg   = 'nlevels must be an positive integer or string.';
err_table(10).err_descr = 'nlevels must be an positive integer or string.';
err_table(10).nargs = 0;

err_table(11).err_id    = 'WMTSA:invalidWaveletTransformFilter';
err_table(11).err_msg   = 'The specified wavelet transofrm (%s) is unknown or unsupported.';
err_table(11).err_descr = 'Invalid wavelet transform filter.';
err_table(11).nargs = 1;

% Template for a new entry
% err_table(n).err_id    = '';
% err_table(n).err_msg   = '';
% err_table(n).err_descr = '';
% err_table(n).nargs = 0;

return
