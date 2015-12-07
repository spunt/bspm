function J0 = modwt_choose_nlevels(choice, wtfname, N)
% modwt_choose_nlevels -- Select J0 based on choice, wavelet filter and data series length.
%
%****f* wmtsa.dwt/modwt_choose_nlevels
%
% NAME
%   modwt_choose_nlevels -- Select J0 based on choice, wavelet filter and data series length.
%
% USAGE
%   J0 = modwt_choose_nlevels(choice, wtfname, N)
%
% INPUTS
%   * choice      -- choice for method for calculating J0 (string)
%                    Valid Values:
%                     'conservative'
%                     'max', 'maximum'
%                     'supermax', 'supermaximum'
%   * wtfname     -- wavelet transform filter name (string)
%                    Valid Values:  see modwt_filter
%   * N           -- number of observations.
%
% OUTPUT
%   * J0          -- number of levels (J0) based selection criteria.
%
% SIDE EFFECTS
%   1.  wtfname is a WMTSA-supported MODWT wtfname, otherwise error.
%   2.  N > 0, otherwise error.
%
% DESCRIPTION
%
%
% EXAMPLE
%   J0 = modwt_choose_nlevels('convservative', 'la8', N)
%
% ERRORS  
%   WMTSA:MODWT:InvalidNumLevels    =  Invalid type/value specified for nlevels.
%
% ALGORITHM
%   for 'conservative':              J0  < log2( N / (L-1) + 1)
%   for 'max', 'maximum':            J0 =< log2(N)
%   for 'supermax', 'supermaximum':  J0 =< log2(1.5 * N)
%
%   For further details, see page 200 of WMTSA.
%
% REFERENCES
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%   Time Series Analysis. Cambridge: Cambridge University Press.
%
% SEE ALSO
%   modwt_filter
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2003-05-24   
%
% COPYRIGHT
%   (c) 2003, 2004, 2005 Charles R. Cornish
%
% REVISION
%   $Revision: 612 $
%
%***


% $Id: modwt_choose_nlevels.m 612 2005-10-28 21:42:24Z ccornish $
  
available_choices = {'conservative', 'max', 'supermax'};

usage_str = ['Usage:  [J0] = ', mfilename,'(choice, wtfname, N)'];

error(nargerr(mfilename, nargin, [3:3], nargout, [0:1], 1, usage_str, 'struct'));

% Check for valid wtfname and get wavelet filter coefficients
try
  [wtf] = modwt_filter(wtfname);
catch
  rethrow(lasterror);
end

L = wtf.L;

error(argterr(mfilename, N, 'int0', [], 1, '', 'struct'));



switch choice
 case 'conservative'
  J0 = floor(log2( (N / (L - 1)) - 1));
 case {'max', 'maximum'}
  J0 = floor(log2(N));
 case {'supermax', 'supermaximum'}
  J0 = floor(log2(1.5 * N));
 otherwise
  display_choices;
  error('WMTSA:invalidNLevelsValue', ...
        encode_errmsg('WMTSA:invalidNLevelsValue', wmtsa_err_table))
end

return


function display_choices

available_choices = {'conservative', 'max', 'supermax'};

disp('Choices available for selecting nlevel:');
for ( i = 1:length(available_choices))
  disp(['  ', available_choices{i}]);
end

return;
