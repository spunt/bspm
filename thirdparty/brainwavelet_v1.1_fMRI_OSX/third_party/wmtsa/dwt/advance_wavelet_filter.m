function nuHj = advance_wavelet_filter(wtfname, j)
% advance_wavelet_filter -- Calculate the advance of the wavelet filter at jth level for a given wavelet.
%
%****f* wmtsa.dwt/advance_wavelet_filter
%
% NAME
%   advance_wavelet_filter -- Calculate the advance of the wavelet filter at jth level for a given wavelet.
%
% SYNOPSIS
%   nuHj = advance_wavelet_filter(wtfname, j)
%
% INPUTS
%   wtfname      = string containing name of WMTSA-supported wavelet filter.
%   j            = jth level (index) of scale or a range of j levels of scales
%                  (integer or Jx1 vector of integers).
%
% OUTPUTS
%   nuHj         = advance of wavelet filter at jth level
%                  (integer or vector of integers).
%
% SIDE EFFECTS
%   wavelet is a WMTSA-supported wavelet filter; otherwise error.
%
% DESCRIPTION
%
%
% EXAMPLE
%
%
% ALGORITHM
%   nuHj = - (2^(j-1) * (L-1) + nu);
%
%   For details, see equation 114b of WMTSA.
%
% REFERENCES
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%     Time Series Analysis. Cambridge: Cambridge University Press.
%
% SEE ALSO
%   advance_time_series_filter, dwt_filter
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2003-05-08
%
% COPYRIGHT
%
%
% REVISION
%   $Revision: 612 $
%
%***

% $Id: advance_wavelet_filter.m 612 2005-10-28 21:42:24Z ccornish $

usage_str = ['Usage:  nuHj = ', mfilename, ...
             ' (wtfname, j)'];

%%  Check input arguments and set defaults.
error(nargerr(mfilename, nargin, [2:2], nargout, [0:1], 1, usage_str, 'struct'));

% Check for valid wavelet and get wavelet filter coefficients
try
  wtf_s = dwt_filter(wtfname);
catch
  rethrow(lasterror);
end
h = wtf_s.h;
g = wtf_s.g;
L = wtf_s.L;
  
error(argterr(mfilename, j, 'int0', [], 1, '', 'struct'));


nuHj = [];

nu = advance_time_series_filter(wtfname);

nuHj = -(2.^(j-1) * (L-1) + nu);

% Return as column vector
nuHj = nuHj(:);

return


