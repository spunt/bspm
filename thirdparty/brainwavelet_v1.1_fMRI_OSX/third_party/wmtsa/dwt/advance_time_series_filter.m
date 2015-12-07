function nu = advance_time_series_filter(wtfname)
% advance_time_series_filter -- Calculate the advance of the time series or filter for a given wavelet.
%
%****f* wmtsa.dwt/advance_time_series_filter
%
% NAME
%   advance_time_series_filter -- Calculate the advance of the time series or filter for a given wavelet.
%
% SYNOPSIS
%   nu = advance_time_series_filter(wtfname)
%
% INPUTS
%   wtfname      -  string containing name of WMTSA-supported wavelet filter.
%
% OUTPUTS
%   nu           -  advance of time series for specified wavelet filter.
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
%  
%  For Least Asymmetric filters, equation 112e of WMTSA:
%   nu =   -L/2 + 1,   for L/2 is even;
%      =   -L/2,       for L = 10 or 18;
%      =   -L/2 + 2,   for L = 14.  
%
%  For Best Localized filter, page 119 of WMTSA.
%   nu =   -5,         for L = 14;
%      =   -11,        for L = 18;
%      =   -9,         for L = 20.
%
%  For Coiflet filters, page 124 and equation 124 of WMTSA:
%   nu =   -2*L/3 + 1
%
% REFERENCES
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%     Time Series Analysis. Cambridge: Cambridge University Press.
%
% SEE ALSO
%   dwt_filter
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

% $Id: advance_time_series_filter.m 612 2005-10-28 21:42:24Z ccornish $


usage_str = ['Usage:  nu = ', mfilename, ...
             ' (wtfname)'];

%%  Check input arguments and set defaults.
error(nargerr(mfilename, nargin, [1:1], nargout, [0:1], 1, usage_str, 'struct'));

% Check for valid wavelet and get wavelet filter coefficients
try
  wtf_s = dwt_filter(wtfname);
catch
  rethrow(lasterror);
end
h = wtf_s.h;
g = wtf_s.g;
L = wtf_s.L;


nu = NaN;

switch lower(wtfname)
  
 % Haar filter
 case {'haar'}
  nu = 0;
 
 % Haar filter
 % Value from Figure 115
 case {'d4'}
  nu = -1;

 % Extremal Phase filters
% case {'haar', 'd4', 'd6', 'd8', 'd12', 'd14', 'd16', 'd18', 'd20'}
 case { 'd6', 'd8', 'd12', 'd14', 'd16', 'd18', 'd20'}
  error('Need to determine nu for Extremal Phase filters  -  Is it -1 for all filters?.');
  
 % Least Asymmetric filters
 case { 'la8', 'la10', 'la12', 'la14', 'la18', 'la16', 'la20'}
  nu = advance_least_asymetric_filter(L);

 % Best Localized filters
 case {'bl14', 'bl18', 'bl20'}
  nu = advance_best_localized_filter(L);

 % Coiflet filters
 case {'c6', 'c12', 'c18', 'c24', 'c30'}
  nu = advance_coiflet_filter(L);

 otherwise
  % Do nothing 
end

return

function nu = advance_least_asymetric_filter(L)
% Equation 112c of WMTSA.
  if (mod(L/2, 2) == 0)
    % L/2 is even, i.e. L = 8, 12, 16, 20
    nu = -(L/2) + 1;
  else
    switch (L)
     case {10, 18}
      nu = -(L/2);
     case {14}
      nu = -(L/2) + 2;
     otherwise
      error(['Invalid filter length (L = ', int2str(L), ') ' ...
                        'specified.']);
    end
  end

return

function nu = advance_best_localized_filter(L)
% Page 119 of WMTSA.
  switch L
   case 14
    nu = -5;
   case 18
    nu = -11;
   case 18
    nu = -9;
   otherwise
    % Do nothing
  end
return
  
function nu = advance_coiflet_filter(L)
%  Page 124 and equation 124 of WMTSA:
  nu = -(2 * L / 3) + 1;
return
  


  
