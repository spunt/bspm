function [wtf] = dwt_filter(wtfname)
% dwt_filter -- Define the DWT filter coefficients.
%
%****f* wmtsa.dwt/dwt_filter
%
% NAME
%   dwt_filter -- Define the DWT filter coefficients.
%
% SYNOPSIS
%  [wtf] = dwt_filter(wtfname)
%
% INPUTS
%   * wtfname    -- name of wavelet transform filter (string, case-insenstive).
%
% OUTPUTS
%   * wtf        -- wavelet tranform filter struct (wtf_s).
%
% SIDE EFFECTS
%   wtfname is a valid wavelet filter name; otherwise error.
%
% DESCRIPTION
%   dwt_filter returns a wtf_s struct containing the DWT wavelet (high-pass)
%   and scaling (low-pass) filter coefficients.
%
% NOTES
%   dwt_filter is deprecated by the wtfilter function.  dwt_filter is a
%   pass thru to wtfilter and maintained for backward compatiablity and convenice.
%
%   The wtf_s struct has fields:
%   * g         -- scaling (low-pass) filter coefficients (vector).
%   * h         -- wavelet (high-pass) filter coefficients (vector).
%   * L         -- filter length (= number of coefficients) (integer).
%   * name      -- name of wavelet filter (string).
%   * wtfclass  -- class of wavelet filters (string).
%   * transform -- name of transform (string).
%
%   Typing dwt_filter('list') displays a list of supported filters.
% 
%   Typing dwt_filter('all') returns a struct array of wtf_s of all 
%   supported filters.
%
% ERRORS  
%   WMTSA:InvalidNumArguments
%
% EXAMPLE
%    [h, g, L, name] = dwt_filter('la8');
%
% NOTES
%   dwt_filter is a wrapper function around the wtfilter function.  
%
% REFERENCES
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%     Time Series Analysis. Cambridge: Cambridge University Press.
%
% SEE ALSO
%   wtfilter, modwt_filter
%
% TOOLBOX
%   wmtsa/dwt
%
% CATEGORY
%   Filters: Filters
%

% AUTHOR
%   Charlie Cornish
%   Brandon Whitcher
%
% CREATION DATE
%   2003-09-18
%
% COPYRIGHT
%   (c) Charles R. Cornish 2005
%
% CREDITS
%   Based on the original function (myfilter.m) by Brandon Whitcher.
%
% REVISION
%   $Revision: 612 $
%
%***

% $Id: dwt_filter.m 612 2005-10-28 21:42:24Z ccornish $

usage_str = ['Usage:  [wtf] = ', mfilename, ...
             ' (wtfname)'];
  
ext_usage_str = [usage_str, '\n', ...
                 'Type ', mfilename, ...
                 '(''list'') for list of available filters.', '\n', ...
                 'Type ', mfilename, ...
                 '(''all'') for struct array of all available filters.', '\n' ...
                ];
%%   Check Input Arguments
error(nargerr(mfilename, nargin, [1:1], nargout, [0:1], ...
              1, sprintf(ext_usage_str), 'struct'));

%% Call wtfilter
try
  wtf = wtfilter(wtfname, 'DWT');
catch
  rethrow(lasterror);
end

return

