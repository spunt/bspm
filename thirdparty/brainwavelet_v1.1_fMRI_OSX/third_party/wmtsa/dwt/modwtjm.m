function [Wtout, Vtout] = modwtjm(Vtin, ht, gt, j)
% modwtjm -- Calculate jth level MODWT coefficients (MATLAB implementation).
%
%****f* wmtsa.dwt/modwtjm
%
% NAME
%   modwtjm -- Calculate jth level MODWT coefficients (MATLAB implementation).
%
% SYNOPSIS
%   [Wtout, Vtout] = modwtjm(Vtin, ht, gt, j)
%
% INPUTS
%   * Vtin        -- Input series for j-1 level (i.e. MODWT scaling coefficients) 
%   * ht          -- MODWT wavelet filter coefficients.
%   * gt          -- MODWT scaling filter coefficients.
%   * j           -- level (index) of scale.
%
% OUTPUTS
%   * Wtout       -- MODWT wavelet coefficients for jth scale.
%   * Vtout       -- MODWT scaling coefficients for jth scale.
%
% SIDE EFFECTS
%
%
% DESCRIPTION
%   modwtjm is an implementation in MATLAB code of the MODWT transform for 
%   the jth level, and is included in the toolkit for illustrative purposes 
%   to demonstrate the pyramid algothrim.
%
%   For speed considerations, the modwt function uses the C implementation of 
%   the MODWT transform, modwtj, which linked in as a MEX function.
%
% EXAMPLE
%   X = wmtsa_data('ecg');
%   wtf = modwt_filter('la8');
%   % Compute the j = 1 level coefficients for ECG time series.
%   j = 1;
%   [Wtout, Vtout] = modwtjm(X, wft.h, wtf.g, j);
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
% BUGS
%
%
% TODO
%
%
% ALGORITHM
%   See page 177-178 of WMTSA for pyramid algorithm.
%
% REFERENCES
%
%
% SEE ALSO
%   modwtj, modwt, modwt_filter
%
% TOOLBOX
%   wmtsa
%
% CATEGORY
%   dwt
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2005-01-12
%
% COPYRIGHT
%   
%
% CREDITS
%
%
% REVISION
%   $Revision: 612 $
%
%***

%   $Id: modwtjm.m 612 2005-10-28 21:42:24Z ccornish $

  usage_str = ['Usage:  [Wtout, Vtout] = ', mfilename, ...
             '(Vtin, ht, gt, j)'];
  
  %%  Check input arguments and set defaults.
  error(nargerr(mfilename, nargin, [4:4], nargout, [0:2], 1, usage_str, 'struct'));

  N = length(Vtin);
  L = length(ht);
  
  Wtout = zeros(N, 1) * NaN;
  Vtout = zeros(N, 1) * NaN;
  
  for (t = 1:N)
    k = t;
    Wtout(t) = ht(1) * Vtin(k);
    Vtout(t) = gt(1) * Vtin(k);
    
    for (n = 2:L)
      k = k - 2^(j-1);
      if (k < 1)
        k = k + N;
      end
      Wtout(t) = Wtout(t) + ht(n) * Vtin(k);
      Vtout(t) = Vtout(t) + gt(n) * Vtin(k);
    end
  end    
      
return
