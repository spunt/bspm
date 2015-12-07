% modwtj -- Compute MODWT coefficients for jth level.
%
%****f* wmtsa.dwt/modwtj
%
% NAME
%   modwtj -- Compute MODWT coefficients for jth level.
%
% USAGE
%   [Wt_j, Vout] = modwtj(Vin, ht, gt, j);
%
% INPUTS
%   Vin             - Initial time series, or scaling coefficients for j-1 level.
%   ht              - MODWT avelet filter coefficients.
%   gt              - MODWT Scaling filter coefficients.
%   j               - Level of decomposition.
%
% OUTPUTS
%   Wt_j            - MODWT wavelet coefficients for jth level.
%   Vout            - MODWT scaling coefficients (residuals) for jth level.
%
% DESCRIPTION
%   modwtj is a Mex-Function written in C, which implements the Pyramid Alogrithm
%   of the MODWT for the jth level.  It is usually used as an internal function
%   called by modwt.
%
%   To compile, type:  mex modwtj.c
%
% EXAMPLE
%   [Wt_j, Vout] = modwtj(Vin, ht, gt, j);
%
% ALGORITHM
%   See pages 177-178 of WMTSA for description of Pyramid Algorithm for
%   the MODWT.
%
% REFERENCES
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%   Time Series Analysis. Cambridge: Cambridge University Press.
%
% SEE ALSO
%   modwt
%

% AUTHOR
%   Charlie Cornish
%   Brandon Whitcher
%
% CREATION DATE
%   2003-05-01
%
% COPYRIGHT
%
%
% CREDITS
%   Based on the original function (modwt.c) by Brandon Whitcher.
%
% REVISION
%   $Revision: 612 $
%
%***

%   $Id: modwtj.m 612 2005-10-28 21:42:24Z ccornish $

