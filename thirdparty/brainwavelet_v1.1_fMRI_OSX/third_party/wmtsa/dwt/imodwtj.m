% imodwtj -- Compute inverse MODWT for jth level.
%
%****f* wmtsa.dwt/imodwtj
%
% NAME
%   imodwtj -- Compute inverse MODWT for jth level.
%
% SYNOPSIS
%   Vout = imodwtj(Wt_j, Vin, ht, gt, j)
%
% INPUTS
%   Wt_j            - MODWT wavelet coefficients for jth level.
%   Vin             - Scaling coefficients at J0, or j-1 level.
%   ht              - MODWT avelet filter coefficients.
%   gt              - MODWT Scaling filter coefficients.
%   j               - Level of decomposition.
%
% OUTPUTS
%   Vout            - MODWT scaling coefficients (residuals) for jth level.
%
% DESCRIPTION
%   modwtj is a Mex-Function written in C, which implements the Pyramid Alogrithm
%   of the MODWT for the jth level.  It is usually used as an internal function
%   called by imodwt or imodwt_mra.
%
%   To compile, type:  mex modwtj.c
%
% EXAMPLE
%   Vout = imodwtj(Wt_j, Vin, ht, gt, j);
%
% ALGORITHM
%
%
% REFERENCES
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%   Time Series Analysis. Cambridge: Cambridge University Press.
%
% SEE ALSO
%   imodwt, imodwt_mra
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
%   Based on the original function (imodwt.c) by Brandon Whitcher.
%
% REVISION
%   $Revision: 612 $
%
%***

%   $Id: imodwtj.m 612 2005-10-28 21:42:24Z ccornish $

  
