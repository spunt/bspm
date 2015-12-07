function X = imodwt(WJt, VJt, att)
% imodwt -- Calculate the inverse (partial) maximal overlap discrete wavelet transform (IMODWT).
%
%****f* wmtsa.dwt.modwt/imodwt
%
% NAME
%   imodwt -- Calculate the inverse (partial) maximal overlap discrete wavelet transform (IMODWT).
%
% USAGE
%   X = imodwt(WJt, VJt, att)
%
% INPUTS
%   * WJt        -- MODWT wavelet coefficents (N x J x NChan array).
%   * VJt        -- MODWT scaling coefficients (N x {1,J} x NChan vector).
%   * att        -- MODWT transform attributes (struct).
%
% OUTPUT
%   * X          -- reconstituted set of observations (vector).
%
% DESCRIPTION
%   imodwt computes the reconstituted time series from the MODWT wavelet
%   and scaling coefficients.
%  
%
% EXAMPLE
%   X = imodwt(WJt, VJt, att);
%
% ERRORS  
%   WMTSA:InvalidNumArguments       =  'Invalid number of arguments specified in function call'
%   WMTSA:InvalidWavelet            =  'Invalid wavelet filter specified'
%
% NOTES
%   1. Tests indicate othat riginal and reconstituted time series agree within a 
%      precision of 10^-11, which is larger than the 10^15 numeric precision.
%  
%   2. If the MODWT was calculated using 'reflection' boundary
%      conditions, which extends the time series, and the  computed coefficients, 
%      have been truncated to the length of the original series.  The inverse MODWT 
%      will yield a reconstituted time series that differs from the original.  
%      If using 'reflection' boundary conditions, compute the MODWT with the 
%      opt.TruncateCoefs = 0 option (the default) to compute IMODWT which
%      yields and exact replica of the original.
% 
% ALGORITHM
%   See pages 177-179 of WMTSA for description of Pyramid Algorithm for
%   the inverse MODWT.
%
% REFERENCES
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%   Time Series Analysis. Cambridge: Cambridge University Press.
%
% SEE ALSO
%   modwt, imodwtj
%
% TOOLBOX
%   wmtsa/dwt
%
% CATEGORY
%   modwt
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2003-05-01
%
% COPYRIGHT
%   (c) 2003, 2004, 2005 Charles R. Cornish
%
% REVISION
%   $Revision: 632 $
%
%***

% $Id: imodwt.m 632 2006-08-02 06:15:14Z ccornish $

usage_str = ['Usage:  [X] = ', mfilename, ...
             ' (WJt, VJt, att)'];
  
%% Check input arguments and set defaults.
error(nargerr(mfilename, nargin, [3:3], nargout, [0:1], 1, usage_str, 'struct'));

%% Get the wt filter coefficients.
if (isstruct(att.WTF))
  wtf_s = att.WTF;
  ht = att.h;
  gt = att.g;
elseif (ischar(att.WTF))
  wtf_s = modwt_filter(att.WTF);
  ht = wtf_s.h;
  gt = wtf_s.g;
end

N             = att.NX;
NW            = att.NW;
J0            = att.J0;
NChan         = att.NChan;
boundary      = att.Boundary;

%% Pre-allocate memory.
Vout = NaN([N,NChan]);
X    = NaN([N,NChan]); 

%% Do the IMODWT.
for (i = 1:NChan)
  if (att.RetainVJ == 1)
    Vin = VJt(:,J0,i);
  else
    Vin = VJt(:,1,i);
  end

  for (j = J0:-1:1)
    Vout = imodwtj(WJt(:,j,i), Vin, ht, gt, j);
    Vin = Vout; 
  end
  
  X(:,i) = Vout(1:N);
    
end
    
return
