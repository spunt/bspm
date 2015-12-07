function [WJt, VJt, att] = modwt(X, wtf, nlevels, boundary, varargin)
%   modwt -- Compute the (partial) maximal overlap discrete wavelet transform (MODWT).
%
%****f* wmtsa.dwt.modwt/modwt
%
% NAME
%   modwt -- Compute the (partial) maximal overlap discrete wavelet transform (MODWT).
%
% SYNOPSIS
%   [WJt, VJt, att] = modwt(X, [wtf], [nlevels], [boundary], [{opts}])
%
% INPUTS
%   * X          -- set of observations 
%                   (vector of length NX or matrix of size NX x Nchan)
%   * wtf        -- (optional) wavelet transform filter name or struct 
%                   (string, case-insensitve or wtf struct).
%                   Default:  'la8'
%   * nlevels    -- (optional) maximum level J0 (integer) 
%                   or method of calculating J0 (string).
%                   Valid values: integer>0 or a valid method name
%                   Default:  'conservative'
%   * boundary   -- (optional) boundary conditions to use (string)
%                   Valid values: 'circular' or 'reflection'
%                   Default: 'reflection'
%   * opts       -- (optional) Additional function options.
%
% OUTPUTS
%   * WJt        -- MODWT wavelet coefficents (NW x J x NChan array).
%   * VJt        -- MODWT scaling coefficients (NW x {1,J} x NChan vector).
%   * att        -- MODWT transform attributes (struct).
%
% SIDE EFFECTS
%   1.  wtf is either a string containing a WMTSA-supported MODWT wavelet filter
%       or a valid wtf struct; otherwise error.
%   2.  nlevels is an integer > 0, or is a string containing valid method for
%       choosing J0; otherwise error.
%
% DESCRIPTION
%   modwt calculates the wavelet and scaling coefficients using the maximal
%   overlap discrete wavelet transform (MODWT).
%
%   The optional input arguments have default values:
%   * wtf      -- 'la8' filter
%   * nlevels  -- 'convservative' --> J0 < log2( N / (L-1) + 1)
%   * boundary -- 'reflection'.
%
%   Optional input arguments are specified as name-value pairs:
%   * RetainVJ -- Boolean flag indicating whether to scaling coefficients
%                   have been retained at all levels.
%                   Values:  1 = true,  all VJ retained, 
%                            0 = false, VJ retained for J0 level.
%                   Default: 0, VJ retained only at J0 level.
%
%   The output argument att is a structure with the following fields:
%   * Transform  -- name of transform ('MODWT')
%   * WTF        -- name of wavelet transform filter or a wtf_s struct.
%   * NX         -- number of observations in original series (= length(X))
%   * NW         -- number of wavelet coefficients
%   * J0         -- number of levels of partial decompsition.
%   * NChan      -- number of channels in a multivariate dataset.
%   * Boundary   -- boundary conditions applied.
%   * Aligned    -- Boolean flag indicating whether coefficients are aligned
%                   with original series (1 = true) or not (0 = false).
%   * RetainVJ -- Boolean flag indicating whether VJ scaling coefficients
%                   at all levels have been retained (1= true) or not (0 = false).
%
% EXAMPLE
%   load_ecg;
%   [WJt, VJt, att] = modwt(ecg, 'la8', 6, 'reflection');
%
% WARNINGS
%   WMTSA:MODWT:LargeJ0  =  'MODWT JO > log2(Number of samples).'
%
% ERRORS  
%   WMTSA:invalidBoundary           =  'Invalid Transform boundary method.'
%
% NOTES
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
%   modwtj, modwt_filter, modwt_choose_nlevels, nargerr, argterr
%
% TOOLBOX
%   wmtsa/wmtsa
%
% CATEGORY
%   Transforms: MODWT
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2003-04-23
%
% COPYRIGHT
%   (c) 2003, 2004, 2005 Charles R. Cornish
%
% CREDITS
%   Based on the original function (modwt_dbp.m) by Brandon Whitcher.
%
% REVISION
%   $Revision: 630 $
%
%***

%   $Id: modwt.m 630 2006-05-02 20:47:17Z ccornish $

%%
%% Purpose:  Compute the (partial) maximal overlap discrete wavelet 
%%           transform
%% -------------------------------------------------------------------------
%% Reference: Percival and Guttorp (1994).  Long-memory processes, the 
%%            Allan variance and wavelets.  In "Wavelets and Geophysics,"
%%            pages 325-344.  Academic Press, Inc, San Diego.
%%            Percival and Mofjeld (1997).  Analysis of Subtidal Coastal 
%%            Sea Level Fluctuations Using Wavelets.  Journal of the 
%%            American Statistical Association 92, pp. 868-880.
%%
%% Input: X         Vector of observations
%%        wavelet   Character string; 'haar', 'd4', 'la8', 'la16'
%%        J0   Level of partial MODWT
%%        boundary  Character string; 'circular' or 'reflection'
%%
%% Output : C  Matrix of MODWT wavelet coefficients
%%

%% Defaults
  
defaults.wtf = 'la8';
defaults.nlevels  = 'conservative';
defaults.boundary = 'reflection';

opts_defaults.RetainVJ = 0;
  
usage_str = ['Usage:  [WJt, VJt, att] = ', mfilename, ...
             '(X, [wtf], [nlevels], [boundary], [opts])'];
  
%%  Check input arguments and set defaults.
error(nargerr(mfilename, nargin, '1:', nargout, [0:3], 1, usage_str, 'struct'));

set_defaults(defaults);

%% Parse and check the options.
opts = parse_opts(varargin{:});
error(validate_opts(opts, opts_defaults, 'struct'));
opts = set_opts_defaults(opts, opts_defaults);

%% Get a valid wavelet transform filter coefficients struct.
if (ischar(wtf))
  try
    [wtf_s] = modwt_filter(wtf);
  catch
    rethrow(lasterror);
  end
elseif (iswtf(wtf))
  wtf_s = wtf;
else
  error('WMTSA:invalidWaveletTransformFilter', ...
        encode_errmsg('WMTSA:invalidWaveletTransformFilter', wmtsa_err_table, 'wtf'));
end
  
wtfname = wtf_s.Name;
gt = wtf_s.g;
ht = wtf_s.h;

% If a vector, make X a column vector
if (wmtsa_isvector(X, 'nonsingleton'))
  X = X(:);
end

%%  NX    = length of original series
%%  NChan = number of channels for multi-variant dataset.
[NX, NChan] = size(X);

%%  If nlevels is an integer > 0, set J0 = nlevels.
%%  otherwise, select J0 based on choice method specified.
if (isa(nlevels, 'char'))
  J0 = modwt_choose_nlevels(nlevels, wtfname, NX);
elseif (isnumeric(nlevels))
  if (nlevels > 0)
    J0 = nlevels;
  else
    error('WMTSA:negativeJ0', ...
          ['nlevels must be an integer greater than 0.']);
  end
else
  error('WMTSA:invalidNLevelsValue', ...
         encode_errmsg('WMTSA:invalidNLevelsValue', wmtsa_err_table));
end


if (J0 < 0)
  error('WMTSA:negativeJ0', ...
        ['J0 must be greater than 0.']);
end

if (2.^J0 > NX)
  warning('WMTSA:MODWT:LargeJ0', 'MODWT JO > log2(Number of samples).');
end

%% Initialize the scale (Vin) for first level by setting it equal to X
%% using specified  boundary conditions
switch boundary
  case 'reflection'
   Xin = cat(1, X, flipdim(X, 1));
  case {'circular', 'periodic'}
   Xin = X;
  otherwise
   error('WMTSA:invalidBoundary', ['Invalid boundary method.']);
end

%% NW = length of the extended series = number of coefficients
NW = size(Xin, 1);


%% Pre-allocate memory.
WJt = NaN([NW, J0, NChan], 'double');
if (opts.RetainVJ)
  VJt = NaN([NW, J0, NChan], 'double');
else
  VJt = NaN([NW, 1, NChan], 'double');
end

%% Do the MODWT.
for (i = 1:NChan)
  Vin = Xin(:,i);
  for (j = 1:J0)
    [Wt_j, Vout] = modwtj(Vin, ht, gt, j);
    WJt(:,j,i) = Wt_j;
    Vin = Vout;
    if (opts.RetainVJ)
      VJt(:,j,i) = Vout;
    end
  end
  if (~opts.RetainVJ)
    VJt(:,1,i) = Vout;
  end
end

%% Update attributes
att.Transform = 'MODWT';
if (isstruct(wtf))
  att.WTF = wtf;
else
  att.WTF = wtfname;
end
att.NX         = NX;
att.NW         = NW;
att.J0         = J0;
att.NChan      = NChan;
att.Boundary   = boundary;
att.Aligned    = 0;
att.RetainVJ = opts.RetainVJ;

return
