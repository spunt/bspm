function bspm_reslice(reference, sources, method, prefix)
% BSPM_RESLICE  Just a wrapper for re-slicing images
%
% USAGE: bspm_reslice(reference, sources, method, prefix)
%
% Options for "method" are:
%   0   - Nearest neighbor
%   1   - Trilinear (default)
%   2:7 - 2nd through 7th degree B-Spline
%

% --------------------------------------- Copyright (C) 2014 ---------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin < 4, prefix = 'r'; end
if nargin < 3, method = 1; end
if nargin < 2, disp('USAGE: bspm_reslice(reference, sources, method, prefix)'); return; end
if ischar(reference); reference = cellstr(reference); end
if ischar(sources); sources = cellstr(sources); end
P = char([reference; sources]);
flags.mean      = false;    % don't write mean image
flags.which     = 1;        % don't reslice reference image
flags.method    = method;   % interpolation method
flags.prefix    = prefix;
flags.wrap      = [0 0 0]; 
spm_reslice(P, flags); 
end 
 
 
 
