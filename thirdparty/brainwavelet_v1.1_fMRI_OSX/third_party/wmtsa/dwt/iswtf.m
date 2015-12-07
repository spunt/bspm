function [tf] = iswtf(wtf)
% iswtf -- Determine if input is a valid wtf struct.
%
%****f* wmtsa.dwt/iswtf
%
% NAME
%   iswtf -- Determine if input is a valid wtf struct.
%
% SYNOPSIS
%   [tf] = iswtf(wtf)
%
% INPUTS
%   * wtf        -- wavelet tranform filter struct (wtf_s).
%
% OUTPUTS
%   * tf         -- flag indicating whether a valid valid wtf struct (Boolean)
%
% SIDE EFFECTS
%   Function call requires a minimum of 1 input arguments; otherwise error.
%
% DESCRIPTION
%   iswtf determines whether the input argument is a valid wtf struct, i.e.
%   having the wtf_s struct fields:
%   * g         -- scaling (low-pass) filter coefficients (vector).
%   * h         -- wavelet (high-pass) filter coefficients (vector).
%   * L         -- filter length (= number of coefficients) (integer).
%   * name      -- name of wavelet filter (character string).
%   * wtfclass  -- class of wavelet filters (character string).
%   * transform -- name of transform (character string).
%
% EXAMPLE
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
%
%
% REFERENCES
%
%
% SEE ALSO
%
%
% TOOLBOX
%   wmtsa/dwt
%
% CATEGORY
%
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%
%
% COPYRIGHT
%   (c) 2005 Charles R. Cornish
%
% MATLAB VERSION
%   7.0
%
% CREDITS
%
%
% REVISION
%   $Revision: 612 $
%
%***

%   $Id: iswtf.m 612 2005-10-28 21:42:24Z ccornish $

%% Set Defaults
  
usage_str = ['[tf] = ', mfilename, '(wtf)'];
  
%% Check arguments.
error(nargerr(mfilename, nargin, [1:1], nargout, [0:1], 1, usage_str, 'struct'));

wtf_s = wtf_s_new;

tf = compare_strut_fieldnames(wtf_s, wtf);

return

  
      
function wtf = wtf_s_new
% wft_new -- Create a new (empty) wtf struct.
  wtf = struct( ...
      'Name', {}, ...
      'g', {}, ...
      'h', {}, ...
      'L', {}, ...
      'WTFClass', {}, ...
      'Transform', {} ...
  );
return
