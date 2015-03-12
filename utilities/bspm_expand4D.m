function fn = bspm_expand4D(fn4D)
% BSPM4D_LEVEL1
%
% USAGE: fn = bspm_expand4D(fn4D)
%

% -------------------------------------- Copyright (C) 2014 --------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Jul_17_2014
if nargin<1, error('USAGE: fn = bspm_expand4D(fn4D)'); end
if iscell(fn4D), fn4D = char(fn4D); end
hdr = spm_vol(fn4D);
append = cellfun(@num2str, num2cell(1:length(hdr))', 'Unif', false);
fn = strcat(fn4D, {','}, append);
end
 
 
 
 
