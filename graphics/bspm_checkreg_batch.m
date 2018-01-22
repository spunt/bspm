function bspm_checkreg_batch(in, addmniref)
% BSPM_CHECKREG_BATCH Wrapper for Checking Registration
%
% USAGE: bspm_checkreg_batch(in)
%
% ARGUMENTS
%   in: an array of cells, with each cell containing paths for images
%   to loop over
%

% ------------------------ Copyright (C) 2014 ------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1, mfile_showhelp; return; end
if nargin<2, addmniref = 0; end
if ~iscell(in), in = cellstr(in); end
nim = length(in);
for i = 1:nim
    bspm_checkreg(in{i}, addmniref)
    input(sprintf('%d of %d -- Press any key to move on.', i, nim));
end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
