function bspm_deface(in)
% BSPM_DEFACE
%
%   USAGE: bspm_deface(in)
%

% ----------------------- Copyright (C) 2014 -----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin < 1, mfile_showhelp; return; end
if ~iscell(in) && strfind(in,'*'); in = files(in); end
if ischar(in), in = cellstr(in); end
fprintf('\nDE-FACING %d files\n', length(in));
job.images = in; 
spm_deface(job);