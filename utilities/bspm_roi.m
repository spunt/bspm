function roi = bspm_roi(searchpat)
% BSPM_TEMPLATE

% ------- Copyright (C) 2014 -------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

bspmdir = fullfile(getenv('HOME'), 'GitHub', 'bspm', 'imagedata', 'rois');
if nargin==0
    opt = files(fullfile(bspmdir, '*'), 'fileonly', 1);
    roi = uicellect(opt);
end
