function bspm_reorient_batch(imdirs)
% BSPM_REORIENT_BATCH
%
%   USAGE: bspm_reorient_batch(imdirs)
%
%   ARGUMENTS
%       imdirs: directories to start from to select images to display
%       
% Created April 14, 2013 - Bob Spunt

% ------------------------ Copyright (C) 2014 ------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, mfile_showhelp; return; end
if ~ischar(imdirs), imdirs = cellstr(imdirs); end
basedir = pwd;
ndir = length(imdirs);
for i = 1:ndir
    if exist(imdirs{i}, 'dir')
        cdir = imdirs{i};
        cd(cdir);
        bspm_display;
    else
        cdir = fileparts(imdirs{i});
        cd(cdir);
        bspm_display(imdirs{i});
    end
    input(sprintf('%d of %d -- Press any key to move on.', i, ndir));
end
cd(basedir)
 
 
 
 
