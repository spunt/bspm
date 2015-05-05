function spmversion = bspm_spmversion(spmversion)
% BSPM_SPMVERSION Check/switch SPM version
%
%  USAGE: version = bspm_spmversion(version) 
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-04-06
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
[~, spmdir] = fileparts(fileparts(which('spm.m')));
if nargin==0, spmversion = spmdir; return; end
dropboxpath = '/Users/bobspunt/Desktop/Dropbox/Bob';
rmversion = 8; 
if spmversion==8, rmversion = 12; end
rmpath([dropboxpath filesep sprintf('Matlab/spm%d/toolbox/rwls', rmversion)]);
rmpath([dropboxpath filesep sprintf('Matlab/spm%d/external/fieldtrip', rmversion)]); 
rmpath(fullfile(dropboxpath, 'Matlab', sprintf('spm%d', rmversion)));
addpath([dropboxpath filesep sprintf('Matlab/spm%d/toolbox/rwls', spmversion)]);
addpath([dropboxpath filesep sprintf('Matlab/spm%d/external/fieldtrip', spmversion)]); 
addpath(fullfile(dropboxpath, 'Matlab', sprintf('spm%d', spmversion)));
