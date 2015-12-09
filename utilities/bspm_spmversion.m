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
toolpath    = fullfile(getenv('HOME'), 'Documents', 'MATLAB', 'Thirdparty');
rmversion = 8; 
if spmversion==8, rmversion = 12; end
rmpath([toolpath filesep sprintf('spm%d/toolbox/rwls', rmversion)]);
rmpath([toolpath filesep sprintf('spm%d/external/fieldtrip', rmversion)]); 
rmpath(fullfile(toolpath, 'Matlab', sprintf('spm%d', rmversion)));
addpath([toolpath filesep sprintf('spm%d/toolbox/rwls', spmversion)]);
addpath([toolpath filesep sprintf('spm%d/external/fieldtrip', spmversion)]); 
addpath(fullfile(toolpath, 'Matlab', sprintf('spm%d', spmversion)));
 
 
 
 
