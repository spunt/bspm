function setup
% 
% FUNCTION:		setup -- compiles mex files and adds paths to file
%
% USAGE:        setup
%
% EXAMPLE:      setup
%
% AUTHOR:       Ameera X Patel
% CREATED:      03-02-2014
%
% CREATED IN:   MATLAB 7.13
%
% REVISION:     4
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: setup.m 4 04-02-2014 BWTv1.1 axpatel

%% add paths and compile mex files

[bwpath,x1,x2] = fileparts(mfilename('fullpath'));

addpath([bwpath filesep 'BWT']);
addpath([bwpath filesep 'third_party/cprintf']);
addpath([bwpath filesep 'third_party/NIfTI']);
addpath([bwpath filesep 'third_party/wmtsa/dwt']);
addpath([bwpath filesep 'third_party/wmtsa/utils']);

basedir=pwd;

cd([bwpath filesep 'third_party/wmtsa/dwt']);

fprintf('\nCompiling...');

mex('modwtj.c');
mex('imodwtj.c');

if not(isempty(which('modwtj'))) && not(isempty(which('imodwtj')))
	fprintf('\nInstallation complete\n');
end

cd(basedir)

clear bwpath x1 x2 basedir

end
