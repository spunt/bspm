function bnii_4dto3d(fname4d, varargin)
% BNII_4DTO3D
%
%  USAGE: bnii_3dto4d(fname3d, varargin)
% __________________________________________________________________________
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-04-03
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
def = { ... 
    'basefname',        [], ...
    'delete4d',         0,  ...
	};
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
if iscell(fname4d), fname4d = char(fname4d); end
[basepath, n] = fileparts(fname4d); 
if isempty(basefname)
    basefname = regexprep(n, '\w4D', '', 'ignorecase');  
end
nii = nii_tool('load', fname4d); 
fprintf('| Saving %d 3D images from: %s', size(nii.img, 4), fname4d);
nii_tool('save', nii, fullfile(basepath, strcat(basefname, '.nii')), 1); 
if exist(fname4d, 'file') && delete4d
    fprintf('\n| Deleting 4D image'); 
    rm(fname4d, 'verbose', 0); 
end
fprintf('\n - Process Complete -\n');