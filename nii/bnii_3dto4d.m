function bnii_3dto4d(fname3d, varargin)
% BNII_3DTO4D_FILENAME Find characters in common across multiple strings
%
%  USAGE: bnii_3dto4d(fname3d, varargin)
% __________________________________________________________________________
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-04-03
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
def = { ... 
    'compress',         0,  ...
    'delete3d',         0,  ...
    'delimiter',        '_' ...
	};
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
if ischar(fname3d), fname3d = cellstr(fname3d); end
fname4d = bnii_3dto4d_filename(fname3d, 'addgz', compress, 'delimiter', '_');
fprintf('| Combing %d 3D images to: %s', length(fname3d), fname4d); 
nii4d   = nii_tool('cat3D', fname3d);
nii_tool('save', nii4d, fname4d);
if exist(fname4d, 'file') && delete3d
    fprintf('\n| Deleting 3D images'); 
    rm(fname3d, 'verbose', 0); 
end
fprintf('\n - Process Complete -\n');