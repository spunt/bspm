function fname4d = bnii_3dto4d_filename(fname3d, varargin)
% BNII_3DTO4D_FILENAME Find characters in common across multiple strings
%
%  USAGE: outname = bnii_3dto4d_filename(in, varargin)
% __________________________________________________________________________
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-04-03
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
def = { ... 
    'addgz',            1,  ...
	'delimiter',     '\W',	...
    'omitzeros',        0,  ...
	};
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
if ischar(fname3d), fname3d = cellstr(fname3d); end

% | Find Them! 
% | ========================================================================
[fpath, fname4d] = cellfun(@fileparts, fname3d, 'unif', false);
fname4d     = regexp(fname4d, delimiter, 'split');
nelem       = length(fname4d{1});  
fname4d     = [fname4d{:}]; 
fname4d     = reshape(fname4d, nelem, length(fname3d))';
uname       = zeros(nelem, 1); 
for i = 1:nelem, uname(i) = length(unique(fname4d(:,i))); end
fname4d     = fname4d(1, uname==1);
fname4d     = strcat(fname4d, '_');
fname4d     = strcat(fname4d{:}, sprintf('T%d_4D.nii', length(fname3d)));
if addgz, fname4d = [fname4d '.gz']; end
fname4d     = fullfile(fpath{1}, fname4d); 
