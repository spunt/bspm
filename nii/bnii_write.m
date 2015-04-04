function bnii_write(img, hdr, varargin)
% BNII_WRITE Uses NII_TOOL to write NifTI images
%
%  USAGE: bnii_write(img, hdr, varargin)
% __________________________________________________________________________
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-04-03
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
def = { ... 
	'outname',		'',     ...
	'compress',		0,      ...
	};
vals = setargs(def, varargin);
if nargin < 2, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
if isempty(outname)
    outname = hdr.file_name;
elseif iscell(outname)
    outname = char(outname);
end
[pat,nam,ext] = fileparts(outname);
outname = fullfile(pat, nam, '.nii'); 
if compress, outname = [outname '.gz']; end
nii = struct('img', img, 'hdr', hdr); 
nii_tool('save', nii, outname); 