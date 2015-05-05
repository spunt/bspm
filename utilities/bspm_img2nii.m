function bspm_img2nii(in, varargin)
% BSPM_IMG2NII
%
%   USAGE: bspm_img2nii(in, varargin)
%

% ----------------------- Copyright (C) 2014 -----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
def = { ... 
	'keeporiginal',     0,  ...
	'compress',         0,  ...
	};
vals = setargs(def, varargin);
if nargin < 1, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
if ~iscell(in) && strfind(in,'*'); in = files(in); end
if ischar(in), in = cellstr(in); end
fprintf('\nIMG -> NII for %d files\n', length(in)); 
for i = 1:length(in)
    
    [pat,nam,ext] = fileparts(in{i});
    if ~exist(fullfile(pat, [nam '.hdr']))
        fprintf('%s: NO HDR FILE FOUND, DELETING: %s\n', printcount(i, length(in)), in{i});
        delete(in{i});
        continue
    end
    fprintf('%s: %s\n', printcount(i, length(in)), in{i});
    nii = nii_tool('load', in{i});
    outname = fullfile(pat, [nam '.nii']); 
    if compress, outname = [outname '.gz']; end
    nii_tool('save', nii, outname);
    if ~keeporiginal
        delete(in{i}); 
        delete(fullfile(pat, [nam '.hdr'])); 
    end
    
end

 
 
 
 
