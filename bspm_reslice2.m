function status = bspm_reslice2(input, runtag)
% BSPM_RESLICE2  Just a wrapper for re-slicing images
%
% USAGE: bspm_reslice2(input, runtag)
%
%   input.reference
%   input.source
%   input.method
%   input.prefix
%   input.checkorient
%
% Options for "method" are:
%
%   0 - Nearest neighbor
%   1 - Trilinear
%   2:7 - 2nd through 7th degree B-Spline
%

% --------------------------------------- Copyright (C) 2014 ---------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin < 2, runtag = 1; end
if nargin < 1, error('No input!'); end
if ~isfield(input, 'method'), input.method = 1; end
if ~isfield(input, 'prefix'), input.prefix = 'r'; end
if ~isfield(input, 'checkorient'), input.checkorient = 0; end
fn = {'reference' 'source' 'method' 'prefix' 'checkorient'};
[status, msg] = checkfields(input, fn);
if ~status, error(msg); end
reference = input.reference; 
source = input.source; 
method = input.method; 
prefix = input.prefix; 
if ischar(reference); reference = cellstr(reference); end
if ischar(source); source = cellstr(source); end
status = 1; 
if input.checkorient
    [status, info] = bspm_check_orient([reference; source], 0);
    if status, return; end
end
    
% run job         
matlabbatch{1}.spm.spatial.coreg.write.ref = reference;
matlabbatch{1}.spm.spatial.coreg.write.source = source;
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = method;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = prefix;

% run job
if runtag
    spm('defaults','fmri');
    spm_jobman('run',matlabbatch);
end
end

function [sts, info] = bspm_check_orient(V, verbose)
% Check the dimensions and orientations of the images
% FORMAT [sts, str] = spm_check_orientations(V [,verbose])
% V       - a struct array as returned by spm_vol
% verbose - [Default: true]
%
% sts     - status (true means OK)
% str     - string describing status, empty if OK
%
% When used without LHS, this function throws an error accordingly.
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_check_orientations.m 5097 2012-12-06 16:08:16Z guillaume $
if nargin < 1, error('USAGE: [sts, str] = bspm_check_orientations(images, verbose)'); end
if nargin < 2, verbose = true; end
if ~isstruct(V)
    if iscell(V), V = char(V); end
    try V = spm_vol(V); catch err, rethrow err; end
end
sts = true;
if nargout==2
    info = [{V.fname}' {V.dim}' {V.mat}']; 
end
dims = cat(1,V.dim);
if any(any(diff(dims,1,1),1))
    sts = false;
    if verbose, disp('The images do not all have the same dimensions.'); end
    info.im 
end
matx = reshape(cat(3,V.mat),[16,numel(V)]);
if any(any(abs(diff(matx,1,2))>1e-4))
    sts = false;
    if verbose, disp('The images do not all have same orientation and/or voxel sizes.'); end
end
end
 
