function matlabbatch = bspm_deformations(in, field, varargin)
% BSPM_DEFORMATIONS
%
%   USAGE: matlabbatch = bspm_deformations(in, field, varargin)
%
%   ARGUMENTS:
%       in      = images to deform
%       field   = deformation field, e.g.
%               iy_*    = forward (native -> MNI)
%               y_*     = backward (MNI -> native)
%               
%

% --------------------- Copyright (C) 2014 ---------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
def = { 'weightimage',      [], ...
        'fwhm',             [0 0 0], ...
        'voxsize',          [1 1 1] ...
        };
vals = setargs(def, varargin);
if nargin < 2, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
if ischar(in), in = cellstr(in); end
if iscell(field), field = char(field); end
if isempty(weightimage), weightimage = {''}; end
if ischar(weightimage), weightimage = cellstr(weightimage); end
if length(fwhm)==1, fwhm = repmat(fwhm,1,3); end
if length(voxsize)==1, voxsize = repmat(voxsize,1,3); end

% | JOB
matlabbatch{1}.spm.util.defs.comp{1}.def                    = cellstr(field); 
matlabbatch{1}.spm.util.defs.out{1}.push.fnames             = in; 
matlabbatch{1}.spm.util.defs.out{1}.push.weight             = weightimage;
matlabbatch{1}.spm.util.defs.out{1}.push.savedir.savesrc    = 1;
matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.bb       = [-78 -112 -50; 78 76 85];
matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.vox      = voxsize;
matlabbatch{1}.spm.util.defs.out{1}.push.preserve           = 0;
matlabbatch{1}.spm.util.defs.out{1}.push.fwhm               = fwhm;

% | RUN
if nargout==0, bspm_runbatch(matlabbatch); end

end
