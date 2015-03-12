function matlabbatch = bspm_norm_estimate(source_img, template_img, source_weighting_img, template_weighting_img)
% BSPM_NORM_ESTIMATE
%
%   ARGUMENTS:
%       source_img = image to warp
%       template_img = template image
%       source_weighting_img = (optional) weighting  images  (consisting
%           of  pixel  values  between  the range  of  zero  to  one)  to
%           be used for registering abnormal or lesioned brains. These
%           images should match the dimensions of the image from which
%           the parameters are estimated, and should contain zeros
%           corresponding to regions of abnormal tissue.
%       template_weighting_img = (optional) this should have the same
%           dimensiosn as the template images, with values from 0 to 1.
%

% ---------------------------------------- Copyright (C) 2014 ----------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 2, error('USAGE: bspm_norm_estimate(source_img, template_img, source_weighting_img, template_weighting_img'); end
if nargin < 3, source_weighting_img = {}; template_weighting_img = {}; end
if ischar(source_img), source_img = cellstr(source_img); end
if ischar(template_img), template_img = cellstr(template_img); end

% build job
% -------------------------------------------------------------------------
matlabbatch{1}.spm.spatial.normalise.est.subj.source = strcat(source_img, {',1'}); % source image
matlabbatch{1}.spm.spatial.normalise.est.subj.wtsrc = source_weighting_img; % source weighting image
matlabbatch{1}.spm.spatial.normalise.est.eoptions.template = strcat(template_img, {',1'}); % template image
matlabbatch{1}.spm.spatial.normalise.est.eoptions.weight = template_weighting_img; % template weighting image
matlabbatch{1}.spm.spatial.normalise.est.eoptions.smosrc = 8;       % source smoothing (FWHM)
matlabbatch{1}.spm.spatial.normalise.est.eoptions.smoref = 0;       % template smoothing (FWHM)
matlabbatch{1}.spm.spatial.normalise.est.eoptions.regtype = 'mni';  % target space
matlabbatch{1}.spm.spatial.normalise.est.eoptions.cutoff = 25;      % nonlinear frequency cutoff
matlabbatch{1}.spm.spatial.normalise.est.eoptions.nits = 16;        % num nonlinear iterations
matlabbatch{1}.spm.spatial.normalise.est.eoptions.reg = 1;          % nonlinear regularisation

% run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end


end
