function template = bspm_template(modality)
% BSPM_TEMPLATE
%   template = bspm_template(modality)
%               t1 (default)
%               epi
%

% ------- Copyright (C) 2014 -------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin==0, modality = 't1'; end
imopt   = {'esbrain_avg152T1.nii' 'etEPI.nii'};
if strcmpi(modality, 'epi'), idx = 2; else idx = 1; end
spmdir  = whichdir('spm');
template = [spmdir filesep 'templates' filesep imopt{idx}];
 
 
 
 
