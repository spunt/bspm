function bspm_init(allowGraphics)
% BSPM_INIT Initialize SPM
%
%   USAGE: bspm_init
%       
%       in  =  array of images OR wildcard pattern for finding them
%       rmtag = flag to delete old image, 0=No (default), 1=Yes 
%

% ----------------------- Copyright (C) 2014 -----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin==0, allowGraphics = usejava('Desktop'); end
spm defaults fmri
spm_jobman initcfg
spm_get_defaults('cmdline', ~allowGraphics)