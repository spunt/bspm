function [] = bspm_smooth(in, fwhm, prefix)
% BSPM_SMOOTH
%
%   USAGE: bspm_smooth(in, fwhm)
%       
%       in  =  array of images to smooth (full path)
%       fwhm  =  smoothing kernel
%       prefix = prefix to add to new filename
%

% --------------- Copyright (C) 2014 ---------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3, prefix = 's'; end
if nargin<2, disp('Give me arguments! bspm_smooth(in, fwhm, prefix)'); return, end
if length(fwhm)==1, fwhm = [fwhm fwhm fwhm]; end

% make sure image names are cell arrays of strings
if ischar(in)
    in = cellstr(in);
end

% add ',1' to images
for i = 1:length(in)
    in(i) = cellstr([in{i} ',1']);
end

% build job variable
matlabbatch{1}.spm.spatial.smooth.data = cellstr(in);
matlabbatch{1}.spm.spatial.smooth.fwhm = fwhm;
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = prefix;

% run job
spm('defaults','fmri'); spm_jobman('initcfg');         
spm_jobman('run',matlabbatch);

end

 
 
 
 
