function [] = bspm_norm_write(param, images, vox_size)
% BSPM_NORM_WRITE
%
%   ARGUMENTS:
%       param = parameter file (*sn.mat)
%       images = images to norm (char or cell array)
%       vox_size = resolution (in mm) at which to re-sample normed volumes
%

% -------------------------- Copyright (C) 2014 --------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1
   disp('bspm_norm_write(param, images, vox_size)');
   return
end

% make sure image names are cell arrays of strings
if ischar(param)
    param = cellstr(param);
end
if ischar(images)
    images = cellstr(images);
end

% fix end of image filename cell array
for i = 1:length(images)
    images(i) = cellstr([images{i} ',1']);
end

% build job
% -------------------------------------------------           
matlabbatch{1}.spm.spatial.normalise.write.subj.matname = param;   
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = images;
matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50; 78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [vox_size vox_size vox_size];                  
matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];

% run job
spm('defaults','fmri'); spm_jobman('initcfg');         
spm_jobman('run',matlabbatch);

end

 
 
 
 
