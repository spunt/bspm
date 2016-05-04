function matlabbatch = bspm_reorient2acpc(images, opt)
% BSPM_REORIENT2ACPC
%
% USAGE:  matlabbatch = bspm_reorient2acpc(images, opt)
%
% ARGUMENTS
%   images
%   opt: 1 for T1, 2 for T2
%

% ---------------------- Copyright (C) 2014 ----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<2, mfile_showhelp; return; end
if ischar(images), images = cellstr(images); end
images = strcat(images, ',1');
spmdir = fileparts(which('spm'));
if opt==1
    ref = fullfile(spmdir, 'templates', 'T1.nii,1');
    for i = 1:length(images)
        matlabbatch{i}.spm.spatial.coreg.estimate.ref =     cellstr(ref);
        matlabbatch{i}.spm.spatial.coreg.estimate.source =  cellstr(images(i));
        matlabbatch{i}.spm.spatial.coreg.estimate.other =   {''};
        matlabbatch{i}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{i}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{i}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{i}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    end
else
    ref = fullfile(spmdir, 'templates', 'EPI.nii,1');
    if length(images) > 1, other = cellstr(images(2:end)); else other = {''}; end
    matlabbatch{1}.spm.spatial.coreg.estimate.ref =     cellstr(ref);
    matlabbatch{1}.spm.spatial.coreg.estimate.source =  cellstr(images(1));
    matlabbatch{1}.spm.spatial.coreg.estimate.other =   other;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
end
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end
end
