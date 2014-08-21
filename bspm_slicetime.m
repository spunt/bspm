function [] = bspm_slicetime(epi_images, nslices, TR, slice_order, reference_slice)
% BSPM_SLICETIME
%
%   ARGUMENTS: 
%       1:  epi_images
%       2:  nslices
%       3:  TR
%       4:  slice_order
%       5:  reference_slice
%

% -------------------------------- Copyright (C) 2014 --------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<5
   disp('Give me arguments! bspm_slicetime(epi_images, nslices, TR, slice_order, reference_slice)');
   return
end

% make sure image names are cell arrays of strings
if ischar(epi_images)
    epi_images = cellstr(epi_images);
end

% fix end of image filename cell array
for i = 1:length(epi_images)
    epi_images(i) = cellstr([epi_images{i} ',1']);
end

% build job variable
matlabbatch{1}.spm.temporal.st.scans{1} = cellstr(epi_images);
matlabbatch{1}.spm.temporal.st.nslices = nslices;
matlabbatch{1}.spm.temporal.st.tr = TR;
matlabbatch{1}.spm.temporal.st.ta = TR-(TR/nslices);
matlabbatch{1}.spm.temporal.st.so = slice_order;
matlabbatch{1}.spm.temporal.st.refslice = reference_slice;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

% run job
spm('defaults','fmri'); spm_jobman('initcfg');         
spm_jobman('run',matlabbatch);

end

 
 
 
 
