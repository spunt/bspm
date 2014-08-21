function matlabbatch = bspm_slicetime2(input)
% BSPM_SLICETIME2
%
%   input structure must contain following fields
%       .epipat
%       .nslices
%       .TR
%       .slice_order
%       .reference_slice
%

% -------------- Copyright (C) 2014 --------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, error('No input!'); end
fn = {'epipat' 'nslices' 'TR' 'slice_order' 'reference_slice'};
[status, msg] = checkfields(input, fn);
if ~status, error(msg); end

% pull out vars
epi_images = files(input.epipat);
nslices = input.nslices;
TR = input.TR;
slice_order = input.slice_order;
reference_slice = input.reference_slice;

% make sure image names are cell arrays of strings
if ischar(epi_images), epi_images = cellstr(epi_images); end

% fix end of image filename cell array
for i = 1:length(epi_images), epi_images(i) = cellstr([epi_images{i} ',1']); end

% build job variable
matlabbatch{1}.spm.temporal.st.scans{1} = cellstr(epi_images);
matlabbatch{1}.spm.temporal.st.nslices = nslices;
matlabbatch{1}.spm.temporal.st.tr = TR;
matlabbatch{1}.spm.temporal.st.ta = TR-(TR/nslices);
matlabbatch{1}.spm.temporal.st.so = slice_order;
matlabbatch{1}.spm.temporal.st.refslice = reference_slice;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

% run job
spm('defaults','fmri');      
spm_jobman('run',matlabbatch);

end

 
 
 
 
