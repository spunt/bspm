function matlabbatch = bspm_slicetime(in, dcmfile)
% BSPM_SLICETIME
%
%   USAGE: matlabbatch = bspm_slicetime(in, dcmfile)
%

% -------------------------------- Copyright (C) 2014 --------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<1, disp('USAGE: matlabbatch = bspm_slicetime(in, dcmfile)'); end
if ischar(in), in = cellstr(in); end
if iscell(dcmfile), dcmfile = char(dcmfile); end
hdr         = spm_dicom_headers(dcmfile); 
TR          = hdr{1}.RepetitionTime/1000; 
slice_time  = hdr{1}.Private_0019_1029';
nslices     = length(slice_time);
slice_order = [slice_time (1:length(slice_time))']; 
slice_order = sortrows(slice_order,1); 
slice_order(:,1) = [];
ref_slice   = slice_order(round(nslices/2));
ref_time    = slice_time(ref_slice); 

% | build job variable
matlabbatch{1}.spm.temporal.st.scans{1} = cellstr(strcat(in, ',1'));
matlabbatch{1}.spm.temporal.st.nslices = nslices;
matlabbatch{1}.spm.temporal.st.tr = TR;
matlabbatch{1}.spm.temporal.st.ta = TR-(TR/nslices);
matlabbatch{1}.spm.temporal.st.so = slice_time;
matlabbatch{1}.spm.temporal.st.refslice = ref_time;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

% | run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end

end
