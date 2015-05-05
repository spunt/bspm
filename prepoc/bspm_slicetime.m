function matlabbatch = bspm_slicetime(in, dcmfile, opt)
% BSPM_SLICETIME
%
%   USAGE: matlabbatch = bspm_slicetime(in, dcmfile, opt)
%
%   opt:    1 = enter as slice order
%           2 = enter as slice times (e.g., for multi-band acquisitions)
%

% -------------------------------- Copyright (C) 2014 --------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<2, disp('USAGE: matlabbatch = bspm_slicetime(in, dcmfile, opt)'); end
if nargin<3, opt = 2; end
if ischar(in), in = cellstr(in); end
if iscell(dcmfile), dcmfile = char(dcmfile); end
if length(in)==1 % - assume 4D
    in = bspm_expand4D(in); 
else
    in = cellstr(strcat(in, ',1'));
end
hdr         = spm_dicom_headers(dcmfile); 
TR          = hdr{1}.RepetitionTime/1000; 
slice_time  = hdr{1}.Private_0019_1029';
nslices     = length(slice_time);
slice_order = [slice_time (1:length(slice_time))']; 
slice_order = sortrows(slice_order,1); 
slice_order(:,1) = [];
order{1}    = slice_order; 
order{2}    = slice_time; 
ref_slice   = slice_order(round(nslices/2));
ref_time    = slice_time(ref_slice);
ref{1}      = ref_slice; 
ref{2}      = ref_time; 

% | build job variable
matlabbatch{1}.spm.temporal.st.scans{1} = in; 
matlabbatch{1}.spm.temporal.st.nslices = nslices;
matlabbatch{1}.spm.temporal.st.tr = TR;
matlabbatch{1}.spm.temporal.st.ta = TR-(TR/nslices);
matlabbatch{1}.spm.temporal.st.so = order{opt};
matlabbatch{1}.spm.temporal.st.refslice = ref{opt};
matlabbatch{1}.spm.temporal.st.prefix = 'a';

% | run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end

end
