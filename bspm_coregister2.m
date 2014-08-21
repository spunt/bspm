function matlabbatch = bspm_coregister2(input, runtag)
% BSPM_COREGISTER2
%
% USAGE:  bspm_coregister2(input)
%
% FIELDS
%   reference = image that stays put
%   source = image that is moved to match the reference
%   other = images to move with the source
%

% ------------------ Copyright (C) 2014 ------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 2, runtag = 1; end
if nargin < 1, error('No input!'); end
if ~isfield(input, 'other'), input.other = {''}; end
fn = {'reference' 'source' 'other'};
[status, msg] = checkfields(input, fn);
if ~status, error(msg); end
reference = input.reference;
source = input.source;
other = input.other;


% make sure images are in cell arrays
if ischar(reference), reference = cellstr(reference); end
if ischar(source), source = cellstr(source); end
if ischar(other), other = cellstr(other); end

%  add ,1 to end of image filenames
reference = cellstr([char(reference) ',1']);
source = cellstr([char(source) ',1']);
if ~isempty(other{1})
    for i = 1:length(other)
        other(i) = cellstr([other{i} ',1']);
    end
end

% build job
matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(reference);
matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(source);
matlabbatch{1}.spm.spatial.coreg.estimate.other = cellstr(other);
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
 
% run job
if runtag
    spm('defaults','fmri');
    spm_jobman('run',matlabbatch);
end
end
 

 
 
 
 
