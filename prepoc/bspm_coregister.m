function matlabbatch = bspm_coregister(reference, source, other, costfun, checkreg_after)
% BSPM_COREGISTER
%
% USAGE:  matlabbatch = bspm_coregister(reference, source, other, costfun, checkreg_after)
%
% ARGUMENTS
%   reference = image that stays put
%   source = image that is moved to match the reference
%   other = images to move with the source
%   costfun =   1 (default) is normalised mutual information (nmi)
%               2 is entropy correlation coefficient (ecc)
%

% ---------------------- Copyright (C) 2014 ----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<2, mfile_showhelp; return; end
if nargin<3, other = {''}; end
if nargin<4, costfun = 1; end
if nargin<5, checkreg_after = 0; end

% make sure images are in cell arrays
if ischar(reference), reference = cellstr(reference); end
if ischar(source), source = cellstr(source); end
if ischar(other), other = cellstr(other); end

%  add ,1 to end of image filenames
if isempty(regexp(char(reference), ',\d$', 'ONCE')), reference = cellstr([char(reference) ',1']); end
if isempty(regexp(char(source), ',\d$', 'ONCE')), source = cellstr([char(source) ',1']); end
if ~isempty(other{1})
    for i = 1:length(other)
        if isempty(regexp(char(other(i)), ',\d+$','ONCE')), other(i) = cellstr([char(other(i)) ',1']); end
    end
end
if costfun==1, cost_fun = 'nmi'; elseif costfun==2, cost_fun = 'ecc'; end

% build job
matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(reference);
matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(source);
matlabbatch{1}.spm.spatial.coreg.estimate.other = cellstr(other);
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = cost_fun;
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

% run job
if nargout==0
    spm_jobman('initcfg'); spm_jobman('run',matlabbatch); 
    % check reg
    if checkreg_after
        checkim = [reference; source]; 
        if ~isempty(other), checkim = [checkim; other(1)]; end
        if length(other) > 1, checkim = [checkim; other(end)]; end
        [~,captions] = cellfun(@fileparts, checkim, 'Unif', false); 
        bspm_checkreg(checkim, 0, captions); 
    end
end

end
