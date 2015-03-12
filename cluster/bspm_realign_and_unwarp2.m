function bspm_realign_and_unwarp2(input)
% BSPM_REALIGN_AND_UNWARP
%
%   input fields
%       .epipat
%       .phasepat
%

% ----------- Copyright (C) 2014 -----------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, error('No input!'); end
if ~isfield(input, 'phasepat'), input.phasepat = []; end
fn = {'epipat' 'phasepat'};
[status, msg] = checkfields(input, fn);
if ~status, error(msg); end
if ischar(input.epipat), input.epipat = cellstr(input.epipat); end
nsess = length(input.epipat);
phaseflag = 1; if isempty(input.phasepat), phaseflag = 0; end
for s = 1:nsess
    

    c_epi_images = files(input.epipat{s});
    if phaseflag
        c_phase_map = files(input.phasepat{s});
    end
    for i = 1:length(c_epi_images)
        c_epi_images(i) = cellstr([c_epi_images{i} ',1']);
    end
    if phaseflag
        c_phase_map = cellstr([char(c_phase_map) ',1']);
    end
    % put images into job
    matlabbatch{1}.spm.spatial.realignunwarp.data(s).scans = cellstr(c_epi_images);           % EPI images
    if phaseflag
        matlabbatch{1}.spm.spatial.realignunwarp.data(s).pmscan = cellstr(c_phase_map);        % Phase Map (vdm* file)
    else
        matlabbatch{1}.spm.spatial.realignunwarp.data(s).pmscan = '';
    end
    
end

% session non-specific paramters
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 4;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

% run job
spm('defaults','fmri');     
spm_jobman('run',matlabbatch);

end

 
 
 
 
