function matlabbatch = bspm_realign_and_unwarp(epi_images, phase_map)
% BSPM_REALIGN_AND_UNWARP
%
%  USAGE: matlabbatch = bspm_realign_and_unwarp(epi_images, phase_map)
%

% --------------------- Copyright (C) 2014 ---------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<1, disp('USAGE: matlabbatch = bspm_realign_and_unwarp(epi_images, phase_map)'); return; end
if nargin<2, phase_map = ''; end
if ischar(epi_images{1})
    nsess = 1; charflag = 1;
else
    nsess   = length(epi_images); charflag = 0;
end
phaseflag = 1;
if isempty(phase_map), phaseflag = 0; end
for s = 1:nsess

    if charflag
        c_epi_images = epi_images;
        if phaseflag
            c_phase_map = phase_map;
        end
    else
        c_epi_images = epi_images{s};
        if phaseflag
            c_phase_map = phase_map{s};
        end
    end

    % fix end of image filename cell array
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
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
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
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end


end
