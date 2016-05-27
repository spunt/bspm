function matlabbatch = bspm_norm_westrite(image2align, images2write, vox_size)
% BSPM_NORM_ESTWRITE
%
%   ARGUMENTS:
%       image2align     = parameter file (*sn.mat)
%       images2write    = images to norm (char or cell array)
%       vox_size        = resolution (in mm) at which to re-sample normed volumes
%

% -------------------------- Copyright (C) 2014 --------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, mfile_showhelp; return; end
if nargin < 3, vox_size = [2 2 2]; end
if length(vox_size)==1, vox_size = [vox_size vox_size vox_size]; end
if ischar(image2align), images2align = cellstr(image2align); end
if ischar(images2write), images2align = cellstr(images2write); end
images2write = bspm_check_filenames(images2write);
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = strcat(image2align, {',1'});
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = images2write;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = cellstr(fullfile(fileparts(which('spm')), 'tpm', 'TPM.nii'));
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70; 78 76 85];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = vox_size;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';

% run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end

end
% ==========================================================================
function mfile_showhelp(varargin)
% MFILE_SHOWHELP
ST = dbstack('-completenames');
if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
eval(sprintf('help %s', ST(2).file));  
end
function fn = bspm_check_filenames(inname)
% BSPM_CHECK_FILENAMES
% 
%   USAGE outname = bspm_check_orientations(inname)
%       imname can be char or cell 
%       output is a cell array
%
if ischar(inname), inname = cellstr(inname); end
inname              = regexprep(inname, ',\d+$', '');
[pcim, ncim, ecim]  = cellfun(@fileparts, inname, 'unif', false);
fn = []; 
for i = 1:length(inname)
    
    if ~ismember(ecim{i}, {'.gz' '.nii' '.img'})
        continue
    end
    
    % | If compressed, try to decompress
    if strcmp(ecim{i}, '.gz')
        try
            pigz(inname{i})
        catch
            gunzip(inname{i});
        end
        inname{i} = cellstr(fullfile(pcim{i}, ncim{i}));  
    end
    
    % | Check to see if 4D
    nii     = nifti(inname{i}); 
    nvol    = size(nii.dat, 4); 
    append  = cellfun(@num2str, num2cell(1:nvol)', 'Unif', false);
    fn      = [fn; strcat(inname{i}, {','}, append)];
end
end