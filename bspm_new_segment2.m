function bspm_new_segment2(input)
% BSPM_NEW_SEGMENT
%
%   ARGUMENTS:
%       input.img = cell array containing paths to all images to segment
%

% ------------------------- Copyright (C) 2014 -------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, error('No input!'); end
fn = {'img'};
[status, msg] = checkfields(input, fn);
if ~status, error(msg); end

in = input.img;

% make sure image names are cell arrays of strings
if ischar(in), in = cellstr(in); end

% fix end of image filename cell array
for i = 1:length(in), in(i) = cellstr([in{i} ',1']); end

% grab spm directory and TPM.nii 
spm_dir = which('spm');
spm_dir = regexprep(spm_dir,'spm.m','');
TPMimg = [spm_dir filesep 'toolbox/Seg/TPM.nii'];

% build job
% -------------------------------------------------                             
matlabbatch{1}.spm.tools.preproc8.channel.vols = in;  % cell array containing paths to all images
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0]; % save [bias_corrected bias_field]
matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[TPMimg ',1']}; % grey matter
matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [0 1];   % [native_space DARTEL_Imported]
matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[TPMimg ',2']}; % white matter
matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [0 1];  % [native_space DARTEL_Imported]
matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[TPMimg ',3']}; % cerebro-spinal fluid (CSF)
matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [0 0];  % [native_space DARTEL_Imported]
matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[TPMimg ',4']}; % bone 
matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [0 0];  % [native_space DARTEL_Imported]
matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[TPMimg ',5']}; % soft tissue
matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [0 0];  % [native_space DARTEL_Imported]
matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[TPMimg ',6']}; % air/background
matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];  % [native_space DARTEL_Imported]
matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.warp.mrf = 0; % MRF parameter (cleanup) (default = 0)
matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
matlabbatch{1}.spm.tools.preproc8.warp.samp = 3; % sampling distance
matlabbatch{1}.spm.tools.preproc8.warp.write = [0 1];  % [inverse_deformation forward_deformation]

% run job
spm('defaults','fmri');    
spm_jobman('run',matlabbatch);

end

 
 
 
 
