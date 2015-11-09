function [spm5 select fparts] = check_spm_setup
%  [spm5 select fparts] = check_spm_setup
% Check for presence and version of spm

if ~exist('spm', 'file')
    error('Can''t find SPM - please add it to you MATLAB path')
end
switch spm('Ver')
    case {'SPM96','SPM99'}
        error('Requires SPM version 2 or above, please upgrade')
    case 'SPM2'
        spm5 = 0;
        select = @(msg) ...
            spm_get(inf, '*.*', msg);
        fparts = @fileparts;
    case {'SPM5', 'SPM8b', 'SPM8'} % SPM8 assumed compatible with SPM5
        spm5 = 1;
        select = @(msg) ...
            spm_select(inf, 'any', msg, '', pwd, '.*(nii|img|mat)$');
        fparts = @spm_fileparts; % (handles 'pth/fnm.ext,number' correctly)
    otherwise
        warning('make_majority_mask:spmversion', ...
            'Unrecognised version of SPM - script might not work...')
        spm5 = 1; % assume newer, but compatible...
        select = @(msg) ...
            spm_select(inf, 'any', msg, '', pwd, '.*(nii|img|mat)$');
        fparts = @spm_fileparts;
end

global defaults
if isempty(defaults)
    spm_defaults
end
