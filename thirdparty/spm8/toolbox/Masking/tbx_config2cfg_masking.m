function tbx_config2cfg_masking
% tbx_config_masking is for SPM5; this function takes that file and
% automatically generates a pair of SPM8 cfg and def files, using
% spm_tbx_config2cfg followed by changing the .vfiles references to .vout
% since my vfiles/vout call-back functions handle both, checking spm('Ver')

if ~( strcmp(spm('Ver'), 'SPM8b') || strcmp(spm('Ver'), 'SPM8') )
    return;
end

if ~exist('spm_tbx_config2cfg', 'file')
    addpath(fullfile(spm('Dir'), 'toolbox'));
end

owd = pwd;
cd(fileparts(which('tbx_config2cfg_masking')));

spm5config = tbx_config_masking;
cfgname = spm_tbx_config2cfg(spm5config);
fid = fopen(sprintf('%s.m', cfgname), 'r');
cfg = fscanf(fid, '%c');
fclose(fid);
cfg = regexprep(cfg, '\.vfiles', '.vout');
% addpath lost from config2cfg, restore:
addpstr = 'addpath(fullfile(spm(''Dir''), ''toolbox'', ''Masking''))';
oldstr = '% This code has been automatically generated.';
newstr = sprintf('%s\n%s', oldstr, addpstr);
cfg = regexprep(cfg, oldstr, newstr);
fid = fopen(sprintf('%s.m', cfgname), 'w');
fprintf(fid, '%s', cfg);
fclose(fid);

cd(owd)
