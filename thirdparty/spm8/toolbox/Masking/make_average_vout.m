function vout = make_average_vout(job)
% handle virtual output for dependencies in subsequent batch jobs

% vfiles for SPM5, or vout for SPM8 matlabbatch system

if strcmp(spm('Ver'), 'SPM5')
    [pth nam ext num] = spm_fileparts(job.outname);
    if isempty(num) % add this, because spm_check_registration looks for it
        num = ',1';
    end
    if isempty(pth)
        if isempty(job.outdir{1})
            pth = pwd;
        else
            pth = job.outdir{1};
        end
    end
    vout{1} = fullfile(pth, [nam ext num]);
else % assume SPM8 or newer (!)
    vout            = cfg_dep;
    vout.sname      = 'Average Image';
    vout.src_output = substruct('.', 'files');
    vout.tgt_spec   = cfg_findspec({{'filter','image', 'strtype','e'}});
end
