function nidm_export_all(path_to, out_path)
    if ~isvarname('out_path')
        out_path = '';
    end

    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_get_defaults('cmdline',true);

    files = dir(path_to);
    subdirs = files([files.isdir]);
    subdirs = setdiff(cellstr(strvcat(subdirs.name)), ...
        {'.', '..', '.git', 'ground_truth'});
    for i = 1:numel(subdirs)
        dname = subdirs{i};
        % FIXME: examples ignored so far
        if ~strcmp(dname, 'spm_explicit_mask') && ...
           ~strcmp(dname, 'spm_2_t_test') && ...
           ~strcmp(dname, 'spm_covariate')
            json_file = fullfile(path_to, dname, 'config.json');
            
            if exist(json_file, 'file') ~= 2
                warning(['No config.json for ' dname])
                continue;
            end
%             json read error to be fixed
%             cfg = spm_jsonread(json_file);
            
%             if strcmp(lower(cfg.software), 'spm')
            if strncmpi(dname,'spm_',4)
                disp(dname)
                nidm_export(fullfile(path_to, dname), out_path)
            end
%             end
        end
    end
end