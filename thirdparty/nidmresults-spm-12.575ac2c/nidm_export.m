function nidm_export(data_path, out_path)
    cwd = pwd;
    cd(data_path)
    test_name = spm_file(data_path, 'filename');
    
    if strcmp(test_name, 'spm_group_wls')
        study_dir = fullfile(pwd, 'mfx');
    else
        study_dir = pwd;
    end    
    % Remove previous nidm exports    
    files = dir(study_dir);
    subdirs = files([files.isdir]);
    for i = 1:numel(subdirs)
        dname = subdirs(i).name;
        if strncmpi(dname,'nidm',4)
            disp(['Removing ' dname])
            rmdir(dname,'s')
        end
        
        nidm_zips = cellstr(strvcat(spm_select('FPList', study_dir, '\.nidm\.zip$')));
        for j = 1:numel(nidm_zips)
            if ~isempty(nidm_zips{j})
                disp(['Deleting ' nidm_zips{j}])
                delete(nidm_zips{j})
            end
        end
    end
    
    if strcmp(test_name, 'spm_full_example001')
        % For SPM full example 001 we use already exported peaks 
        % and clusters list to get exactly the same graph
        load(fullfile(data_path, 'nidm_example001.mat'));
        [SPM, xSPM] = set_study_path(SPM, xSPM, pwd);
        % FIXME: should be extracted from json (when reader fixed)
        subjects.subject = 1;
        modality = 'FMRI';
        space = 'ixi';
        spm_results_nidm(SPM,xSPM,TabDat,subjects,modality,space);
    else
        run(fullfile(pwd, 'batch.m'))
        result_batch = matlabbatch(end);
        
        result_batch{1}.spm.stats.results.spmmat = {fullfile(study_dir, 'SPM.mat')};
        
        % FIXME: this will have to be read from json when json reader is
        % fixed
        if isempty(findstr(test_name, 'group'))
            group_analysis = false;
        else
            group_analysis = true;
        end
        
        if ~group_analysis
            % single-subject analysis (fMRI / subject-space)         
            result_batch{1}.spm.stats.results.export.nidm.subjects.subject = 1;
            result_batch{1}.spm.stats.results.export.nidm.modality = 2;
            result_batch{1}.spm.stats.results.export.nidm.refspace = 1;
        else
            % group-analysis (fMRI / segment space)
            result_batch{1}.spm.stats.results.export.nidm.subjects.group.label = 'Control';
            % FIXME: this will have to be read from json when json reader is
            % fixed
            result_batch{1}.spm.stats.results.export.nidm.subjects.group.numsubjects = 14;
            result_batch{1}.spm.stats.results.export.nidm.modality = 2;
            result_batch{1}.spm.stats.results.export.nidm.refspace = 2;
        end
        
        try
            spm_jobman('run', result_batch)
        catch ME
            switch ME.identifier
                case 'matlabbatch:run:jobfailederr'
                    % Voxel-wise FDR requires topoFDR to be disabled
                    global defaults;
                    defaults.stats.topoFDR = 0;
                    spm_jobman('run', result_batch)
                    defaults.stats.topoFDR = 1;
                otherwise
                    rethrow(ME)
            end
        end
    end
    
    unzip('spm_0001.nidm.zip', 'nidm')
    
    if ~isempty(out_path)
        test_name = spm_file(data_path, 'basename');
        
        target_dir = fullfile(out_path, ['ex_' test_name]);
        if isdir(target_dir)
            disp(['Removing ' target_dir])
            rmdir(target_dir,'s')
        end
        movefile('nidm', target_dir)

%         error on mac to be fixed
%         spm_jsonwrite(fullfile(target_dir, 'config.json'), json_cfg)
        json_file = fullfile(data_path, 'config.json');
%         aa = spm_jsonread(json_file)
%         aa=1
        copyfile(json_file, fullfile(target_dir, 'config.json'));
             
        fname = json_file;
        fid = fopen(fname);
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);

        expression = '\[".*\.ttl"\]';
        gts = regexp(str,expression,'match'); 
        gts = strrep(strrep(strrep(gts{1}, '[', ''), ']', ''), '"', '');
        gts = strsplit(gts, ', ');
        
        for i = 1:numel(gts)
            gt = gts{i};
    %         disp(gt)
    %         gt = json_cfg.ground_truth;
            version = regexp(str,'"versions": \[".*"\]','match');
            version = strrep(strrep(strrep(version{1},'"versions": ["', ''), '"', ''), ']', '');
            gt_file = fullfile(data_path, '..', '_ground_truth', version, gt);

    %         target_gt_dir = fullfile(out_path, 'ground_truth', version, spm_file(gt,'path'));
    %         disp(gt)
            % FIXME: version should be extracted from json        
    %         gt_file = fullfile(path_to_script_folder, '..', 'ground_truth', '1.2.0', gt);

            target_gt_dir = fullfile(out_path, 'ground_truth', version, spm_file(gt,'path'));
            if isdir(target_gt_dir)
                disp(['Removing ' target_gt_dir])
                rmdir(target_gt_dir,'s')
            end
            mkdir(target_gt_dir)
            copyfile(gt_file, target_gt_dir);
        end
    end
    
    cd(cwd);
end

function [SPM, xSPM] = set_study_path(SPM, xSPM, new_dir)
    SPM.swd = new_dir;
    xSPM.swd = new_dir;
end