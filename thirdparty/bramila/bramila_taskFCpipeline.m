function [adj,roits,cfg,group_cfg] = bramila_taskFCpipeline( cfg , group_cfg )
% Pipeline for fmri task-data. Otherwise similar to bramila_restFCpipeline,
% but takes the whole group and does ISC-based tissue signal regression.
% 

[cfg,group_cfg] = bramila_checkparameters_group(cfg,group_cfg);

for subj = 1:length(cfg)    
    
    fprintf('\nSubject %i (FileID ''%s'')\n',subj,cfg{subj}.fileID);
    fprintf('EPI path: %s\n',cfg{subj}.infile);
    fprintf('Save folder: %s\n',cfg{subj}.outpath);
    
    % create EPI mask based on signal quality
    [cfg{subj}.vol,cfg{subj}.mask] = bramila_makemask(cfg{subj});
    
    % nullify all voxels outside the EPI mask
    cfg{subj}.vol = bramila_maskdata(cfg{subj});
    
    % remove trends from the data
    cfg{subj}.vol = bramila_detrend(cfg{subj});
    
    % regress out motion inside EPI mask
    cfg{subj}.vol = bramila_motionregress(cfg{subj});
    
    if cfg{subj}.write==0
        % we still need to save the EPI
        if group_cfg.save_memory==1            
            cfg{subj}.infile = bramila_savevolume(cfg{subj},cfg{subj}.vol,'masked, detrended and motion regressed EPI','mask_detrend_motreg.nii');
            cfg{subj}.vol=[];
        elseif group_cfg.do_spatial_ISC==1 && ( (isfield(group_cfg,'ISC_mask') && isempty(group_cfg.ISC_mask)) || (~isfield(group_cfg,'ISC_mask') ))           
            cfg{subj}.infile = bramila_savevolume(cfg{subj},cfg{subj}.vol,'masked, detrended and motion regressed EPI','mask_detrend_motreg.nii');
        end
    end
    
end

fprintf('\n');

% regress out WM and CSF signals inside EPI mask, takes into account
% inter-subject correlations
[cfg,group_cfg] = bramila_tissueregress_task(cfg,group_cfg);

for subj = 1:length(cfg)
    
    if subj==1
        group_cfg.mask = cfg{subj}.mask;
        group_cfg.analysis_mask = cfg{subj}.analysis_mask;
    else
        group_cfg.mask = group_cfg.mask.*cfg{subj}.mask;
        group_cfg.analysis_mask = group_cfg.analysis_mask.*cfg{subj}.analysis_mask;        
    end
    
    if group_cfg.save_memory==1
        nii=load_nii(cfg{subj}.infile);
        cfg{subj}.vol=nii.img;
    end
    
    fprintf('\nSubject %i (FileID ''%s'')\n',subj,cfg{subj}.fileID);
    % compute quality control measures        
    cfg{subj}.dvars = bramila_dvars(cfg{subj});    
    [cfg{subj}.fDisplacement,cfg{subj}.fDisplacement_rms]=bramila_framewiseDisplacement(cfg{subj});
    
    % temporally filter the data inside EPI mask
    cfg{subj} = bramila_filter(cfg{subj});
    
    % load node definition file
    load(cfg{subj}.network_nodes,'rois');
    cfg{subj}.rois = rois;
    
    % create network
    [adj{subj},roits{subj}] = bramila_makenet(cfg{subj});
    
    if group_cfg.save_memory==1 && cfg{subj}.write==0
       % delete temporary EPI file
       cfg{subj}.vol = [];
       if ~strcmp(cfg{subj}.infile_orig,cfg{subj}.infile)
            delete(cfg{subj}.infile);
       else
            warning('!! For some reason, temporary EPI file is the original EPI file (did you skip some steps?), cannot delete !!')
       end
    end
    
end

end
