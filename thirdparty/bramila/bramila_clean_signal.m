function [cfg,group_cfg] = bramila_clean_signal(varargin)
% bramila_clean_signal.m, pipeline for cleaning fMRI data with the following operations:
% (1) masking
% (2) detrending
% (3) motion regression
% (4) tissue regression (with ISC functionality)
% finally computes subject-wise DVARS and FD
if nargin==1 % Single subject way
    cfg = varargin{1};
    group_cfg = 0;
    cfg = bramila_checkparameters(cfg);
    
    fprintf('\nFileID ''%s''\n',cfg.fileID);
    fprintf('EPI path: %s\n',cfg.infile);
    fprintf('Save folder: %s\n',cfg.outpath);
    
    % create EPI mask based on signal quality
    [cfg.vol,cfg.mask] = bramila_makemask(cfg);
    
    % nullify all voxels outside the EPI mask
    cfg.vol = bramila_maskdata(cfg);
    
    % remove trends from the data
    cfg.vol = bramila_detrend(cfg);

    % regress out motion inside EPI mask
    cfg.infile=[];
    cfg = bramila_motionregress(cfg);
    
    % regress out WM and CSF signals inside EPI mask
    cfg = bramila_tissueregress(cfg);
    
else % Group way
    cfg = varargin{1};
    group_cfg = varargin{2};
    
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
        cfg{subj}.infile=[];
        cfg{subj} = bramila_motionregress(cfg{subj});
        
        if group_cfg.save_memory==1
            % we still need to save the EPI
            if cfg{subj}.write==0
                cfg{subj}.infile = bramila_savevolume(cfg{subj},cfg{subj}.vol,'masked, detrended and motion regressed EPI','mask_detrend_motreg.nii');
            end
            cfg{subj}.vol=[];
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
        
    end
    
    if length(cfg)>1
        fprintf('\nGroup EPI mask %i voxels and analysis mask %i voxels\n',nnz(group_cfg.mask),nnz(group_cfg.analysis_mask));
        fprintf('Saving group masks\n');        
        bramila_savevolume(group_cfg,group_cfg.mask,'group EPI mask','group_mask.nii');
        bramila_savevolume(group_cfg,group_cfg.analysis_mask,'group analysis mask','group_analysis_mask.nii');
    end
end
end
