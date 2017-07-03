function [cfg,group_cfg] = bramila_tissueregress_task(cfg,group_cfg)
% INPUT
%  cfg = subject-wise configuration file (see bramila_checkparameters.m)
%  group_cfg = group-wise configuration file (see bramila_checkparameters_group.m)
% OUTPUT
%  modified cfg and group_cfg files

if isfield(group_cfg,'do_spatial_ISC') && group_cfg.do_spatial_ISC==1
        
    if ~isfield(group_cfg,'ISC_mask')
        
        if group_cfg.write==0 && group_cfg.save_memory==0
            for i=1:length(cfg)
                cfg{subj}.infile = bramila_savevolume(cfg{subj},cfg{subj}.vol,'masked, detrended and motion regressed EPI','mask_detrend_motreg.nii');
            end
        end        
        % compute new standard ISC mask with threshold p<0.05 FDR,
        % using folder of the first subject!
        fprintf('-----------Creating ISC mask (ISC toolbox)-------------\n');
        group_cfg.ISC_path = [cfg{1}.outpath,cfg{1}.separator,'bramila',cfg{1}.separator,'spatial_ISC',cfg{1}.separator];
        group_cfg.ISC_mask = bramila_compute_ISC_mask(cfg,group_cfg);
        fprintf('\n-----------Bramila pipeline continues-------------\n');
    end
    
    for i=1:length(cfg)
        cfg{i}.ISC_mask = group_cfg.ISC_mask;
    end
    
end

for i=1:length(cfg)
    if isfield(cfg{i},'white_mask') && isfield(cfg{i},'csf_mask') && (cfg{i}.max_tissue_pca_count>0 || cfg{i}.remove_global>0)
        if ~isfield(cfg{i},'vol') || (isfield(cfg{i},'vol') && isempty(cfg{i}.vol))
            nii=load_nii(cfg{i}.infile);
            cfg{i}.vol=double(nii.img);
        end
        fprintf('Creating tissue regressors for subject %i (FileID ''%s'')\n',i,cfg{i}.fileID);
        [reg_cfg{i},cfg{i}] = bramila_maketissueregress(cfg{i});
        if group_cfg.save_memory==1
            cfg{i}.vol=[];
            reg_cfg{i}.vol=[];
        end
    else
        reg_cfg{i}=[];
    end
end

if group_cfg.do_temporal_tissue_ISC>0
    
    if group_cfg.do_temporal_tissue_ISC==1
       pval = bramila_mvcorrtest(reg_cfg,5000);     
    else
       fprintf('\nTissue ISC test skipped\n');
       pval = zeros(1,3);
    end
    method = group_cfg.temporal_tissue_ISC_method;
    [reg_cfg,var_expl,cv_data] = bramila_make_group_regressors(reg_cfg,method,pval<0.05);
    group_cfg.tissue_ISC_var = var_expl;
    if ~isempty(cv_data)           
        group_cfg.tissue_ISC_cvdata=cv_data;
    end

end

fprintf('\n');

for i=1:length(cfg)
    
    if ~isempty(reg_cfg{i})
        
        if isfield(cfg{i},'mask')
            reg_cfg{i}.mask = cfg{i}.mask;
        end
        if isempty(reg_cfg{i}.vol)
            nii=load_nii(cfg{i}.infile);
            reg_cfg{i}.vol=double(nii.img);
        end
        fprintf('Running tissue regression for subject %i (FileID ''%s'')\n',i,cfg{i}.fileID);
        [clean_vol,r2] = bramila_regress(reg_cfg{i});        
        reg_cfg{i}.vol=[];
        
        if group_cfg.save_memory==0
            cfg{i}.vol = clean_vol;
        end        
        
        if cfg{i}.write==0 && group_cfg.save_memory==1
            if ~strcmp(cfg{i}.infile_orig,cfg{i}.infile)
                delete(cfg{i}.infile);
            else
                warning('!! For some reason, temporary EPI file is the original EPI file (did you skip some steps?), cannot delete !!')
            end
        end
        
        if cfg{i}.write==1 || group_cfg.save_memory==1
            cfg{i}.infile = bramila_savevolume(cfg{i},clean_vol,'EPI volume after masking, detrending, motion and tissue regression','mask_detrend_fullreg.nii');
        end
        
        bramila_savevolume(cfg{i},r2,'tissue regression R2','tissue_R2.nii');
        bramila_savedesign(cfg{i},reg_cfg{i}.reg,reg_cfg{i}.labels,'tissue_reg_design.mat');
        bramila_savevolume(cfg{i},cfg{i}.analysis_mask,'Analysis mask','analysis_mask.nii');
        
    else
        
        fprintf('Skipping tissue regression for subject %i (fileID %s)\n',i,cfg{i}.fileID);
        cfg{i}.analysis_mask=cfg{i}.mask;
        bramila_savevolume(cfg{i},cfg{i}.analysis_mask,'Analysis mask','analysis_mask.nii');
        
    end
    
end

end


