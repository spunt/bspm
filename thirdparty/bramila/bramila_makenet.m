function [cfg,adj,roits] = bramila_makenet(cfg,group_cfg)

fprintf('\nComputing network data\n');

if iscell(cfg)
    for subj=1:length(cfg)
        
        fprintf('Subject %i (FileID ''%s'')\n',subj,cfg{subj}.fileID);
        load(cfg{subj}.network_nodes,'rois');
        cfg{subj}.rois = rois;
        
        fprintf('..creating rois\n');
        temp_cfg.infile = cfg{subj}.infile;
        temp_cfg.vol = cfg{subj}.vol;
        temp_cfg.mask = group_cfg.analysis_mask; % limit network inside group analysis mask
        temp_cfg.rois = rois;
        [roits{subj},perc] = bramila_roiextract(temp_cfg);
        
        % Making the network
        temp_cfg.roits=roits{subj};
        temp_cfg.conn_type='pearson';
        cfg{subj}.conn_type=temp_cfg.conn_type;
        
        % add cfg.toi
        fprintf('..creating network\n');
        adj{subj}=bramila_funcconn(temp_cfg);
        adj{subj}(find(isnan(adj{subj})))=0;
        
        fprintf('..saving network\n');
        bramila_savematrix(cfg{subj},adj{subj});
        
        if cfg{subj}.write==0 && group_cfg.save_memory==1
            if ~strcmp(cfg{subj}.infile_orig,cfg{subj}.infile)
                delete(cfg{subj}.infile);
            else
                warning('!! For some reason, temporary EPI file is the original EPI file (did you skip some steps?), cannot delete !!')
            end
        end
        
    end
else
    
    load(cfg.network_nodes,'rois');
    cfg.rois = rois;
    
    fprintf('..creating rois\n');
    temp_cfg.infile = cfg.infile;
    temp_cfg.vol = cfg.vol;
    temp_cfg.mask = cfg.analysis_mask;
    temp_cfg.rois = rois;
    [roits,perc] = bramila_roiextract(temp_cfg);
    
    % Making the network
    temp_cfg.roits=roits;
    temp_cfg.conn_type='pearson';
    cfg.conn_type=temp_cfg.conn_type;
    
    % add cfg.toi
    fprintf('..creating network\n');
    adj=bramila_funcconn(temp_cfg);
    adj(find(isnan(adj)))=0;
    if(cfg.write==1)
    fprintf('..saving network\n');
    bramila_savematrix(cfg,adj);
    end
    
end