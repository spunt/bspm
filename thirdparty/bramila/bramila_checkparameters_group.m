function [cfg,group_cfg] = bramila_checkparameters_group(cfg,group_cfg)
        
    fprintf('Checking group parameters...')

    do_spatial_ISC = 0;
    if ~isfield(group_cfg,'do_spatial_ISC')
        group_cfg.do_spatial_ISC = do_spatial_ISC;
    end 
      
    if group_cfg.do_spatial_ISC>0 && isfield(group_cfg,'ISC_mask')
       if exist(group_cfg.ISC_mask,'file')==0
           error('Given ISC mask file not found!');
       end
    end
       
    a = fieldnames(group_cfg);    
    for j=1:length(cfg)
        for i=1:length(a)        
           if ~isfield(cfg{j},a{i})
                cfg{j}.(a{i})=group_cfg.(a{i});
           end
        end
        cfg{j} = bramila_checkparameters(cfg{j});
    end   
    
    group_cfg.separator = cfg{1}.separator;
    
    if ~isfield(group_cfg,'voxelsize')
        group_cfg.voxelsize=cfg{1}.voxelsize;
    end
    
    if ~isfield(group_cfg,'write')
       group_cfg.write = cfg{1}.write;
    end    
            
    if ~isfield(group_cfg,'save_path')
        group_cfg.outpath = cfg{1}.outpath;
    else
       if exist(group_cfg.save_path,'dir')~=7
           mkdir(group_cfg.save_path);
       end
       group_cfg.outpath = group_cfg.save_path;
    end
    
    group_cfg.fileID = 'bramila_groupdata';
        
    save_memory=1;
    if ~isfield(group_cfg,'save_memory')
        group_cfg.save_memory = save_memory;
    end    
    
    if is_output_overlap(cfg)==1
        if group_cfg.save_memory==1
           error('!! Overlap detected between save paths and FileIDs, either change filenames/folders/fileIDs or set ''group_cfg.save_memory=0'' !!');
        else
           warning('!! Overlap detected between save paths and FileIDs, some or all output files will be overwritten !!');
        end
    end    
    
    ISC_toolbox_path = [group_cfg.bramilapath,filesep,'external',filesep,'ISC_toolbox',filesep];
    if  ~isfield(group_cfg,'ISC_toolbox_path')
        group_cfg.ISC_toolbox_path=ISC_toolbox_path;
    end    
    
    do_temporal_tissue_ISC=0;
    if ~isfield(group_cfg,'do_temporal_tissue_ISC')
        group_cfg.do_temporal_tissue_ISC = do_temporal_tissue_ISC;
    end            
    
    temporal_tissue_ISC_method='pls';
    if ~isfield(group_cfg,'temporal_tissue_ISC_method')
        group_cfg.temporal_tissue_ISC_method = temporal_tissue_ISC_method;
    end    
    
    if ~ismember(group_cfg.temporal_tissue_ISC_method,{'full','pls','pca'})
        error('Temporal tissue ISC regression type must be {full,pls,pca}');
    end
    
    frames = zeros(1,length(cfg));
    siz = zeros(length(cfg),3);
    TRs = zeros(1,length(cfg));
    for j=1:length(cfg)
        frames(j)=get_nii_frame(cfg{j}.infile);
        nii = load_nii(cfg{j}.infile,1);        
        siz(j,1:3)=size(nii.img);
        TRs(j)=cfg{j}.TR;
    end      
    
    if length(unique(frames))>1
        for i=1:length(frames)
            fprintf('Subject %i: %i volumes\n',i,frames(i))
        end
        error('The number of volumes is not same for all subjects (use rest pipeline)!')
    end    
    if length(unique(siz(:,1)))>1 || length(unique(siz(:,2)))>1 || length(unique(siz(:,3)))>1
        error('Volume dimensions are not same for all subjects (use rest pipeline)!')
    end
    if length(unique(TRs))>1
        warning('TR values are not same for all subjects, skipping all ISC functionalities!')
        group_cfg.do_spatial_ISC=0;
        group_cfg.do_temporal_tissue_ISC=0;
    end     
           
    fprintf(' OK\n');
    
end

function res = is_output_overlap(cfg)

res=0;
for i=1:length(cfg)
    for j=(i+1):length(cfg)
        
        fileid1 = cfg{i}.fileID;
        path1 = cfg{i}.outpath;
        fileid2 = cfg{j}.fileID;
        path2 = cfg{j}.outpath;
        
        if strcmp(fileid1,fileid2) && strcmp(path1,path2)
            res = 1;
            return;
        end
    end
end

end