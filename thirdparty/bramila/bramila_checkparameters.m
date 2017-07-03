function cfg = bramila_checkparameters(cfg)
    
    % defaults:
    voxelsize=[2,2,2];
    mask = 'MNI152_ENLARGED';
    tissue_derivatives = 0;
    min_tissue_var_expl = 50;
    max_tissue_pca_count = 12;
    remove_global = 0;
    save_path = [];
    motion_reg_type = 'friston';
    mot_derivatives = 1;
    white_mask_th = 0.90;
    csf_mask_th = 0.90;     
    detrend_type='linear-nodemean';
    filter_type = 'butter';
    filter_limits = [0,0.01,0.08,0.09];
    write = 0;      
    network_nodes = [];
    do_spatial_ISC=0;
    do_spatial_smooth=0;
    smooth_FWHM = 6;
    smooth_method = 'SPM';
    
    print_cfg = 0;    
    
%    fprintf('Checking parameters...');

    % set defaults
    
    cfg.separator = filesep();  
    
    addpath([cfg.bramilapath,cfg.separator,'external',cfg.separator,'niftitools']);
    addpath([cfg.bramilapath,cfg.separator,'external',cfg.separator,'bin']);
        
    if ~isfield(cfg,'voxelsize')
        cfg.voxelsize=voxelsize;
    end

    if ~isfield(cfg,'FSLDIR')
        cfg.FSLDIR='';
    end	
    
    if ~isfield(cfg,'write')
        cfg.write=write;
    end    
        
    if ~isfield(cfg,'detrend_type')
        cfg.detrend_type = detrend_type;
    end   
    
    if ~isfield(cfg,'filter_type')
        cfg.filter_type = filter_type;
    end   

    if ~isfield(cfg,'mask')
        cfg.mask  = mask;
    end     
    
    if ~isfield(cfg,'do_spatial_ISC')
        cfg.do_spatial_ISC  = do_spatial_ISC;
    end         
       
    if ~isfield(cfg,'white_mask_th')
        cfg.white_mask_th  = white_mask_th;
    end   
    
    if ~isfield(cfg,'csf_mask_th')
        cfg.csf_mask_th = csf_mask_th;
    end
    
    if ~isfield(cfg,'white_mask')
        cfg.white_mask = [cfg.bramilapath,cfg.separator,'external',cfg.separator,'white.nii'];
    end   
    
    if ~isfield(cfg,'csf_mask')
        cfg.csf_mask = [cfg.bramilapath,cfg.separator,'external',cfg.separator,'csf.nii'];
    end       
    
    if ~isfield(cfg,'filter_limits')
        cfg.filter_limits = filter_limits;
    end      
    
    if ~isfield(cfg,'motion_reg_type')
        cfg.motion_reg_type=motion_reg_type;
    end    

    if ~isfield(cfg,'mot_derivatives')
        if strcmp(cfg.motion_reg_type,'standard')
            cfg.mot_derivatives=mot_derivatives;
        else
            cfg.mot_derivatives=0;
        end
    end
    
    if ~isfield(cfg,'tissue_derivatives')
        cfg.tissue_derivatives=tissue_derivatives;
    end

    if ~isfield(cfg,'min_tissue_var_expl')
        cfg.min_tissue_var_expl=min_tissue_var_expl;
    end

    if ~isfield(cfg,'max_tissue_pca_count')
        cfg.max_tissue_pca_count = max_tissue_pca_count;
    end
    
    if ~isfield(cfg,'remove_global') 
        cfg.remove_global=remove_global;
    end    
    
    if ~isfield(cfg,'network_nodes')
       cfg.network_nodes = [cfg.bramilapath,cfg.separator,'external',cfg.separator,'rois6mm_final_plus_IF2.mat'];
    end
    
    % check parameters
    
    if ~isfield(cfg,'TR') || (isfield(cfg,'TR') && cfg.TR<=0 && cfg.TR>50)
        error('TR value missing or invalid!')
    end    
    
    if strcmp(cfg.infile(end-3:end),'.nii')
        if exist(cfg.infile,'file')==0
            error('Input NIFTI file not found\n%s!',cfg.infile)
        end         
        if isfield(cfg,'fileID') && ~isempty(cfg.fileID)
            a = cfg.fileID;
            cfg = createID(cfg);
            cfg.fileID = a;
        else
           cfg = createID(cfg); 
        end
        cfg.infile_orig = cfg.infile; % only used for fail-safe reasons
    else
        error('Input EPI file must be 4D nifti!\n%s!',cfg.infile)
    end
    
    if exist(cfg.motionparam,'file')==0
        error('Motion parameters file not found!\n%s',cfg.motionparam)
    end            
    
    if (~isfield(cfg,'save_path') ) || (isfield(cfg,'save_path')  && isempty(cfg.save_path) )
       cfg.outpath=cfg.inpath;
    else
        if exist(cfg.save_path,'dir')~=7
           mkdir(cfg.save_path);
        else
           %warning('Using custom output path. EPIs from different subjects with the same filename (e.g. EPI.nii) will be overwritten!');
        end
        cfg.outpath = cfg.save_path;
    end
    
    if isfield(cfg,'white_mask') && isfield(cfg,'csf_mask')
        nii=load_nii(cfg.infile,1);
        s = size(nii.img);
        
        nii=load_nii(cfg.white_mask);
        if ~all(size(nii.img)==s(1:3))
            error('Given WM mask size does not match EPI size!\n%s',cfg.white_mask)
        end
        vals=nii.img(:);
        vals(isnan(vals))=[];
        if min(vals)<0 || max(vals)>1
            warning('Values are not between [0,1] for a given WM mask ([%f,%f])!',min(vals),max(vals));
        end
        
        nii=load_nii(cfg.csf_mask);
        if ~all(size(nii.img)==s(1:3))
            error('Given CSF mask size does not match EPI size!\n%s',cfg.csf_mask)
        end        
        vals=nii.img(:);
        vals(isnan(vals))=[];
        if min(vals)<0 || max(vals)>1
            warning('Values are not between [0,1] for a given CSF mask ([%f,%f])!',min(vals),max(vals));
        end        
    end
    
    if ~ismember(cfg.motion_reg_type,{'standard','friston','volterra'})
        error('Motion regression type must be {standard,friston,volterra}');
    end
    
    if ~all(cfg.filter_limits==sort(cfg.filter_limits))
        error('filter limits must be in increasing order')
    end
    
    if strcmp(cfg.mask,'MNI152')
        nii = load_nii([cfg.bramilapath,cfg.separator,'external',cfg.separator,'MNI152_T1_2mm_brain_mask.nii']);
        cfg.mask = nii.img;
    elseif strcmp(cfg.mask,'MNI152_ENLARGED')
        nii = load_nii([cfg.bramilapath,cfg.separator,'external',cfg.separator,'MNI152_T1_2mm_brain_mask_ENLARGED.nii']);
        cfg.mask = nii.img;        
    elseif ischar(cfg.mask)
        nii=load_nii(cfg.mask);
        cfg.mask = nii.img;
    elseif (isnumeric(cfg.mask) || islogical(cfg.mask)) && length(size(cfg.mask))==3
        % assume that proper mask-matrix is entered
    elseif isempty(cfg.mask)
        % no mask
    else
        error('Input mask has unknown format!')
    end

    if ~exist(cfg.network_nodes,'file')
           error('Network node ROI definition file not found!');
    end
    
     if ~isfield(cfg,'do_spatial_smooth')
         cfg.do_spatial_smooth=do_spatial_smooth;
     end
     
     if ~isfield(cfg,'smooth_FWHM')
         cfg.smooth_FWHM=smooth_FWHM;
     else
         if isscalar(cfg.smooth_FWHM) && cfg.smooth_FWHM>=0
             % ok
         else
            error('FWHM must be a scalar (isotropic kernel)');
         end
     end
         
     if ~isfield(cfg,'smooth_method')
         cfg.smooth_method='SPM';
     end        
     
    if ~ismember(cfg.smooth_method,{'SPM','FSL','AFNI','FSLgauss'}) 
        error('Unknown smoothing function (must be SPM, FSL, FSLgauss or AFNI)')
    end
        
    if print_cfg==1
       s=fieldnames(cfg);
       fprintf('---parameters---\n');
       for i=1:length(s)
           fprintf('%s = %s\n',s{i},num2str(cfg.(s{i})));
       end
       fprintf('----------------\n');
    end
        

end

function cfg = createID(cfg)
    
    ind = strfind(cfg.infile,cfg.separator);
    if isempty(ind)   
       cfg.fileID = cfg.infile(1:(end-4));
       cfg.inpath='';
    else
       cfg.fileID = cfg.infile((ind(end)+1):(end-4));
       cfg.inpath = cfg.infile(1:(ind(end)-1));
    end
    
end
