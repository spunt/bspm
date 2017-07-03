function cfg = bramila_spatiotemporal_filtering(varargin)
% bramila_spatiotemporal_filtering.m, does temporal and spatial filtering
% for cleaned EPI data.

if nargin==1 % Single subject way
    cfg = varargin{1};
    
    fprintf('\nRunning spatiotemporal processing (FileID ''%s'')\n',cfg.fileID);
    
    cfg = bramila_filter(cfg);
    
    if cfg.do_spatial_smooth==1
        cfg = bramila_smooth(cfg);
    end
    
else % Group way
    cfg = varargin{1};
    group_cfg = varargin{2};
    
    for subj=1:length(cfg)
        
        fprintf('\nRunning spatiotemporal processing for subject %i (FileID ''%s'')\n',subj,cfg{subj}.fileID);
        
        if isfield(cfg{subj},'vol') && ~isempty(cfg{subj}.vol)
            % OK
        elseif isfield(cfg{subj},'infile') && ~isempty(cfg{subj}.infile)
            nii=load_nii(cfg{subj}.infile);
            cfg{subj}.vol=nii.img;
        else
            error('No volume or infile found!')
        end
        
        if cfg{subj}.write==0 && group_cfg.save_memory==1
            if ~strcmp(cfg{subj}.infile_orig,cfg{subj}.infile)
                delete(cfg{subj}.infile);
            else
                warning('!! For some reason, temporary EPI file is the original EPI file (did you skip some steps?), cannot delete !!')
            end
        end
           
        cfg{subj}.infile=[];
        
        cfg{subj} = bramila_filter(cfg{subj});
        
        if cfg{subj}.do_spatial_smooth==1
            cfg{subj} = bramila_smooth(cfg{subj},group_cfg);
        end
                
        if group_cfg.save_memory==1
            if cfg{subj}.write==0
                % we still need to save the EPI
                cfg{subj}.infile = bramila_savevolume(cfg{subj},cfg{subj}.vol,'masked, detrended and motion regressed EPI','mask_detrend_fullreg_filtered_smoothed.nii');
            end
            cfg{subj}.vol=[];
        end
        
    end
end