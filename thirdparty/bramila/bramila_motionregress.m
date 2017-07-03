function cfg=bramila_motionregress(cfg)
% BRAMILA_MOTIONREGRESS - Regresses out motion parameters as described in 
% Power et al. (2014) doi:10.1016/j.neuroimage.2013.08.048
% 	- Usage:
%	cfg = bramila_motionregress(cfg);
%	- Input:
%	cfg is a struct with following parameters
%		cfg.infile = string with location where the subject NII file (4D)
%   	cfg.motionparam = location where the subject motion param file (txt)
%   	cfg.vol = input can be a 4D volume, if this exists, then infile is not loaded but used as ID (optional)
%  		cfg.write = 0 (default 0, set to 1 if you want to store the volume)
%		cfg.save_path = optional, must be set if cfg.write=1. If this is set it will also save the regressors used
%		cfg.motion_reg_type = one of 'standard','friston','volterra' (default 'friston')    
%		cfg.mot_derivatives = 0-2  number of motion derivatives (default 1)
%	- Output:
%		cfg = similar to the input, plus added are the motion regressors, the regressed volume and other default values     
%		cfg.vol = regressed volume
% 		cfg.motionR2 = R2 values of voxels affected by motion regression
%    	cfg.motion_regressors = regressors used
%		cfg.motion_regressors_labels = labels for the regressors
%
% 	Last edit: EG 2014-04-24

    %% Default values, same as bramila_checkparams
	motion_reg_type = 'friston';
    mot_derivatives = 1;
	write = 0;
	
	%% input validation

	% data: 
    if ~isfield(cfg,'vol') || (isfield(cfg,'vol') && isempty(cfg.vol))
    	if(strcmp(class(cfg.infile),'char') )    
			nii=load_nii(cfg.infile);
			cfg.vol=double(nii.img);
		else
			error(['bramila_motionregress >> Variable cfg.infile is not a valid string']);
		end
    end
    
    % motion type
	if ~isfield(cfg,'motion_reg_type')
        cfg.motion_reg_type=motion_reg_type;
    end
	if ~ismember(cfg.motion_reg_type,{'standard','friston','volterra'})
        error('Motion regression type must be {standard,friston,volterra}');
    end

	% motion derivative
    if ~isfield(cfg,'mot_derivatives')
        if strcmp(cfg.motion_reg_type,'standard')
            cfg.mot_derivatives=mot_derivatives;
        else
            cfg.mot_derivatives=0; % i.e. for the friston or volterra case derivatives are there already
        end
    end

	% motion parameter
	if exist(cfg.motionparam,'file')==0
        error('Motion parameters file not found!\n%s',cfg.motionparam)
    end

	%% BEGIN 
	% first we need to create the motion regressors

    fprintf('Creating motion regressors...');
    
    reg_cfg = bramila_makemotionregress(cfg);
    
    fprintf(' done\n');
    
    if isfield(cfg,'mask')
        reg_cfg.mask = cfg.mask;
    end
    [clean_vol,r2] = bramila_regress(reg_cfg);
    cfg.vol=clean_vol;
    cfg.motionR2=r2;
	cfg.motion_regressors=reg_cfg.reg;
	cfg.motion_regressors_labels=reg_cfg.labels; 
	
    %% write to disk if it's what we want
	% check that outpath exists and is not empty, otherwise use save_path
	outpathOK=0;
	if(isfield(cfg,'outpath') && ~isempty(cfg.outpath))
		outpathOK=1;
	end

	savepathOK=0;
	if(isfield(cfg,'save_path') && ~isempty(cfg.save_path))
		savepathOK=1;
	end	

	if(outpathOK==0 && savepathOK==1)
		cfg.outpath=cfg.save_path;
		outpathOK=1;
	end
		
	if (outpathOK==1 )
	% do the saving only if we have a save_path field set
		if cfg.write==1 || nargout<1
			cfg.infile = bramila_savevolume(cfg,clean_vol,'masked, detrended and motion regressed EPI','mask_detrend_motreg.nii');
			cfg.outfile = cfg.infile;
		else
			if(isfield(cfg,'outfile'))
				rmfield(cfg,'outfile');
			end
		end
		bramila_savevolume(cfg,r2,'motion regression R2','mot_reg_R2.nii');
		bramila_savedesign(cfg,reg_cfg.reg,reg_cfg.labels,'mot_reg_design.mat');
	else
		disp('bramila_motionregress >> No save_path set, I will not save volumes or regressors');
	end
	
end
