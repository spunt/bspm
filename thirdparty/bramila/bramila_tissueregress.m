function cfg = bramila_tissueregress(cfg)
% INPUT
%   cfg.infile = location where the subject NII file (4D)
%   cfg.vol = input can be a 4D volume, if this exists, then infile is not loaded but used as ID (optional)
%   cfg.write = 0 (default 0, set to 1 if you want to store the volume)
%   Janne to invent other parameters related to masks (e.g. individual masks or standard masks, % of prob for each mask
% OUTPUT
%   vol = a 4D volume tissue regressed

if ~isfield(cfg,'vol') || (isfield(cfg,'vol') && isempty(cfg.vol))
    nii=load_nii(cfg.infile);
    cfg.vol=double(nii.img);
end

if isfield(cfg,'white_mask') && isfield(cfg,'csf_mask') && (cfg.max_tissue_pca_count>0 || cfg.remove_global>0)
    
    fprintf('Creating tissue regressors\n');
    
    [reg_cfg,cfg] = bramila_maketissueregress(cfg);    
    
    if isfield(cfg,'mask')
        reg_cfg.mask = cfg.mask;
    end
    [clean_vol,r2] = bramila_regress(reg_cfg);
    cfg.vol = clean_vol;
    
    if cfg.write==1
        cfg.outfile=bramila_savevolume(cfg,clean_vol,'EPI volume after masking, detrending, motion and tissue regression','mask_detrend_fullreg.nii');
	else
		if(isfield(cfg,'outfile'))
			rmfield(cfg,'outfile');
		end
	end
    bramila_savevolume(cfg,r2,'tissue regression R2','tissue_R2.nii');
    bramila_savedesign(cfg,reg_cfg.reg,reg_cfg.labels,'tissue_reg_design.mat');
    bramila_savevolume(cfg,cfg.analysis_mask,'Analysis mask','analysis_mask.nii');
    
else
    
    fprintf('Skipping tissue regression\n');
    cfg.analysis_mask=cfg.mask;
    bramila_savevolume(cfg,cfg.analysis_mask,'Analysis mask','analysis_mask.nii');
    
end

end

