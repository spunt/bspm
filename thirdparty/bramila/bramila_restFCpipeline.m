function [adj,roits,cfg] = bramila_restFCpipeline( cfg )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%----Run bramila pipeline code (should not be modified by the user) -----------------------

cfg = bramila_checkparameters(cfg);

fprintf('\nFileID: ''%s''\n',cfg.fileID);
fprintf('EPI path: %s\n',cfg.infile);
fprintf('Save folder: %s\n',cfg.outpath);

% create EPI mask based on signal quality
[cfg.vol,cfg.mask] = bramila_makemask(cfg);

% nullify all voxels outside the EPI mask
cfg.vol = bramila_maskdata(cfg);

% remove trends from the data
cfg.vol = bramila_detrend(cfg);

% regress out motion inside EPI mask
cfg = bramila_motionregress(cfg);

% regress out WM and CSF signals inside EPI mask
cfg = bramila_tissueregress(cfg);

% compute quality control measures
cfg.dvars = bramila_dvars(cfg);
[cfg.fDisplacement,cfg.fDisplacement_rms]=bramila_framewiseDisplacement(cfg);

% temporally filter the data inside analysis mask
cfg = bramila_filter(cfg);

% load node definition file
load(cfg.network_nodes,'rois');
cfg.rois = rois;

% create network
[cfg,adj,roits] = bramila_makenet(cfg);

end

