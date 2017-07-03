function bramila_savedesign(cfg,mat,lab,filename)
if(~isfield(cfg,'separator'))
    cfg.separator='/';  % defaults to *nix separator
end

if(~isfield(cfg,'inpath'))
	warning('bramila_savedesign>> cfg.inpath is empty');
	cfg.inpath='';
end

if(~isfield(cfg,'fileID'))
    warning('bramila_savedesign >> fileID is empty');
    cfg.fileID='';
end


if ~exist([cfg.outpath,cfg.separator,'bramila'],'dir')
    mkdir([cfg.outpath,cfg.separator,'bramila']);
end
design_matrix.mat = mat;
design_matrix.labels = lab;
save([cfg.outpath,cfg.separator,'bramila',cfg.separator,cfg.fileID,'_',filename],'design_matrix');

end
