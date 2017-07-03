function filepath = bramila_savevolume(cfg,vol,description,filename)
if(~isfield(cfg,'separator'))
	cfg.separator='/';	% defaults to *nix separator
end

if(~isfield(cfg,'fileID'))
	warning('bramila_savevolume >> fileID is empty');
	cfg.fileID='';
end

if(~isfield(cfg,'voxelsize'))
	warning('bramila_savevolume >> assuming voxelsize of 2mm');
	cfg.voxelsize=2;
end

if ~exist([cfg.outpath,cfg.separator,'bramila'],'dir')
    mkdir([cfg.outpath,cfg.separator,'bramila']);
end

filepath = [cfg.outpath,cfg.separator,'bramila',cfg.separator,cfg.fileID,'_',filename];
% fix origin
ref = load_nii(cfg.StdTemplate); % might make sense to keep image loaded throughout the whole pipeline?
ref.hdr.dime.bitpix=64;
ref.hdr.dime.datatype=64;
ref.img = double(vol);
ref.hdr.hist.descrip=description;
ref.hdr.dime.cal_max=1000;
ref.hdr.dime.cal_min=0;
siz = size(vol);
if length(siz)==4
    ref.hdr.dime.dim(1)=4;
    ref.hdr.dime.dim(5)=siz(4);
end
save_nii(ref,filepath);

end


