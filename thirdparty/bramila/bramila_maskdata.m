function vol = bramila_maskdata(cfg)
% INPUT
%   cfg.infile = location where the subject NII file (4D)
%   cfg.vol = input can be a 4D volume, if this exists, then infile is not loaded but used as ID (optional)
%   cfg.write = 0 (default 0, set to 1 if you want to store the volume)
%   cfg.mask = optional, mask to use (it can be a 3D vol or a NII file or a MAT file)
% OUTPUT
%   vol = a 4D volume masked

if(isfield(cfg,'vol'))
	data=cfg.vol;
	% add check that it's a 4D vol
elseif(isfield(cfg,'infile'))
	nii=load_nii(cfg.infile);
	data=nii.img;
end

data=double(data);
kk=size(data);
T=kk(4);

% add here code to check if we have a mask
% make also sure that the mask has the same size as our data

if(isfield(cfg,'mask'))
	mask=cfg.mask;
else
	% if we don't have a mask, compute the mask
	mask = bramila_makemask(cfg);
end

fprintf('Applying mask...');

vol=zeros(kk);
for t=1:T
    vol(:,:,:,t)=mask.*data(:,:,:,t); % this can be made faster with reshape
end

fprintf(' done\n');

% if cfg.write==1 || nargout<1
%     bramila_savevolume(cfg,vol,'EPI volume after masking','masked.nii');
% end