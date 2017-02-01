function [data,mask] = bramila_makemask(cfg)
    % INPUT
    %   cfg.vol = input can be a 4D volume, if this exists, then infile is not loaded but used as ID (optional)
    % OUTPUT
    %   vol = a 3D volume with the mask
    % by default the mask is saved, add the save volume function 
    
    % we need to specify some options if we want the mask based on Power of time series or on amplitude at each volume
    % to check the source to reference for the 'spm' stuff. SPM recent code seems different too            
    
if(isfield(cfg,'vol'))
	data=cfg.vol;
	% add check that it's a 4D vol
elseif(isfield(cfg,'infile'))
	nii=load_nii(cfg.infile);
	data=nii.img;
end

kk=size(data);
T=kk(4);
% for now we implement only the 'spm?' case. Enrico to copy the code from funpsy for the other case

fprintf('Computing EPI mask...');

mask =  ones(kk(1:3));
for t=1:T  
    temp=squeeze(data(:,:,:,t));
    mask=mask.*(temp>0.1*quantile(temp(:),.98));
end

fprintf(' done\n');

if isfield(cfg,'mask') && ~isempty(cfg.mask)
    if all(size(cfg.mask)==size(mask))        
        N1=nnz(cfg.mask);
        mask = mask.*double(cfg.mask);
        N2=nnz(mask);

        fprintf('...EPI mask size %i voxels (from %i)\n',N2,N1);
    end
else
   fprintf('...EPI mask size %i voxels\n',nnz(mask)); 
end

if nnz(mask)==0
   warning('Empty EPI mask!');
end

if nnz(mask)>0.5*numel(mask)
   warning('EPI mask has over 50% of total FOV volume!');
end

%if cfg.write==1 || nargout<1
bramila_savevolume(cfg,mask,'EPI mask','EPI_mask.nii');
%end
