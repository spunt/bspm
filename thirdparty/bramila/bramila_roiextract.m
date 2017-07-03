function [nodeTS,perc,total_comp] = bramila_roiextract(cfg)
% INPUT
%   cfg.rois % see previous roi files
%   cfg.infile = location where the subject NII file (4D)
%   cfg.vol = input can be a 4D volume, if this exists, then infile is not loaded but used as ID (optional)
%   cfg.write = 0 (default 0, set to 1 if you want to store the output)
%   cfg.verbose = 1 (default 0, set 1 if you want feedback on the matlab window
%   cfg.max_PCA_count = number of PCA components requested (OPTIONAL)
% 	cfg.usemean = 1, does the mean instead of PC (optional)
% OUTPUT
%   nodeTS = matrix of time series for each roi (time in 1st dimension)
%   perc = percentage of variance explained (OPTIONAL)



if isfield(cfg,'vol') && ~isempty(cfg.vol)
    data=cfg.vol;
    % add check that it's a 4D vol
elseif isfield(cfg,'infile') && ~isempty(cfg.infile)
    nii=load_nii(cfg.infile);
    data=nii.img;
else
    error('No volume or infile found!')
end

hasmask=0;
if isfield(cfg,'mask') && ~isempty(cfg.mask)
	mask=cfg.mask;
	hasmask=1;
end

if(~isfield(cfg,'verbose'))
    cfg.verbose=0;
end

if(~isfield(cfg,'usemean'))
    cfg.usemean=0;
end



% Extracting ROI timeseries
rois=cfg.rois;
R=length(rois);
T=size(data,4);

if isfield(cfg,'max_PCA_count') && cfg.max_PCA_count>1
    
    nodeTS=cell(R);
    perc = cell(R);
    total_comp = cell(R);
    kk=0;
    for r=1:R % go through each ROI
        if(cfg.verbose) disp(['Extracting node ' num2str(r)]); end
        map=rois(r).map; % store the voxel coordinates
        ts=zeros(T,size(map,1));
        k=0;
        for m=1:size(map,1) % for each voxel in ROI
            c=map(m,:); % take the coordinates of this voxel
            if hasmask==0 || (hasmask==1 && mask(c(1),c(2),c(3))~=0)
                k=k+1;
                ts(:,k)=squeeze(data(c(1),c(2),c(3),:)); % take the time series of this voxel and store it in ts
            end
        end
        ts=ts(:,1:k);
        if k>0
            kk=kk+1;
            [nodeTS{r},perc{r},total_comp{r}]=bramila_princomp(ts,cfg.max_PCA_count);
        else
            nodeTS{r}=nan;
            perc{r}=nan;
        end
        
    end
    
    if kk<R && R>1
       fprintf('..total %i of %i rois (partially) inside mask\n',kk,R);
    end
    
else % take max_PCA_count==1
    
    perc = [];
    total_comp = [];
    nodeTS=zeros(T,R);
    kk=0;
    for r=1:R % go through each ROI
        if(cfg.verbose) disp(['Extracting node ' num2str(r)]); end
        map=rois(r).map; % store the voxel coordinates
        ts=zeros(T,size(map,1));
        k=0;
        for m=1:size(map,1) % for each voxel in ROI
            c=map(m,:); % take the coordinates of this voxel
			if hasmask==0 || (hasmask==1 && mask(c(1),c(2),c(3))~=0)
                k=k+1;
                ts(:,k)=squeeze(data(c(1),c(2),c(3),:)); % take the time series of this voxel and store it in ts
            end
        end
        ts=ts(:,1:k);
        
        if k>0
            kk=kk+1;
			if(cfg.usemean>0)
				nodeTS(:,r)=mean(ts,2);
			else
				if(k>125) disp('>>> I would consider using cfg.usemean = 1, rather than 1st principal component'); end
				nodeTS(:,r)=bramila_princomp(ts); % take the 1st PC from ROI
			end
        else
			%if we have no voxels
            nodeTS(:,r)=nan; 
        end
    end
    
    if kk<R && R>1
       fprintf('..total %i of %i rois (partially) inside mask\n',kk,R);
    end
    
end

end
