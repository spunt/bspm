function [reg_cfg,cfg] = bramila_maketissueregress(cfg)

max_dim = 25;
labels = [];
X_white = [];
X_csf =[];
X_grey = [];
id = [];

if cfg.remove_global == 1
    fprintf('...global mean signal\n')
    X_grey = compute_mean(cfg.vol,cfg.mask);
    labels{1} = 'mean brain signal inside EPI mask';
    id = cat(2,id,1);
end

if cfg.max_tissue_pca_count>0
    
    new_cfg.vol  = cfg.vol;
        
    [white,csf,cfg] = bramila_get_tissue_rois(cfg);
    
    fprintf('...extracting PCA timeseries ');
    
    dim=0;
    perc=0;
    if ~isempty(white.map) && size(white.map,1)>4
  	 new_cfg.max_PCA_count = min(size(white.map,1),max_dim);
        new_cfg.rois = white;
        [X_white,perc,comps] = bramila_roiextract(new_cfg);
        
        dim = min([comps{1},cfg.max_tissue_pca_count,find(perc{1}>=cfg.min_tissue_var_expl,1,'first')]);
        cfg.PCA_variance_WM = perc{1}(dim);
        perc=cfg.PCA_variance_WM;
        X_white = X_white{1}(:,1:dim);
        labels = cat(2,labels,get_tissue_labels('WM',dim));
        id = cat(2,id,2*ones(1,dim));
    end
    fprintf('( WM=%i,%3.1f%%, ',dim,perc);
    
    dim=0;
    perc=0;
    if ~isempty(csf.map) && size(csf.map,1)>4
        new_cfg.max_PCA_count = min(size(csf.map,1),max_dim);
        new_cfg.rois = csf;
        [X_csf,perc,comps]  = bramila_roiextract(new_cfg);
        
        dim = min([comps{1},cfg.max_tissue_pca_count,find(perc{1}>=cfg.min_tissue_var_expl,1,'first')]);
        cfg.PCA_variance_CSF = perc{1}(dim);
        perc=cfg.PCA_variance_CSF;
        X_csf = X_csf{1}(:,1:dim);
        labels = cat(2,labels,get_tissue_labels('CSF',dim));
        id = cat(2,id,3*ones(1,dim));
    end
    fprintf('CSF=%i,%3.1f%% )\n',dim,perc);
    
else
   cfg.analysis_mask = cfg.mask; 
end

X_nuisance = [X_grey,X_white,X_csf];

X=zscore(X_nuisance);  % 6 columns x T timepoints
labels_der = [];
id_der=[];
for i=1:cfg.tissue_derivatives
    X_nuisance=cat(2,X_nuisance,bramila_derivative(X,i));
    labels_der = cat(2,labels_der,get_deriv_labels(labels,i));
    id_der = cat(2,id_der,id);
end

labels_total = cat(2,labels,labels_der);
id = cat(2,id,id_der);

reg_cfg.reg = X_nuisance;
reg_cfg.id = id;
reg_cfg.vol = cfg.vol;
reg_cfg.labels = labels_total;

end

function labels = get_deriv_labels(labels,order)

dot = repmat('''',1,order);
for i=1:length(labels)
    labels{i}=[labels{i},dot];
end
end

function labels = get_tissue_labels(tissuename,dim)
for i=1:dim
    labels{i}=sprintf('%s PCA-%i',tissuename,i);
end
end

% function ts = get_roi_ts(cfg)
%
%     ts = zeros(size(cfg.vol,4),size(cfg.roi.map,1));
%
%     for i=1:size(cfg.roi.map,1)
%         a= cfg.vol(cfg.roi.map(i,1),cfg.roi.map(i,2),cfg.roi.map(i,3),:);
%         ts(:,i)=squeeze(a);
%     end
%
% end

function [white,csf,cfg] = bramila_get_tissue_rois(cfg)

volsize=size(cfg.vol);

%--Load masks-----------------

if cfg.white_mask_th<0.6
    warning('!! Very large WM mask (th<0.6), a high risk of good signal removal !!')
end
nii=load_nii(cfg.white_mask);
white_mask=0*double(nii.img);
white_mask(nii.img>cfg.white_mask_th)=1;
white_mask(isnan(white_mask))=0;
white_mask=double(white_mask~=0); % convert to mask (if not already so)
if nnz(white_mask)>numel(white_mask)*0.20
    warning('!! Very large WM mask (>20% of volume), a high risk of good signal removal !!')
end

if cfg.csf_mask_th<0.6
    warning('!! Very large CSF mask (th<0.6), a high risk of good signal removal !!')
end
nii=load_nii(cfg.csf_mask);
csf_mask=0*double(nii.img);
csf_mask(nii.img>cfg.csf_mask_th)=1;
csf_mask(isnan(csf_mask))=0;
csf_mask=double(csf_mask~=0); % convert to mask (if not already so)
if nnz(csf_mask)>numel(csf_mask)*0.10
    warning('!! Very large CSF mask (>10% of volume), a high risk of good signal removal !!')
end

%--Double check mask size-----------------

if ~all(size(white_mask)==volsize(1:3))
    error('WHITE mask size does not match EPI size!')
end
if ~all(size(csf_mask)==volsize(1:3))
    error('CSF mask size does not match EPI size!')
end

csf_mask=csf_mask.*(~white_mask); % remove any overlap between WHITE and CSF masks

%--Apply ISC mask-----------------

if cfg.do_spatial_ISC>0 && isfield(cfg,'ISC_mask') % only for task-related data
    fprintf('...ISC mask ')
    nii=load_nii(cfg.ISC_mask);
    ISC_mask = double(nii.img);
    if ~all(size(ISC_mask)==volsize(1:3))
        error('ISC mask size does not match EPI size!')
    end
    fprintf('(%i voxels)\n',nnz(ISC_mask));
    
    N1=nnz(white_mask);
    white_mask=white_mask.*(~ISC_mask);
    N11=nnz(white_mask);
    N2=nnz(csf_mask);
    csf_mask=csf_mask.*(~ISC_mask);
    N22=nnz(csf_mask);
    
    if N11<N1 || N22<N2
        fprintf('...ISC mask overlap: WM reduced by %i voxels (%i->%i) and CSF by %i voxels (%i->%i)\n',N1-N11,N1,N11,N2-N22,N2,N22);
    end
    
end

if ~isfield(cfg,'mask')
    cfg.mask = ones(volsize(1:3));
    warning('No mask found (you should have an initial EPI mask at this point!). Creating all zero mask.')
end

N1=nnz(cfg.mask);

% nullify masks outside EPI mask
csf_mask(cfg.mask==0)=0;
white_mask(cfg.mask==0)=0;

% remove nuisance regions from EPI mask, creates analysis mask
cfg.analysis_mask=cfg.mask;
cfg.analysis_mask(csf_mask~=0)=0;
cfg.analysis_mask(white_mask~=0)=0;

% create an empty layer between analysis mask and nuisance mask ("eroding")
cfg.analysis_mask = bramila_erodemask(cfg.analysis_mask,csf_mask+white_mask);

N2=nnz(cfg.analysis_mask);

fprintf('...analysis mask size %i voxels (from %i)\n',N2,N1);
fprintf('...WM mask size %i voxels\n',nnz(white_mask));
fprintf('...CSF mask size %i voxels\n',nnz(csf_mask));

white = makeroi(white_mask);
white.label = 'white matter mask';
bramila_savevolume(cfg,white_mask,'white matter mask','wm_mask.nii');

csf = makeroi(csf_mask);
csf.label = 'cerebrospinal fluid mask';
bramila_savevolume(cfg,csf_mask,'csf mask','csf_mask.nii');

end

function roi = makeroi(mask)

[x,y,z] = ind2sub(size(mask),find(mask));
if isempty(x)
    roi.map = [];
else
    roi.map = [x,y,z];
end

end

function X = compute_mean(data,mask)
ind = mask~=0;
T=size(data,4);
X = zeros(T,1);
for t=1:T
    vol = squeeze(data(:,:,:,t));
    X(t)=mean(vol(ind));
end
end
