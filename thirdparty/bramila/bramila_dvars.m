function [dvars,img]=bramila_dvars(cfg)
% BRAMILA_DVARS - Computes Derivative VARiance across voxels as defined in
% Power et al. (2012) doi:10.1016/j.neuroimage.2011.10.018
%   - Usage:
%   dvars=bramila_dvars(cfg) Returns a time series 'dvars' with a value of
%   RMS for each time point. First time point is set to 0.
%   - Input:
%   cfg is a struct with following parameters
%       Possible input formats
%       cfg.infile = 'path/to/a/nifti/file' - insert the full path to a nifti
%           file with 4D fMRI data
%       cfg.vol = vol - a matlab 4D volume with fMRI data, time on the
%           4th dimension
%       cfg.ts = ts - a two dimensional vector of time series, time on the
%           1st dimension
%       cfg.plot = 0 or 1 - set to 1 if you want to output a plot like in
%           Power et al. (2014) doi:10.1016/j.neuroimage.2013.08.048 (defult 0)
%       cfg.mask = a 3D matlab volume mask of voxels to consider for RMS computation
%           (the parameter is ignored if cfg.ts is specified)
%   - Note:
%   if more than one input format is specified matlab will give
%   priority to cfg.ts > cfg.vol > cfg.nii
%
%   Last edit: EG 2014-01-10

%% loading the data

fprintf('Computing DVARS...');

% load data
if isfield(cfg,'vol') && ~isempty(cfg.vol)
    img=double(cfg.vol);
    reshape_vol=1;
elseif isfield(cfg,'infile') && ~isempty(cfg.infile)
    nii=load_nii(cfg.infile);
    img=double(nii.img);
    reshape_vol=1; 
else
    error('No input EPI data found!');
end

plotIt=0;
if(isfield(cfg,'plot'))
    plotIt=cfg.plot;
end

% if we have a mask and if we have 4D data, then apply the mask
hasmask=0;
if (isfield(cfg,'analysis_mask') && size(img,4)>0)
    mask=double(cfg.analysis_mask);
    sz=size(img);
    if(~any(size(mask) ==sz(1:3)))
        error(['The specified mask has a different size than the fMRI data. Quitting.'])
    end
    for t=1:size(img,4)
        img(:,:,:,t)=mask.*img(:,:,:,t);
    end
    hasmask=1;
elseif (isfield(cfg,'mask') && size(img,4)>0)
    mask=double(cfg.mask);
    sz=size(img);
    if(~any(size(mask) ==sz(1:3)))
        error(['The specified mask has a different size than the fMRI data. Quitting.'])
    end
    for t=1:size(img,4)
        img(:,:,:,t)=mask.*img(:,:,:,t);
    end
    hasmask=1;
end

% doing a reshape if needed
if(reshape_vol==1)
    sz=size(img);
    T=sz(4);
    img=reshape(img,[],T);
    img=img';   % Time in first dimension
end

if(isfield(cfg,'ts'))
    img=double(cfg.ts);
    T=size(img,1);
end

img=bramila_bold2perc(img); % convert bold to percentage signal
di=diff(img);
di=[
    zeros(1,size(di,2)) % adding a zero as first sample of the derivate
    di
    ];
    
if hasmask==1
    mask_ids=find(mask>0);
    dvars=sqrt(mean(di(:,mask_ids).^2,2)); % Root Mean Square across voxels
else
    dvars=sqrt(mean(di.^2,2)); % Root Mean Square across voxels
end

if(plotIt~=0)
    % plot could be improved by separating gray matter, white matter, etc
    % as in Power et al 2014
    figure;
    imagesc(img',[-2 2]);
    colorbar
    xlabel('Time [samples]')
    ylabel('Voxels')
    colormap(gray)
end

fprintf(' done\n');