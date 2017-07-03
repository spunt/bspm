function [rois] = bramila_makeRoiStruct(cfg)
% BRAMILA_MAKEROISTRUCT - Computes the ROI (Regions of interest) struct
% variable, used to define regions of interests in the BraMiLa and FUNPSY
% toolboxes.
%
%   - Usage: 
%       rois = bramila_makeRoiStruct(cfg)
%   - Input:
%       cfg.roimask  = Clustered ROI mask filename with full path.
%           This is a nifti files with integer indexes, zeros for where no
%           data is available
%       cfg.roi     = list ROI coordinates in MNI format. Used if cfg.roimask 
%           is not defined 
%       cfg.radius  = ROI sphere radius in voxels (default == 2.5 voxels)
%       cfg.labels  = a cell vector with each label for each roi
%       cfg.mirror  = mirror ROIs to opposing hemisphere and additionally
%           calculate averages for each mirrored pair (default: false)
%       cfg.imgsize = MRI image XYZ dimensions in voxels (default:[91 109 91])
%   - Output:
%       A struct vector with all the 


%       Example 1:
%           cfg.roi(1,:) = [10 0 0];
%           cfg.roi(2,:) = [-10 0 0];
%           cfg.radius = 3;
%           cfg.imgsize = [91 109 91];
%           makerois(cfg);
%
%       Example 2:
%           cfg.infile{1} = '/inputfiles/inputfile1.nii';
%           cfg.infile{2} = '/inputfiles/inputfile2.nii';
%           cfg.outpath = '/outputfiles/';
%           cfg.roi(1,:) = [10 0 0];
%           cfg.radius = 3;
%           cfg.mirror = true;
%           cfg.imgsize = [91 109 91];
%           makerois(cfg);
%
% Last edit: EG 2014-01-23



%% ROI & Image Data
% ROIs
if ~isfield(cfg,'roimask')
    if ~isfield(cfg,'roi')
        error('No ROIs defined (cfg.roi or cfg.roimask required).');
    end
end
% ROI sphere radius
if ~isfield(cfg,'radius')
    cfg.radius = 2.5;
    warning('BRAMILA:noROIRadius',...
        'No ROI sphere radius (cfg.radius) defined. \nUsing default value of 2.5 voxels');
end

% Image dimensions
if ~isfield(cfg,'imgsize')
    cfg.imgsize = [91 109 91];
    warning('BRAMILA:noImgSize',...
        'No image size (cfg.imgsize) defined. \nUsing default value of [91 109 91]');
end



% Mirror mode
if ~isfield(cfg,'mirror')
    cfg.mirror = false;
end


%% Process files
%

if isfield(cfg,'roimask')
    mask = load_nii(cfg.roimask);
    % make sure that mask matches the data size
    if(any(~(size(mask.img) == cfg.imgsize)))
        error('Mask and imgsize do not match');
    end
    
    % Get ROIs coordinates
    k = 0;
    uID=unique(mask.img(:));
    if(uID(1)==0)
        uID(1)=[];
    end
    for r = 1:length(uID)
        hits = find(mask.img==uID(r));
        [x, y, z] = ind2sub(cfg.imgsize,hits);
        rois(r).map = [x y z];
        rois(r).centroid=round(mean(rois(r).map,1));
    end
else    
    roicoo=cfg.roi;
    if cfg.mirror
        ids=find(roicoo(:,1)~=0);
        roicoo=[roicoo;
            -roicoo(ids,1) roicoo(ids,2) roicoo(ids,3) ];
    end
    for r=1:size(roicoo,1)
        temp=roicoo(r,:);
		cfgtemp=[];
		cfgtemp.coordinates=[temp(1) temp(2) temp(3)];
		cfgtemp.type='MNI';
		cfgtemp.imgsize=cfg.imgsize;
        %[x, y, z]=bramila_MNI2space(temp(1), temp(2), temp(3));
        [x, y, z]=bramila_MNI(cfgtemp);
        rois(r).centroid=[x y z];
    end
        
    A=ones(cfg.imgsize(1),cfg.imgsize(2),cfg.imgsize(3));
    [x y z]=ind2sub(size(A),find(A==1));
    v=[x y z];

    for n=1:length(rois);
        vt=v-repmat(rois(n).centroid,size(v,1),1);
        vt=vt';
        mapID=find(sum((vt).^2)<=cfg.radius^2);
        map=v(mapID',:);
        rois(n).map=map;
    end
end
for n=1:length(rois);
    rois(n).label='';
    if(n<=length(cfg.labels))
        rois(n).label=cfg.labels{n};
    end
	cfgtemp=[];
	cfgtemp.coordinates=[rois(n).centroid(1),rois(n).centroid(2),rois(n).centroid(3)];
	cfgtemp.type='space';
	cfgtemp.imgsize=cfg.imgsize;
    %[mnix mniy mniz]=bramila_space2MNI(rois(n).centroid(1),rois(n).centroid(2),rois(n).centroid(3));
    [mnix mniy mniz]=bramila_MNI(cfgtemp);
    rois(n).centroidMNI=[mnix mniy mniz];
end
