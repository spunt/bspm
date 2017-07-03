function [scores labels]=bramila_compareModules(cfg)
%
% Usage:    cfg.reference='power' or 'cole' or 'yeo' or 'gordon'
%           cfg.infile=filename

switch lower(cfg.reference)
    case {'power'}
        refimg='power2011_modules.nii';
        power_labels;
    case {'cole'}
        refimg='cole2013_modules.nii';
        cole_labels;
    case {'yeo'}
        refimg='YEO_parc2mm_fsl.nii';
        yeo_labels;
    case {'gordon'}
        refimg='Parcels_MNI_222_modules.nii';
        gordon_labels;
end


img=load_nii(cfg.infile);

ref=load_nii(refimg);
ref.img=round(ref.img); % to make sure we have module ids
ref.img(find(ref.img<0))=0; % ignore negative indeces

voi=find(ref.img>0); %voxels of interests

%% go through each ground truth module

% first identify the positive module IDs
modids=unique(ref.img(:));
modids(find(modids<=0))=[];

for i=1:length(modids)
    temp=zeros(size(ref.img));
    temp(find(ref.img==modids(i)))=1;
    scores(i,1)=corr(double(temp(voi)),double(img.img(voi)));
    scores(i,2)=modids(i);
end

% normalized scores of overlap
allids=histc(ref.img(:),modids);
temp=img.img; % it could be masked?
scores(:,3)=histc(ref.img(find(temp(:)>0)),modids)./allids;



