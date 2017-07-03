function [blacklist,roimask,new_rois] = bramila_maskrois(rois,mask,limit_percentage)
% bramila_create_blacklist.m, returns a list of ROIs that are outside given mask.
% IN:
% rois = roi structure (contains at least a field 'map')
% mask = 3D mask of interest (values 0 are outside mask)
% limit_percentage = roi must contain at least this many voxels compared to roi definition, range (0,1] (OPTIONAL)
% OUT:
% blacklist = list of bad rois
% roimask = mask of good rois with each roi uniquely numbered (OPTIONAL)
% new_rois = modified roi structure with only good rois remaining (OPTIONAL)

% has much voxels must remain inside roi (0.50 = 50%), otherwise roi is marked as bad
if nargin<3
    limit_percentage = 0.5;
else
   if limit_percentage<=0 || limit_percentage>1
      error('Limit percentage must be between (0,1]') 
   end
end

N=length(rois);

clustersizes = zeros(1,N);
clustersizes_full = zeros(1,N);

roimask = 0*mask;

for clust = 1:N    
    good = 0;        
    coord = rois(clust).map;    
    clustersizes_full(clust)=size(coord,1);    
    for i=1:size(coord,1)        
        roimask(coord(i,1),coord(i,2),coord(i,3))=clust;
        if mask(coord(i,1),coord(i,2),coord(i,3))~=0
            good = good + 1;
        end
    end    
    clustersizes(clust)=good;        
end

bad_clusters = find(clustersizes<floor(clustersizes_full*limit_percentage));
good_clusters = find(clustersizes>=floor(clustersizes_full*limit_percentage));

if length(bad_clusters)+length(good_clusters)~= length(clustersizes)
    error('Array size mismatch (BUG FOUND)!')
end

blacklist = bad_clusters;

if length(blacklist)>N/2
    warning('Over 50%% bad ROIs, is your mask valid?')
end

new_rois=rois;
new_rois(blacklist)=[];

end

