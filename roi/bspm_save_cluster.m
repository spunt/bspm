function success = bspm_save_cluster(in, outname, heightThresh, sizeThresh, mask, nameflag)
% BSPM_SAVE_CLUSTER
%
% USAGE: bspm_save_cluster(in, out, heightThresh, sizeThresh, mask)
%
%   ARGUMENTS
%       in: input image name (full path if not in current directory)
%       outname: base name for output image (written in same directory
%       as the input image)
%       heightThresh: intensity threshold for defining clusters
%       sizeThresh: extent threshold for defining clusters
%       mask: mask image to use
%       nameflag: if 1, number of voxels will be appended to written roi
%       
% Created January 1, 2013 - Bob Spunt

% ------------------------------------ Copyright (C) 2014 ------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<6, nameflag = 0; end
if nargin<5, mask = []; end
if nargin<4, disp('USAGE: bspm_save_cluster(in, outname, heightThresh, sizeThresh, mask, nameflag)'); return; end
    
% make sure image names are character arrays
% ------------------------------------------------------
if iscell(in), in = char(in); end;
maskflag = 1;
if isempty(mask) || length(mask)<1, maskflag = 0; mask = 'no mask';
elseif iscell(mask), mask = char(mask); end;

% write current images to command window
% ------------------------------------------------------
bspm_display_message('Looking for Cluster');
fprintf('\nSource Image: %s\nMask Image: %s\nnHeight Threshold: %d\nSizeThreshold: %d\n', in, mask, heightThresh, sizeThresh);

% load images
% ------------------------------------------------------
in_hdr = spm_vol(in);
in = spm_read_vols(in_hdr);
imdims = size(in);
% if necessary, calculate critical t
if ismember(heightThresh,[.10 .05 .01 .005 .001 .0005 .0001]);
    tmp = in_hdr.descrip;
    idx1 = regexp(tmp,'[','ONCE');
    idx2 = regexp(tmp,']','ONCE');
    df = str2num(tmp(idx1+1:idx2-1));
    heightThresh = bspm_p2t(heightThresh, df);
end
if maskflag
    mask = bspm_reslice(mask,in_hdr.fname,1,1);
else
    mask = in>0;
end


% apply mask and height threshold
% ------------------------------------------------------
in(mask==0) = NaN;
in(in<heightThresh) = NaN;

% grab voxels
% ------------------------------------------------------
[X Y Z] = ind2sub(size(in), find(in > 0));
voxels = sortrows([X Y Z])';

% get cluster indices of voxels
% ------------------------------------------------------
cl_index = spm_clusters(voxels);
if isempty(cl_index)
	disp('No clusters meet extent threshold.')
    success = 0;
	return
end

% find index of voxels of sufficient size
% ------------------------------------------------------
for i = 1:max(cl_index)
    a(cl_index == i) = sum(cl_index == i);
end
which_vox = (a >= sizeThresh);
numClusters = sum(length(unique(cl_index(find(a >= sizeThresh)))));
if numClusters == 0
	disp('No clusters meet extent threshold.')
    success = 0;
	return
end

% get voxels for the largest cluster and convert to mask
% ------------------------------------------------------
cluster_vox = voxels(:,a==max(a));
cluster_vox = cluster_vox';
roi_mask = zeros(imdims);
for i = 1:size(cluster_vox,1)
    roi_mask(cluster_vox(i,1),cluster_vox(i,2),cluster_vox(i,3)) = 1;
end
roi_mask = double(roi_mask);

% write the roi as an image
% ------------------------------------------------------
roi_hdr = in_hdr;
[path name ext] = fileparts(in_hdr.fname);
[path2 name2 ext2] = fileparts(outname);
if nameflag==1
    fulloutname = [name2 '_' num2str(max(a)) 'voxels.nii'];
else
    fulloutname = [name2 '.nii'];
end
roi_hdr.fname = [path filesep fulloutname];
spm_write_vol(roi_hdr, roi_mask);

% write message to command window
fprintf('ROI of %d voxels written to:\n%s\n\n', max(a), [path filesep fulloutname]);
success = 1;


end






 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
