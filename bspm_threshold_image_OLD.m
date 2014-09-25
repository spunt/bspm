function out = bspm_threshold_image_OLD(in, heightThresh, sizeThresh, binflag, outname)
% BSPM_THRESHOLD_IMAGE
%
% USAGE: out = bspm_threshold_image(in, heightThresh, sizeThresh, binflag, outname)
%
%   ARGUMENTS
%       in:                     input image name (full path if not in current directory)
%       heightThresh:   intensity threshold f(if < .10, will assume it is
%                               an alpha level and will covert to critical t)
%       sizeThresh:      extent threshold for defining clusters
%       binflag:            flag to binarize output image (default = 0)
%       outname:         name for file to write (default = no file written)
%       
% Created January 1, 2013 - Bob Spunt

% ---------------------------------- Copyright (C) 2014 ----------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<5, outname = []; end; 
if nargin<4, binflag = 0; end;
if nargin<3, disp('USAGE: out = bspm_threshold_image(in, heightThresh, sizeThresh, binflag, outname)'); return; end
if iscell(in), in = char(in); end;

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
    heightThresh = bob_p2t(heightThresh, df);
end
in(in<heightThresh) = NaN;
in(in==0)=NaN;

% grab voxels
% ------------------------------------------------------
[X Y Z] = ind2sub(size(in), find(in > 0));
voxels = sortrows([X Y Z])';

% get cluster indices of voxels
% ------------------------------------------------------
cl_index = spm_clusters(voxels);

% find index of clusters of sufficient size
% ------------------------------------------------------
for i = 1:max(cl_index)
    a(cl_index == i) = sum(cl_index == i);
end
which_vox = (a >= sizeThresh);
cluster_vox = voxels(:,which_vox);
cluster_vox = cluster_vox';
roi_mask = zeros(imdims);
for i = 1:size(cluster_vox,1)
    roi_mask(cluster_vox(i,1),cluster_vox(i,2),cluster_vox(i,3)) = in(cluster_vox(i,1),cluster_vox(i,2),cluster_vox(i,3));
end
out = double(roi_mask);
out(out==0) = NaN;
if binflag, out = out>0; else out(out==0) = NaN; end

if ~isempty(outname)
    
    h = in_hdr;
    h.fname = outname;
    tmp = h.descrip;
    h.descrip = [tmp sprintf(' - THRESHOLDED: %2.3f, %d', heightThresh, sizeThresh)];
    spm_write_vol(h,out);
    
end


end






 
 
 
 
