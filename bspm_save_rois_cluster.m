function success = bspm_save_rois_cluster(in, heightThresh, sizeThresh, mask)
% BSPM_SAVE_CLUSTER
%
% USAGE: success = bspm_save_rois_cluster(in, heightThresh, sizeThresh, mask)
%
%   INPUTS
%       in:             image filename (full path if not in current dir)
%       heightThresh:   intensity threshold for defining clusters
%       sizeThresh:     extent threshold for defining clusters
%       mask:           optional mask filename (full path if not in 
%                       current directory)
%
%   OUTPUTS
%       success:        returns a 0 if no clusters could be identified
%                       after thresholding, 1 otherwise
%       

% ----------------------------- Copyright (C) 2014 -----------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3, disp('USAGE: success = bspm_save_rois_cluster(in, heightThresh, sizeThresh, mask)'); return; end
if nargin<4, mask = []; end
    
% make sure image names are character arrays
% ------------------------------------------------------
if iscell(in), in = char(in); end;
maskflag = 1;
if isempty(mask) || length(mask)<1, maskflag = 0; mask = 'no mask';
elseif iscell(mask), mask = char(mask); end;

% write current images to command window
% ------------------------------------------------------
bob_display_message('Looking for Cluster');
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
    heightThresh = bob_p2t(heightThresh, df);
end
if maskflag
    mask = bob_reslice(mask,in_hdr.fname,1,1);
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
	disp('No clusters found.')
    success = 0;
	return
end

% find index of voxels of sufficient size
% ------------------------------------------------------
cidx = unique(cl_index);
count = 0;
base_roi = zeros(imdims);
for i = cidx
    cluster_vox = voxels(:,cl_index==cidx(i));
    cluster_vox = cluster_vox';
    if length(cluster_vox)>=sizeThresh
        count = count + 1;
        cc = base_roi;
        for ii = 1:size(cluster_vox,1)
            cc(cluster_vox(ii,1),cluster_vox(ii,2),cluster_vox(ii,3)) = 1;
        end
        all_roi(:,:,:,count) = cc;
        all_size(count) = size(cluster_vox,1);
    end
end
if count==0
	disp('No clusters meet extent threshold.')
    success = 0;
	return
end
all_roi = double(all_roi);

% write the roi as an image
% ------------------------------------------------------
for i = 1:size(all_roi,4)
    roi_hdr = in_hdr;
    path = fileparts(in_hdr.fname);
    fulloutname = [path filesep 'ROI_Cluster' num2str(i) '_k=' num2str(all_size(i)) '.nii'];
    roi_hdr.fname = fulloutname;
    spm_write_vol(roi_hdr, all_roi(:,:,:,i));
    fprintf('\nROI of %d voxels written to:\n%s\n', all_size(i), [path filesep fulloutname]);
end
success = 1;
end
% -------------------------------------------------------------------------
% SUBFUNCTIONS
% -------------------------------------------------------------------------
function [out, outmat] = bob_reslice(in, ref, int, nowrite)
% BOB_RESLICE 
%
% USAGE: [out M] = bob_reslice(in, ref, int, nowrite)
%
% ARGUMENTS
%   in: path to image to reslice
%   ref: path to reference image (image to reslice to)
%   int: interpolation method, 0=Nearest Neighbor, 1=Trilinear(default)
%   nowrite: option to not write new volume (default = 0)
%
% OUTPUT
%   out: the resliced image volume
%
% Most of the code is adapted from rest_Reslice in REST toolbox:
% Written by YAN Chao-Gan 090302 for DPARSF. Referenced from spm_reslice.
% State Key Laboratory of Cognitive Neuroscience and Learning 
% Beijing Normal University, China, 100875
% --------------------------------------------------------------------------
if nargin<4, nowrite = 0; end
if nargin<3, int = 1; end
if nargin<2, display('USAGE: out = bob_reslice(in, ref, int, nowrite)'); return; end
if iscell(in); in = char(in); end
if iscell(ref); ref = char(ref); end

% read in reference image
RefHead = spm_vol(ref); 
RefData = spm_read_vols(RefHead);
mat=RefHead.mat;
dim=RefHead.dim;

% read in image to reslice
SourceHead = spm_vol(in);
SourceData = spm_read_vols(SourceHead);

% do the reslicing
[x1,x2,x3] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
d     = [int*[1 1 1]' [1 1 0]'];
C = spm_bsplinc(SourceHead, d);
v = zeros(dim);
M = inv(SourceHead.mat)*mat; % M = inv(mat\SourceHead.mat) in spm_reslice.m
y1   = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
y2   = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
y3   = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
out    = spm_bsplins(C, y1,y2,y3, d);

%Revised by YAN Chao-Gan 121214. Apply a mask from the source image: don't extend values to outside brain.
tiny = 5e-2; % From spm_vol_utils.c
Mask = true(size(y1));
Mask = Mask & (y1 >= (1-tiny) & y1 <= (SourceHead.dim(1)+tiny));
Mask = Mask & (y2 >= (1-tiny) & y2 <= (SourceHead.dim(2)+tiny));
Mask = Mask & (y3 >= (1-tiny) & y3 <= (SourceHead.dim(3)+tiny));
out(~Mask) = 0;
outmat = mat;
if ~nowrite
    OutHead=SourceHead;
    OutHead.mat      = mat;
    OutHead.dim(1:3) = dim;
    [p n e] = fileparts(SourceHead.fname);
    newname = sprintf('%s_%dx%dx%d%s',n,dim,e);
    OutHead.fname = [p filesep newname];
    spm_write_vol(OutHead,out);
end
end
function t = bob_p2t(alpha, df)
% BOB_P2T Get t-value from p-value + df
%
%   USAGE: t = bob_p2t(alpha, df)
%       
%   OUTPUT
%       t = crtical t-value
%
%   ARGUMENTS
%       alpha = p-value
%       df = degrees of freedom
%
% =========================================
if nargin<2, disp('USAGE: bob_p2t(p, df)'); return, end
t = tinv(1-alpha, df);
end








 
 
 
 
