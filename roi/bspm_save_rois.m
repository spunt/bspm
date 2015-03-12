function bspm_save_rois(in, thresh, roi, name)
% BSPM_SAVE_ROIS
%
%   USAGE: bspm_save_rois(in, thresh, roi, name)
%
%   ARGUMENTS
%
%       in = reference image
%
%       thresh.cluster = [intensity extent];
%       thresh.peak = [intensity extent];
%       thresh.separation = minimum peak separation
%
%       roi.shape = 'Sphere' or 'Box'
%       roi.size = radius of ROI
%
%       name = append to file (default = name of input filename)
%       
% Created April 8, 2013 - Bob Spunt
% Uses code authored by:
% Dr. Robert Welsh (SimpleROIBuilder.m)
% Drs. Donald McLaren & Aaron Schultz (peak_nii.m)

% --------------------- Copyright (C) 2014 ---------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<3, error('USAGE: bspm_save_rois(in, thresh, roi, name)'); end
if iscell(in), in = char(in); end
if nargin<4, [p, name] = fileparts(in); name = upper(regexprep(name,' ','_')); end

% pull out fields
refimage = in;

% get peak information from peak_nii
% ------------------------------------------------------
% default structure input to peak_nii
peaknii.thresh = thresh.peak(1);
peaknii.cluster = thresh.peak(2);
peaknii.out = '';
peaknii.sign = 'pos';
peaknii.type = 'T';
peaknii.voxlimit = [];
peaknii.separation = thresh.separation;
peaknii.SPM = 1;
peaknii.conn = [];
peaknii.mask = [];
peaknii.df1 = [];
peaknii.df2 = [];
peaknii.nearest = 1;
peaknii.label = 'aal_MNI_V4';

% run peak_nii
[voxels] = peak_nii(refimage,peaknii);

% get names and coordinates of peaks
names = voxels{2};
names = regexprep(names,',','');
names = regexprep(names,' ','_');
roi_coords = voxels{1};
roi_coords = roi_coords(:,3:5);


% load reference image and threshold
% ------------------------------------------------------
P_HDR = spm_vol(refimage);
refMASK = bspm_threshold_image(refimage,thresh.cluster(1),thresh.cluster(2),1);

% create empty image and header for ROIs
% ------------------------------------------------------
maskIMG = zeros(P_HDR.dim(1:3));
maskHDR = P_HDR;
maskHDR.dim = [P_HDR.dim];     % Make it uint8 since it is binary.
maskHDR.pinfo = [1;0;0];
maskHDR.mat = P_HDR.mat;
roiINFO = {};
nROIS = length(names);

for iROI = 1:nROIS
    roiINFO{iROI}.center_mm = roi_coords(iROI,:)';
    tmp = inv(P_HDR.mat)*([roiINFO{iROI}.center_mm; 1]);
    roiINFO{iROI}.center_vox = tmp(1:3);
    roiINFO{iROI}.type = roi.shape;
    roiINFO{iROI}.size = roi.size;
end

% build an array of coordinates for each and every voxel in the mask
xOrds = (1:P_HDR.dim(1))'*ones(1,P_HDR.dim(2));
yOrds = ones(P_HDR.dim(1),1)*(1:P_HDR.dim(2));
xOrds = xOrds(:)';
yOrds = yOrds(:)';
Coords = zeros(3,prod(P_HDR.dim(1:3)));
for iZ = 1:P_HDR.dim(3)
    zOrds = iZ*ones(1,length(xOrds));
    Coords(:,(iZ-1)*length(xOrds)+1:iZ*length(xOrds)) = [xOrds; yOrds; zOrds];
end

% Now put them into mm's
mmCoords = P_HDR.mat*[Coords;ones(1,size(Coords,2))];
mmCoords = mmCoords(1:3,:);
boxBIT = zeros(4,size(mmCoords,2));

% Now loop on the ROI definitions and drop them
% into the mask image volume matrix.
[refpath, refname, refext] = fileparts(maskHDR.fname);
if isempty(refpath), refpath = pwd; end

for iROI = 1:nROIS
    % Found the center of this ROI in voxels
    % and then build it.
    cx = roiINFO{iROI}.center_mm(1);
    cy = roiINFO{iROI}.center_mm(2);
    cz = roiINFO{iROI}.center_mm(3);
    xs = mmCoords(1,:) - roiINFO{iROI}.center_mm(1);
    ys = mmCoords(2,:) - roiINFO{iROI}.center_mm(2);
    zs = mmCoords(3,:) - roiINFO{iROI}.center_mm(3);
    cROI = zeros(size(maskIMG));
    cHDR = maskHDR;
    switch lower(roiINFO{iROI}.type)
        case 'sphere'
            radii = sqrt(xs.^2+ys.^2+zs.^2);
            VOXIdx = find(radii<=roiINFO{iROI}.size);
        case 'box'
            xsIDX = find(abs(xs)<=roiINFO{iROI}.size(1));
            ysIDX = find(abs(ys)<=roiINFO{iROI}.size(1));
            zsIDX = find(abs(zs)<=roiINFO{iROI}.size(1));
            boxBIT  = 0*boxBIT;
            boxBIT(1,xsIDX) = 1;
            boxBIT(2,ysIDX) = 1;
            boxBIT(3,zsIDX) = 1;
            boxBIT(4,:) = boxBIT(1,:).*boxBIT(2,:).*boxBIT(3,:);
            VOXIdx = find(boxBIT(4,:));
    end
    cROI(VOXIdx) = 1;
    cROI = cROI.*refMASK; % intersect with thresholded image
    
    csize = sum(cROI(:));
    
    cHDR.fname = [refpath filesep 'ROI_' name '_' names{iROI} '_' num2str(cx) '_' num2str(cy) '_' num2str(cz) '_k='  ... 
        num2str(csize) '_' upper(roi.shape) num2str(roi.size) '_SEP' num2str(thresh.separation) '.nii'];

    
    
    
    spm_write_vol(cHDR,cROI);
end
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
function out = bspm_threshold_image(in, heightThresh, sizeThresh, binflag, outname)
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
% ========================================================================%
if nargin<5, outname = []; end; 
if nargin<4, binflag = 0; end;
if nargin<3, disp('USAGE: out = bspm_threshold_image(in, heightThresh, sizeThresh, binflag, outname)'); return; end
if iscell(in), in = char(in); end;


% load images
% ------------------------------------------------------
try
    in_hdr = spm_vol(in);
    in = spm_read_vols(in_hdr);
catch
    in = in;
end
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


 
 
 
 
