function ts_roi = bspm_harvest_timeseries(images, roi, method)
% BSPM_HARVEST_TIMESERIES
%
%   USAGE: bspm_harvest_timeseries(images, roi, method)
% 
%   images = filenames from which to harvest
%   roi = filename of region of interest
%   method = 0: mean, >0: pca (value specifies number of dimensions)
%

% ----------------------- Copyright (C) 2014 -----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3, mfile_showhelp; return; end
if iscell(images); images = char(images); end;
if iscell(roi); roi = char(roi); end;

%% load in raw timeseries %%
hdr = spm_vol(images); ts = spm_read_vols(hdr);
dims = size(ts);

%% load in roi and reslice if necessary %%
roi_hdr = spm_vol(roi); roi = spm_read_vols(roi_hdr);
if dims(1:3)~=size(roi), roi = bspm_reslice(roi_hdr.fname, hdr(1).fname, 1, 1); end

%% extract signal from roi %%
roi = roi(:);
ts_rs = reshape(ts, prod(dims(1:3)), dims(4));
ts_roi = ts_rs(roi > 0,:);
if method
    % first eigenvariate from PCA
    [coeff, score, eigval] = pca(ts_roi');
    ts_roi = score(:,1:method);
else
    ts_roi = nanmean(ts_roi)';
end

end

%% SUBFUNCTION %%

function [out outmat] = bspm_reslice(in, ref, int, nowrite)
% BOB_RESLICE 
%
% USAGE: [out M] = bspm_reslice(in, ref, int, nowrite)
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
if nargin<2, display('USAGE: out = bspm_reslice(in, ref, int, nowrite)'); return; end
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
 
 
 
 
 
 
 
 
