function [out, outmat] = bspm_reslice(in, ref, int, nowrite)
% BSPM_RESLICE 
%
% USAGE: [out, outmat] = bspm_reslice(in, ref, int, nowrite)
%
% ARGUMENTS
%   in:         path to image to reslice
%   ref:        path to reference image (image to reslice to)
%   int:        interpolation, 0=Nearest Neighbor, 1=Trilinear(default)
%   nowrite:    option to not write new volume (default = 0)
%
% OUTPUT
%   out:        the resliced image volume
%
% Most of the code is adapted from rest_Reslice in REST toolbox:
% Written by YAN Chao-Gan 090302 for DPARSF. Referenced from spm_reslice.
% State Key Laboratory of Cognitive Neuroscience and Learning 
% Beijing Normal University, China, 100875
% --------------------------------------------------------------------------
if nargin<2, display('USAGE: out = bob_reslice(in, ref, int, nowrite)'); return; end
if nargin<4, nowrite = 0; end
if nargin<3, int = 1; end
if iscell(in); in = char(in); end
if iscell(ref); ref = char(ref); end

% read in reference image
RefHead     = spm_vol(ref); 
mat         = RefHead.mat;
dim         = RefHead.dim;
SourceHead  = spm_vol(in);

% do the reslicing
[x1,x2,x3]  = ndgrid(1:dim(1),1:dim(2),1:dim(3));
d           = [int*[1 1 1]' [1 1 0]'];
C       = spm_bsplinc(SourceHead, d);
M       = SourceHead.mat\mat; % M = inv(mat\SourceHead.mat) in spm_reslice.m
y1      = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
y2      = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
y3      = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
out     = spm_bsplins(C, y1,y2,y3, d);

%Revised by YAN Chao-Gan 121214. Apply a mask from the source image: don't extend values to outside brain.
tiny = 5e-2; % From spm_vol_utils.c
Mask = true(size(y1));
Mask = Mask & (y1 >= (1-tiny) & y1 <= (SourceHead.dim(1)+tiny));
Mask = Mask & (y2 >= (1-tiny) & y2 <= (SourceHead.dim(2)+tiny));
Mask = Mask & (y3 >= (1-tiny) & y3 <= (SourceHead.dim(3)+tiny));
out(~Mask) = 0;
outmat = mat;
if ~nowrite
    OutHead             = SourceHead;
    OutHead.mat         = mat;
    OutHead.dim(1:3)    = dim;
    [p, n, e] = fileparts(SourceHead.fname);
    newname = sprintf('%s_%dx%dx%d%s',n,dim,e);
    OutHead.fname = [p filesep newname];
    spm_write_vol(OutHead,out);
end
end

