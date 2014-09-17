function XYZmm = bspm_vol_xyzmm(im)
% BSPM_VOL_XYZMM
%
%   USAGE: XYZmm = bspm_vol_xyzmm(im)
%
%   ARGUMENTS
%
%       im =  image
%       

% --------- Copyright (C) 2014 ---------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1, error('USAGE: XYZmm = bspm_vol_xyzmm(im)'); end
if iscell(im), im = char(im); end
hdr = spm_vol(im);
[R,C,P]  = ndgrid(1:hdr.dim(1),1:hdr.dim(2),1:hdr.dim(3));
RCP      = [R(:)';C(:)';P(:)'];
RCP(4,:) = 1;
XYZmm    = hdr.mat(1:3,:)*RCP;   


 
 
 
 
