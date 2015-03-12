function bspm_3dto4d(in, out)
% BSPM_3DTO4D
%
% USAGE: bspm_3dto4d(in, out)
%
% ARGUMENTS
%   in: array of 3D volumes to convert (cell or char)
%   out: name to give 4D volume [default = 4D.nii]
%
% Bob Spunt, November 18, 2012
% Based on function by Tor Wager

% ---------------- Copyright (C) 2014 ----------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, error('USAGE: bspm_3dto4d(in, out)'); end
if nargin < 2, out = '4D.nii'; end
if iscell(in), in = char(in);end
V    = spm_vol(in);
ind  = cat(1,V.n);
N    = cat(1,V.private);
mx   = -Inf;
mn   = Inf;
for i=1:numel(V),
    dat      = V(i).private.dat(:,:,:,ind(i,1),ind(i,2));
    dat      = dat(isfinite(dat));
    mx       = max(mx,max(dat(:)));
    mn       = min(mn,min(dat(:)));
end;
dat(dat==0) = NaN;
sf         = max(mx,-mn)/32767;
ni         = nifti;
ni.dat     = file_array(out,[V(1).dim numel(V)],'INT16-BE',0,sf,0);
ni.mat     = N(1).mat;
ni.mat0    = N(1).mat;
ni.descrip = '4D image';
create(ni);
for i=1:size(ni.dat,4),
    ni.dat(:,:,:,i) = N(i).dat(:,:,:,ind(i,1),ind(i,2));
    spm_get_space([ni.dat.fname ',' num2str(i)], V(i).mat);
end;
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
