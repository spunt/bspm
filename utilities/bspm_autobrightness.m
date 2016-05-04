function out = bspm_autobrightness(in, nowrite)
% AUTOBRIGHTNESS Auto-adjust brightness of image volume
%
%  USAGE: out = bspm_autobrightness(in, *nowrite)	*optional input
% __________________________________________________________________________
%  INPUTS
%	in:  image volume filename 
%   nowrite: option to not write outfile
%

% ---------------------- Copyright (C) 2014 Bob Spunt ----------------------
%	Created:  2014-10-04
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 1, mfile_showhelp; return; end
if nargin < 2, nowrite = 0; end
lim = 0.5;
if iscell(in), in = char(in); end
h       = spm_vol(in);
in      = spm_read_vols(h);
iszero  = in==0;
in(iszero) = NaN;
dat.min = nanmin(in(:)); 
dat.max = nanmax(in(:));
dat.dim = size(in);
% in = double(in)./255;
in = cat(4,in,in,in); 
for p = 1:3
    for s = 1:dat.dim(p)
        if p==1
            in(s,:,:,p) = in(s,:,:,p) + (lim - nanmean(nanmean(in(s,:,:,p))))*(1 - in(s,:,:,p));
        elseif p==2
            in(:,s,:,p) = in(:,s,:,p) + (lim - nanmean(nanmean(in(:,s,:,p))))*(1 - in(:,s,:,p));
        else
            in(:,:,s,p) = in(:,:,s,p) + (lim - nanmean(nanmean(in(:,:,s,p))))*(1 - in(:,:,s,p));
        end
    end
end
in = nanmean(in, 4); 
in(iszero) = 0; 
in(~iszero) = scaledata(in(~iszero), [dat.min dat.max]);
if nargout>0, out = in; end
if ~nowrite
    [fpath, fname, fext] = fileparts(h.fname); 
    h.fname = fullfile(fpath, ['e' fname fext]); 
    spm_write_vol(h, in);
end
