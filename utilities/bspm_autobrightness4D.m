function out = bspm_autobrightness4D(in, nowrite)
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



dat.min = min(in(in>0)); 
dat.max = max(in(in>0));
dat.dim = size(in); 
in = double(in)./255;
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
in(in>0) = scaledata(in(in>0), [dat.min dat.max]);
out = in;


if ~nowrite
    [fpath, fname, fext] = fileparts(h.fname); 
    h.fname = fullfile(fpath, ['e' fname fext]); 
    spm_write_vol(h, in);
end
