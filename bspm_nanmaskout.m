function bspm_nanmaskout(in, mask, outprefix)
% BSPM_NANMASKOUT
%
%   USAGE: bspm_nanmaskout(in, mask, outprefix)
%       
%       in  =  array of images OR wildcard pattern for finding them
%       rmtag = flag to delete old image, 0=No (default), 1=Yes 
%

% ----------------------- Copyright (C) 2014 -----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<2, error('USAGE: bspm_nanmaskout(in, mask, outprefix)'); end
if nargin<3, outprefix = 'nanout'; end
if ischar(in), in = cellstr(in); end
if iscell(mask), mask = char(mask); end
mask = bob_reslice(mask, in{1}, 1, 1); 
mask = mask > 0;
for i = 1:length(in)
    
    hdr = spm_vol(in{i});
    img = spm_read_vols(hdr);
    img(~mask) = NaN; 
    oldname = hdr.fname;
    [p, n, e] = fileparts(oldname);
    hdr.fname = [p filesep outprefix n e];
    spm_write_vol(hdr, img);
    
end

 
 
 
 
