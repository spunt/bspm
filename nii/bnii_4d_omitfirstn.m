function bnii_4d_omitfirstn(fname4d, nomit, outname)
% BNII_3DTO4D_FILENAME Find characters in common across multiple strings
%
%  USAGE: bnii_3dto4d(fname3d, varargin)
% __________________________________________________________________________
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-04-03
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
if nargin<2, mfile_showhelp; return; end
if nargin<3, outname = []; end
if ischar(fname4d), fname4d = cellstr(fname4d); end
for i = 1:length(fname4d)
    [img, hdr] = bnii_read(fname4d{i});
    img(:,:,:,1:nomit) = [];
    fp = fileparts(fname4d{i}); 
    if isempty(outname)
        outfn = fname4d{i}; 
    else
        outfn = fullfile(fp, outname); 
    end
    bnii_write(img, hdr, 'outname', outfn); 
end
fprintf('\n - Process Complete -\n');