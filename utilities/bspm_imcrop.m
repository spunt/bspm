function bspm_imcrop(imfile, x, y, z)
% BSPM_IMSHOW Show transverse images in a Matlab figure
%
% USAGE: bspm_imcrop(imfile, x, y, z)
%

% ------------------ Copyright (C) 2014 ------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<1, mfile_showhelp; return; end
if nargin<2, x = []; end
if nargin<3, y = []; end
if nargin<4, z = []; end
if ischar(imfile), imfile = cellstr(imfile); end
for i = 1:length(imfile)
    
    im = imfile{i};
    [p,n] = fileparts(imfile{i});
    hdr = spm_vol(im); 
    im = spm_read_vols(hdr);
    imsize = size(im);
    prc = floor(prctile(1:imsize(3), [1 10 20 40 60 80 90 99]));
    figure('color', 'white');
    for ii = 1:length(prc)
        subplot(2,4,ii);
        imagesc(squeeze(imrotate(im(:,:,prc(ii)),180))); colormap('gray'); title(num2str(prc(ii))); 
    end
    if isempty(x)
        x = input(sprintf('KEEP IDX:\tx [1:%d] : ', imsize(1)));
        if isempty(x), x = 1:imsize(1); end
    end
    if isempty(y)
        y = input(sprintf('KEEP IDX:\ty [1:%d] : ', imsize(2)));
        if isempty(y), y = 1:imsize(2); end
    end
    if isempty(z)
        z = input(sprintf('KEEP IDX:\tz [1:%d] : ', imsize(3)));
        if isempty(z), z = 1:imsize(3); end
    end
    im = im(x, y, z); 
    hdr.fname = fullfile(p, ['crop_' n '.nii']);
    hdr.dim = size(im);
    spm_write_vol(hdr, im);
    
end
   
 
 
 
 
