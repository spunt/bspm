function bspm_imcrop(allim, x, y, z)
% BSPM_IMSHOW Show transverse images in a Matlab figure
%
% USAGE: bspm_imshow(allim, x, y, z)
%

% ------------------ Copyright (C) 2014 ------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1, error('USAGE: bspm_imshow(allim)'); end
if nargin<2, x = []; end
if nargin<3, y = []; end
if nargin<4, z = []; end
if ischar(allim), allim = cellstr(allim); end
for i = 1:length(allim)
    
    im = allim{i}; 
    hdr = spm_vol(im); 
    im = spm_read_vols(hdr);
    im2 = imrotate(im, 90);
    imsize = size(im);
    prc = floor(prctile(1:imsize(1), [1 10 20 40 60 80 90 99]));
    figure('color', 'white');
    for ii = 1:length(prc)
        subplot(2,4,ii);
        imagesc(imrotate(squeeze(im(prc(ii),:,:)), 90)); colormap('gray'); title(num2str(prc(ii))); 
    end
    if isempty(x)
        x = input(sprintf('x (%d) : ', imsize(1)));
        if isempty(x), x = 1:imsize(1); end
    end
    if isempty(y)
        y = input(sprintf('y (%d) : ', imsize(2)));
        if isempty(y), y = 1:imsize(2); end
    end
    if isempty(z)
        z = input(sprintf('z (%d) : ', imsize(3)));
        if isempty(z), z = 1:imsize(3); end
    end
    im = im(x, y, z); 
    [p,n] = fileparts(allim{i});
    hdr.fname = fullfile(p, ['crop_' n '.nii']);
    hdr.dim = size(im);
    spm_write_vol(hdr, im);
    
end
   
 
 
 
 
