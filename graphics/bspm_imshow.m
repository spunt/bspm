function bspm_imshow(im, slices)
% BSPM_IMSHOW Show transverse images in a Matlab figure
%
% USAGE: bspm_imshow(im, slices)
%

% ------------------ Copyright (C) 2014 ------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1, error('USAGE: bspm_imshow(im, slices)'); end
if nargin<2, slices = []; end
if iscell(im), im = im{1}; end
if ischar(im), im = spm_vol(im); im = spm_read_vols(im); end

imsize = size(im);
if isempty(slices), 
    msize = [3 5];
    nim = prod(msize);
    slices = 1:floor(imsize(3)/nim):imsize(3); 
else
    msize = [1 length(slices)];
    nim = prod(msize);
end

if length(slices)>nim
    extra = length(slices) - nim;
    rm1 = floor(extra/2);
    rm2 = ceil(extra/2);
    slices([1:rm1 end-rm2+1:end]) = [];
end

figure('Color','white');
for i = 1:nim
    subplot(msize(1),msize(2),i);
    imagesc(imrotate(im(:,:,slices(i)), 90));
    colormap('gray');
    title(sprintf('Slice %d',slices(i))); 
    axis off
end
    
 
 
 
 
