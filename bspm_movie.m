function bspm_movie(in, slice)
% BSPM_MOVIE
%
% USAGE: bspm_movie(in, slice)
%

% ------ Copyright (C) 2014 ------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, error('USAGE: bspm_movie(in, slice)'); end
if iscell(in), in = char(in); end
h       = spm_vol(in);
nvol    = length(h);
dimvol  = h(1).dim; 
if nargin < 2, slice = round(dimvol(1))/2; end
fprintf('\nReading data for %d image volumes... ', nvol);
d       = spm_read_vols(h);
d       = imrotate(squeeze(d(slice, :, :, :)), 90);
fprintf('DONE\n');
imagesc(d(:,:,1)); colormap('gray'); 
axis tight; set(gca,'nextplot','replacechildren','visible','off');
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,nvol) = 0;
for i = 2:nvol
  imagesc(d(:,:,i)); colormap('gray');
  f = getframe;
  im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
end
 
 
 
 
