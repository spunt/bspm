function bspm_checkreg(images, captions, xyz, title)
% A visual check of image registration quality.
%
% USAGE: bspm_checkreg(images, captions, xyz, title)
%
% FORMAT spm_check_registration
% FORMAT spm_check_registration(images, captions)
% Orthogonal views of one or more images are displayed. Clicking in
% any image moves the centre of the orthogonal views. Images are
% shown in orientations relative to that of the first selected image.
% The first specified image is shown at the top-left, and the last at
% the bottom right. The fastest increment is in the left-to-right
% direction (the same as you are reading this).
%__________________________________________________________________________
% Copyright (C) 1997-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_check_registration.m 4330 2011-05-23 18:04:16Z ged $

if ~nargin
    [images, sts] = spm_select([1 15],'image','Select images');
    if ~sts, return; end
end
if ~iscell(images) & strfind(images,'*'); images = files(images); end
if nargin < 2, captions = cellfun(@num2str, num2cell(1:length(images)), 'Unif', false); end
if ischar(images), images = cellstr(images); end
if iscell(images{1}), images = spm_vol(vertcat(images{:})); end
spm_figure('GetWin','Graphics');
spm_figure('Clear','Graphics');
spm_orthviews('Reset');
mn = length(images);
n  = round(mn^0.4);
m  = ceil(mn/n);
w  = 1/n;
h  = 1/m;
ds = (w+h)*0.02;
if nargin==4,
    if mn==2, 
        set(get(gca,'Title'), 'String', title, 'FontSize', 20, 'FontName', 'Arial', 'FontWeight', 'Bold', 'Position', [0.6 0.5]); 
    elseif mn>2,
        set(get(gca,'Title'), 'String', title, 'FontSize', 20, 'FontName', 'Arial', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center'); 
    end
    axis off; 
end
for ij=1:mn
    i = 1-h*(floor((ij-1)/n)+1);
    j = w*rem(ij-1,n);
    handle = spm_orthviews('Image', images{ij},...
        [j+ds/2 i+ds/2 w-ds h-ds]);
    if ij==1, spm_orthviews('Space'); end
    spm_orthviews('AddContext',handle);
    if ~isempty(captions)
        captions = cellstr(captions);
        mn = numel(captions);
        if ij <= mn
            spm_orthviews('Caption', ij, captions{ij}, 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
        end
    end
end

if exist('xyz', 'var')
    spm_orthviews('Reposition', xyz);
else
    spm_orthviews('Reposition',[0 0 0]);
end

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
