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
if nargin < 1, error('USAGE: bspm_movie(in, [slice])'); end
if iscell(in), in = char(in); end
h       = spm_vol(in);
nvol    = length(h);
dimvol  = h(1).dim; 
if nargin < 2, slice = round(dimvol(1))/2; end
fprintf('\n | - Reading data for %d image volumes', nvol);
d       = spm_read_vols(h);
d       = imrotate(squeeze(d(slice, :, :, :)), 90);
fprintf(' - DONE\n');
playagain = 1; 
while playagain==1
    ddim = size(d);
    pos = [100 100 ddim([2 1])*10]; 
    figure('units', 'pixels', 'pos', pos, 'resize', 'off', 'color', 'black', 'menu', 'none', 'name', 'Timeseries Movie');    
    h1  = imagesc(d(:,:,1),[min(d(:)) max(d(:))]); colormap('gray'); axis off;
    f   = sprintf('%%0%dd', length(num2str(nvol))); 
    t   = title(sprintf(['Volume ' f ' of ' f], 1, nvol),  'FontName', 'Arial', 'Color', [1 1 1], 'FontSize', 20, 'FontWeight', 'bold');
    hold on; 
    for i=1:nvol
        set(t, 'String', sprintf(['Volume ' f ' of ' f], i, nvol)); 
        set(h1,'CData',d(:,:,i));
        pause(1/15);
    end;
    playagain = input('Play it again? [1=Yes, 2=No] ');
    close gcf
end
end
function vol = uigetvol(message, multitag)
% UIGETVOL Dialogue for selecting image volume file
%
%   USAGE: vol = uigetvol(message, multitag)
%       
%       message = to display to user
%       multitag = (default = 0) tag to allow selecting multiple images 
%
if nargin < 2, multitag = 0; end
if nargin < 1, message = 'Select Image File'; end
if ~multitag
    [imname, pname] = uigetfile({'*.img; *.nii', 'Image File'; '*.*', 'All Files (*.*)'}, message);
else
    [imname, pname] = uigetfile({'*.img; *.nii', 'Image File'; '*.*', 'All Files (*.*)'}, message, 'MultiSelect', 'on');
end
if isequal(imname,0) || isequal(pname,0)
    vol = [];
else
    vol = fullfile(pname, strcat(imname));
end
end

 
 
 
