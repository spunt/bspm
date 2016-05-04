function bspm_movie_compare(in1, in2, slice)
% BSPM_SHOW4DMOVIE
%
%   USAGE: bspm_show4Dmovie
%
if nargin==0, mfile_showhelp; return; end
if iscell(in1), in1 = char(in1); end
if iscell(in2), in2 = char(in2); end
fprintf('Reading data... ');
h1          = spm_vol(in1);
nvol        = length(h1);
dimvol      = h1(1).dim;
if nargin==2, slice       = round(dimvol(1)/2); end
[p1,n1]     = fileparts(in1);
[p2,n2]     = fileparts(in2); 
d1          = spm_read_vols(h1);
d1          = imrotate(squeeze(d1(slice, :, :, :)), 90);
h2          = spm_vol(in2); 
d2          = spm_read_vols(h2); 
d2          = imrotate(squeeze(d2(slice, :, :, :)), 90);
fprintf('DONE\n');
playagain = 1; 
while playagain==1
    figure('units', 'normal', 'position', [.05 .05 .75 .45], 'resize', 'off', 'color', 'white', 'menu', 'none', 'name', 'Timeseries Movie');
    
    subplot(1,2,1); h1 = imagesc(d1(:,:,1),[min(d1(:)) max(d1(:))]); colormap('gray'); axis off;
    title(n1, 'FontName', 'Arial', 'FontSize', 20, 'FontWeight', 'bold');
    hold on; 
    subplot(1,2,2); h2 = imagesc(d2(:,:,1),[min(d2(:)) max(d2(:))]); colormap('gray'); axis off;
    title(n2, 'FontName', 'Arial', 'FontSize', 20, 'FontWeight', 'bold');
    hold on; 
    for i=1:nvol
        set(h1,'CData',d1(:,:,i));
        set(h2,'CData',d2(:,:,i));
        pause(1/15);
    end;
    playagain = input('Play it again? [1=Yes, 2=No] ');
    close gcf
end
fprintf('Catch Ya Later!\n\n');
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

