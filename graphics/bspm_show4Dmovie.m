function bspm_show4Dmovie
% BSPM_SHOW4DMOVIE
%
%   USAGE: bspm_show4Dmovie
%
in1 = uigetvol('Select uncorrected 4D image file');
fprintf('\nUncorrected file... %s\n', in1);
in2 = uigetvol('Select corrected 4D image file');
fprintf('Corrected file..... %s\n', in2);
[rpfile, pathname] = uigetfile({'rp*txt', 'Realignment Parameter File'}, 'Select realignment parameter file'); 
rpfile = fullfile(pathname, rpfile);
fprintf('Realignment Parameter file..... %s\n', rpfile);
h1       = spm_vol(in1);
nvol    = length(h1);
dimvol  = h1(1).dim; 
slice = round(dimvol(1)/2);
fprintf('Reading data... ');
d1       = spm_read_vols(h1);
d1       = imrotate(squeeze(d1(slice, :, :, :)), 90);
h2 = spm_vol(in2); 
d2 = spm_read_vols(h2); 
d2 = imrotate(squeeze(d2(slice, :, :, :)), 90);
rp = load(rpfile);
rp(:,4:6) = rp(:,4:6)*57.3;
fprintf('DONE\n');
playagain = 1; 
while playagain==1
    figure('units', 'normal', 'position', [.05 .05 .75 .9], 'resize', 'off', 'color', 'white', 'menu', 'none', 'name', 'Timeseries Movie');
    
    subplot(3,2,1); h1 = imagesc(d1(:,:,1),[min(d1(:)) max(d1(:))]); colormap('gray'); axis off;
    title('Before Motion Correction', 'FontName', 'Arial', 'FontSize', 20, 'FontWeight', 'bold');
    hold on; 
    
    subplot(3,2,2); h2 = imagesc(d2(:,:,1),[min(d2(:)) max(d2(:))]); colormap('gray'); axis off;
    title('After Motion Correction', 'FontName', 'Arial', 'FontSize', 20, 'FontWeight', 'bold');
    hold on; 
    
    subplot(3,2,3:4); h3 = plot(rp(:,1:3));
    ylabel('Translation (mm)', 'FontName', 'Arial', 'FontSize', 20);
    legendh1 = legend({'x' 'y' 'z'});
    set(legendh1, 'Orientation', 'vertical', 'Location', 'best', 'FontSize', 13);
    ax = axis; 
    l(1) = line([1 1], ax(3:4), 'color', [.4 .4 .4], 'linewidth', 6);
    hold on; 
    
    subplot(3,2,5:6); h4 = plot(rp(:,4:6));
    xlabel('Image Number', 'FontName', 'Arial', 'FontSize', 20); 
    ylabel('Rotation (degrees)', 'FontName', 'Arial', 'FontSize', 20);
    legendh2 = legend({'pitch' 'roll' 'yaw'});
    set(legendh2, 'Orientation', 'vertical', 'Location', 'best', 'FontSize', 13);
    ax = axis; 
    l(2) = line([1 1], ax(3:4), 'color', [.1 .1 .1], 'linewidth', 6);
    hold on;  
    
    for i=1:nvol
        set(h1,'CData',d1(:,:,i));
        set(h2,'CData',d2(:,:,i));
        set(l, 'XData', [i i]);
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

