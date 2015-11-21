function nii_viewer(fname, overlayName)
% Basic tool to visualize NIfTI images.
% 
%  NII_VIEWER('/data/subj2/fileName.nii.gz')
%  NII_VIEWER('background.nii', 'overlay.nii')
% 
% If no input is provided, the viewer will ask for background NIfTI.
% 
% Here are some features and usage.
% 
% The basic use is to open a NIfTI file to view. When a NIfTI (background) is
% open, the display always uses the image plane close to xyz axes (voxel space)
% even for tilted acquisition. The possible confusion comes if the acquisition
% was tilted with a large angle, and then the orientation labeling will be less
% accurate. The benefit is that no interpolation is needed for the background
% image. The display is always in correct scale at three axes even with
% non-isotropic voxels. The displayed IJK always correspond to left -> right,
% posterior -> anterior and inferior -> superior directions, although the NIfTI
% data may not be saved in this order or along these directions. The I index is
% increasing from left to right even when the display is flipped as radiological
% convention (right on left side).
% 
% Navigation in 4D can be done by mouse click, dialing XYZ numbers, or using
% keys (arrow keys and [] for 3D, and <> for volume). 
% 
% After the viewer is open, dragging and dropping a NIfTI file into the viewer
% will open it as background.
% 
% By default, the viewer shows full view of the image data. The zoom-in always
% applies to three views together, and enlarges around the location of the
% crosshair. To center at a different location, set the crosshair to the
% interested location, and apply zoom again either by View -> Zoom in, or
% pressing Ctrl (Cmd) and +/-.
% 
% Overlays are always transformed into the space of background image, so
% interpolation (nearest/linear/cubic/spline) is usually involved. The overlay
% makes sense only when it has the same coordinate system as the background
% image, while the resolution and dimension can be different. The viewer tries
% to match any of sform and qform between the images. If there is no match, a
% warning message will show up.
% 
% A special overlay feature "Add aligned overlay" can be used to check the
% effect of registration. It will ask for a NIfTI and registration matrix which
% aligns the NIfTI to the background image. Here is an example to check FSL
% alignment. From a .feat/reg folder, Open "highres" as background image. "Add
% overlay" for "example_func". If there is head movement between the highres and
% the functional scans, the overlap will be off. Now "Add aligned overlay" for
% "example_func", and use "example_func2highres.mat" as the matrix. The two
% dataset should overlap well if the alignment matrix is accurate.
% 
% When the mouse is moved onto a voxel, the voxel indices, corresponding x/y/z
% coordinates and voxel value will show on the panel. If there is an overlay,
% the overlay voxel value will also show up, unless its display is turned off.
% When the mouse is outside an image, the information for the voxel at crosshair
% will be displayed. The display format is as following:
%  (i,j,k)=(x,y,z): val_background val_overlay1 val_overlay2 ...
%
% Note that although the x/y/z coordinates are shared by background and overlay
% images, indices are always for background image. 
% 
% Image display can be smoothed for background and overlays in 3D (default is
% off). The smooth is slow when the image dimension is large, even when the
% current implementation of smooth does not consider voxel size.
% 
% Background image and overlays are listed in a popup menu at the left side of
% the panel. All parameters at the right side of the list are for the selected
% file. This feature is indicated by a frame grouping these parameters. Each
% NIfTI file has its own set of parameters (display min and max value, LUT,
% alpha, whether to smooth, interpolation method, and number of volumes) to
% control its display. Moving the mouse onto a parameter will show its meaning.
% 
% If more than one overlays are added, the last added overlay will be on the top
% of display, while it is at the bottom of the file list. The overlay order can
% be changed easily from Overlay -> Move overly ...
% 
% Each NIfTI display can be turned on/off by clicking the small checkbox next to
% the file list. This provides a way to turn on/off an overlay back and forth to
% view the overlap. Most operations are applied to the selected NIfTI in the
% list, such as Show NIfTI hdr/ext under Window menu, Move/Close overlay under
% Overlay menu, and most operations under File menu.
% 
% A NIfTI mask can be applied to the selected image. Ideally, the mask should be
% binary, and only the non-zero part of the mask will be displayed. In case
% non-binary mask is detected, a threshold to binarize will be asked. If the
% effect is not satisfied with a threshold, one can apply the same mask with a
% different threshold without re-loading anything. The way to remove a mask is
% to Close, then Add overlay again.
% 
% For multi-volume data, one can change the Volume Number (the parameter at
% rightmost of the panel) to check the head motion. Click in the number dialer,
% and hold up or down arrow key, or press < or > key, to simulate movie play. It
% is better to open the 4D data as background, since it will be slower to map it
% to background image.
% 
% Several simple LUT options are implemented. The color coding can be shown by
% color bar (View -> Show colorbar). The last two LUT options are special. The
% "two-sided" allows to show both positive and negative data in one view. For
% example, if the display range is 3 to 10 for a t-map, positive T above 3 will
% be coded as red-yellow, and T below -3 will be coded as blue-green. This means
% the absolute display range values are used.
% 
% The other special LUT is "lines". This is for diffusion vector display. Under
% this LUT, all other parameters for the display are ignored. Also this is
% normally overlaid onto FA map from the same acquisition, so it requires the
% background with the same resolution and dimension. The color of the "lines" is
% the max color of previous LUT. For example, if one likes to show blue vector
% lines, choose LUT blue first, then change it to "lines".
% 
% The figure can be copied into clipboard (not available for Linux) or saved as
% variety of image format. For high quality picture, one can increase the output
% resolution from Help -> Preferences -> Resolution. Higher resolution will take
% longer time to copy or save, and result in large file. If needed, one can
% change to white background for picture output. With white background, the
% threshold for the background image needs to be carefully selected to avoid
% black/white strips. White background is intended only for picture output.
% 
% The selected NIfTI can also be saved into different format from File -> Save
% NIfTI as. For example, a file can be saved as a different resolution. With a
% transformation matrix, a file can also be saved into a different template. The
% latter is needed for FSLview since it won't allow overlay with different
% resolution or dimension.
% 
% See also NII_TOOL, DICM2NII, NII_XFORM

% By Xiangrui Li (xiangrui.li@gmail.com)
% History(yymmdd):
% 151021 Almost ready to publish.
% 151104 Include drag&drop by Maarten van der Seijs.
% 151105 Bug fix for Show NIfTI hdr/ext.
% 151106 Use hs.q.interp to avoid unnecessary interp for overlay;
%        Use check mark for colorbar/crosshair/background menu items.
% 151111 Take care of slope/inter of img; mask applies to DTI lines.
% 151114 Avoid see-thu for white background.
% 151119 Make smooth faster by using only 3 slices; Allow flip L/R.
% 151120 Implement key navigation and zoom in/out.
% End of history. Don't edit this line!

if nargin>1 && ischar(overlayName) && strcmp(overlayName, 'gui_callback')
    nii_viewer_cb(fname);
    return;
end

% load preference file if exist
if ispc, home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else home = getenv('HOME');
end
fnameH = fullfile(home, 'dicm2nii_gui_para.mat');
pref_file = fnameH;
if ~exist(pref_file, 'file')
    pref_file = [fileparts(which('dicm2nii')) '/dicm2nii_gui_para.mat'];
    if ~exist(pref_file, 'file')
        fid = fopen(pref_file, 'w'); % check permission
        if fid<1
            pref_file = fnameH;
        else
            fclose(fid); delete(pref_file);
        end
    end
end
pf = struct('openPath', pwd, 'addPath', pwd, 'interp', 'linear', ...
    'extraV', NaN, 'dpi', '0', 'pref_file', pref_file);
if exist(pref_file, 'file')
    para = load(pref_file); para = para.para;
    try pf = para.nii_viewer; catch, end
else
    para.nii_viewer = pf;
    try save(pref_file, 'para'), catch, end
end

if nargin<1
    [fname, pName] = uigetfile([pf.openPath '/*.nii; *.hdr; *.nii.gz;*.hdr.gz'], ...
        'Select NIfTI to view');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    para.nii_viewer.openPath = pName; %#ok
    try save(pref_file, 'para'), catch, end
end

[hs.q, hs.form_code, rg, dim, hs.pixdim] = read_nii(fname);
nVol = size(hs.q.nii.img, 4);

res = get(0, 'ScreenSize'); res = res(3:4)-res(1:2)+1; % screen resolution
mm = dim .* hs.pixdim; % FOV
siz = [sum(mm(1:2)) sum(mm(2:3))]; % image area width/height
x0 = mm(1) / siz(1); % normalized width of left images
y0 = mm(2) / siz(2); % normalized height of bottom image (tra)
z0 = mm(3) / siz(2); % normalized height of top images

scr = [res(1)-40 res(2)-156]; % room for system bar and figure menu/title/panel
siz = siz * min(scr ./ siz); % almost max size
if siz(2)>800, siz = siz*0.8; end % make it smaller on big screen

pos = round((res-siz)/2);
if pos(1)+siz(1) > res(1), pos(1) = 1; end
if pos(2)+siz(2) > res(2)-130, pos(2) = min(pos(2), 50); end

p.show = true;
p.lb = rg(1); 
p.ub = rg(2);
p.lut = 1;
p.alpha = 1;
p.smooth = false;
p.interp = 1;
p.volume = 1;

hs.dim = single(dim); % single may save a lot memory for meshgrid
hs.siz = siz; % image panel size
hs.gap = min(hs.pixdim) ./ hs.pixdim * 3; % units of smallest pixdim

fname = hs.q.nii.hdr.file_name;
[pName, niiName, ext] = fileparts(fname);
if strcmpi(ext, '.gz'), [~, niiName] = fileparts(niiName); end

set(0, 'ShowHiddenHandles', 'on');
n = 'ni'*256.^(1:2)'; % start with a big number for figure
while 1
    n = n+1;
    fh = figure(n);
    if ~strcmp(get(fh,'Tag'), 'nii_viewer'), break; end
end
set(0, 'ShowHiddenHandles', 'off');
set(fh, 'Toolbar', 'none', 'Menubar', 'none', 'UserData', {fname}, 'Render', 'opengl', ...
    'Name', ['nii_viewer - ' fname], 'NumberTitle', 'off', 'Tag', 'nii_viewer', ...
    'Position', [pos siz+[0 64]], 'DockControls', 'off');
hs.fig = fh;
cb = @(cmd) ['nii_viewer(''' cmd ''', ''gui_callback'');']; % callback shortcut

%% menus
h = uimenu(fh, 'Label', '&File');
uimenu(h, 'Label', 'Open', 'Accelerator', 'O', 'UserData', pName, 'Callback', cb('open'));
uimenu(h, 'Label', 'Apply mask', 'Callback', cb('mask'));
h_savefig = uimenu(h, 'Label', 'Save figure as');
h_saveas = uimenu(h, 'Label', 'Save NIfTI as');
uimenu(h, 'Label', 'Close window', 'Accelerator', 'W', 'Callback', 'close gcf');

uimenu(h_saveas, 'Label', 'SPM 3D NIfTI (one file/pair per volume)', 'Callback', cb('saveas'));
uimenu(h_saveas, 'Label', 'NIfTI standard RGB (for AFNI, later mricron, SPM)', ...
    'Callback', cb('saveas'), 'Separator', 'on');
uimenu(h_saveas, 'Label', 'FSL style RGB (RGB saved in dim 4)', 'Callback', cb('saveas'));
uimenu(h_saveas, 'Label', 'Old mricron style RGB (RGB saved in dim 3)', 'Callback', cb('saveas'));
uimenu(h_saveas, 'Label', 'file with a new resolution', 'Callback', cb('saveas'), 'Separator', 'on');
uimenu(h_saveas, 'Label', 'file matching background', 'Callback', cb('saveas'));
uimenu(h_saveas, 'Label', 'file in aligned template space', 'Callback', cb('saveas'));

fmt = {'pdf' 'eps' 'png' 'jpg' 'tif' 'bmp'};
if ispc, fmt = [fmt 'emf']; end
for i = 1:numel(fmt)
    uimenu(h_savefig, 'Label', fmt{i}, 'Callback', cb('save'));
end

if ispc || ismac
    h = uimenu(fh, 'Label', '&Edit');
    uimenu(h, 'Label', 'Copy figure', 'Callback', cb('copy'));
end

h_over = uimenu(fh, 'Label', '&Overlay');
hs.add = uimenu(h_over, 'Label', 'Add overlay', 'Accelerator', 'A', ...
    'UserData', pf.addPath, 'Callback', cb('add'));
uimenu(h_over, 'Label', 'Add aligned overlay', 'Callback', cb('add'));

h = uimenu(h_over, 'Label', 'Move overlay', 'Enable', 'off');
uimenu(h, 'Label', 'to top',         'Callback', cb('stack'));
uimenu(h, 'Label', 'to bottom',      'Callback', cb('stack'));
uimenu(h, 'Label', 'one level up',   'Callback', cb('stack'));
uimenu(h, 'Label', 'one level down', 'Callback', cb('stack'));
hs.overlay = h;

hs.overlay(2) = uimenu(h_over, 'Label', 'Remove overlay', 'Accelerator', 'R', ...
    'Callback',  cb('close'), 'Enable', 'off');
hs.overlay(3) = uimenu(h_over, 'Label', 'Remove overlays', 'Accelerator', 'Q', ...
    'Callback', cb('closeAll'), 'Enable', 'off');

h_view = uimenu(fh, 'Label', '&View');
h = uimenu(h_view, 'Label', 'Zoom in by');
for i = [1 1.2 1.5 2 3 4 5 8 10 20]
    uimenu(h, 'Label', num2str(i), 'Callback',  cb('zoom'));
end
uimenu(h_view, 'Label', 'White background', 'Callback',  cb('background'));
uimenu(h_view, 'Label', 'Right on left side', 'Callback',  cb('flipLR'));
uimenu(h_view, 'Label', 'Show colorbar', 'Callback',  cb('colorbar'));
uimenu(h_view, 'Label', 'Show crosshair', 'Separator', 'on', ...
    'Checked', 'on', 'Callback',  cb('cross'));
uimenu(h_view, 'Label', 'Crosshair color', 'Callback',  cb('color'));
h = uimenu(h_view, 'Label', 'Crosshair gap');
for i = [1 2 3 4 5 6 8 10 20 40]
    str = num2str(i); if i==6, str = [str ' (default)']; end %#ok
    uimenu(h, 'Label', str, 'Callback', cb('gap'));
end
h = uimenu(h_view, 'Label', 'Crosshair thickness');
uimenu(h, 'Label', '0.5 (default)', 'Callback', cb('thickness'));
for i = [0.75 1 2 4 8]
    uimenu(h, 'Label', num2str(i), 'Callback', cb('thickness'));
end

h = uimenu(fh, 'Label', '&Window');
uimenu(h, 'Label', 'Show NIfTI hdr', 'Callback', cb('hdr'));
uimenu(h, 'Label', 'Show NIfTI ext', 'Callback', cb('ext'));
uimenu(h, 'Label', 'DICOM to NIfTI converter', 'Callback', 'dicm2nii', 'Separator', 'on');

h = uimenu(fh, 'Label', '&Help');
hs.pref = uimenu(h, 'Label', 'Preferences', 'Callback', cb('pref'), 'UserData', pf);
uimenu(h, 'Label', 'Key shortcut', 'Callback', cb('keyHelp'));
uimenu(h, 'Label', 'Show help text', 'Callback', 'doc nii_viewer');
uimenu(h, 'Label', 'About', 'Callback', cb('about'));

%% Three views: sag, cor, tra
% this panel only makes resize easy: subplot relative to the panel
ph = uipanel(fh, 'Units', 'pixels', 'Position', [1 1 siz], ...
    'BorderType', 'none', 'BackgroundColor', 'k');
hs.im_panel = ph;
c = round(hs.q.R \ [0 0 0 1]'); c = c(1:3)' + 1;
ind = c<=1 | c>=dim;
c(ind) = round(dim(ind)/2);

pos = [x0 y0 1-x0 z0;  0 y0 x0 z0;  0 0 x0 y0];
for i = 1:3
    j = 1:3; j(j==i) = [];
    hs.ax(i) = subplot('Position', pos(i,:), 'Parent', ph);
    hs.q.hsI(i) = image(zeros(dim(j([2 1])))); hold(hs.ax(i), 'on');
    
    x = c(j(1))+[-1 1 0 0]*hs.gap(j(1)); u = [-dim(j(1))-1 dim(j(1))+1 0 0];
    y = c(j(2))+[0 0 -1 1]*hs.gap(j(2)); v = [0 0 -dim(j(2))-1 dim(j(2))+1];
    hs.cross(i) = quiver(x, y, u, v, 'ShowArrowHead', 'off', 'AutoScale', 'off');
end

labls='ASLSLP';
pos = [0.95 0.5; 0.47 0.96;  0 0.5; 0.47 0.96; 0 0.5; 0.47 0.05]; 
for i = 1:numel(labls)
    hs.ras(i) = text(pos(i,1), pos(i,2), labls(i), 'Units', 'normalized', ...
        'Parent', hs.ax(ceil(i/2)), 'FontSize', 12, 'FontWeight', 'bold');
end

% early matlab's colormap works only for axis, so ax(4) is needed.
hs.ax(4) = subplot('Position', [x0 0 1-x0 y0], 'Parent', ph, 'Ylim', [0 1]);
colorbar('Units', 'Normalized', 'Position', [x0+0.2 y0*0.15 0.03 y0*0.7]);
hs.colorbar = findobj(fh, 'tag', 'Colorbar'); % trick for early matlab
set(hs.colorbar, 'Visible', 'off', 'UIContextMenu', '', 'EdgeColor', [1 1 1]);
% hs.colorbar = colorbar(hs.ax(4), 'YTicks', [0 0.5 1], 'Color', [1 1 1], ...
%     'Location', 'west', 'PickableParts', 'none', 'Visible', 'off');

% image makes YDir reversed. Turn off ax and ticks
set(hs.ax, 'YDir', 'normal', 'Visible', 'off');
set([hs.ras hs.cross], 'Color', 'b', 'UIContextMenu', ''); %, 'PickableParts', 'none');

%% control panel
pos = get(fh, 'Position'); pos = [1 pos(4)-64 pos(3) 64];
ph = uipanel(fh, 'Units', 'pixels', 'Position', pos, 'BorderType', 'none');
hs.panel = ph;
clr = get(ph, 'BackgroundColor');

% IJK java spinners
labls = 'XYZ';
cmd = {'spin_x' 'spin_y' 'spin_z'};
str = {'Left to Right' 'Posterior to Anterior' 'Inferior to Superior'};
pos = [38 44 22]; posTxt = [38 12 20];
for i = 1:numel(labls)
    loc = [(i-1)*64+14 pos];
    txt = sprintf('%s, 1:%g', str{i}, dim(i));
    hs.spinner(i) = java_spinner(loc, [c(i) 1 dim(i) 1], ph, cb(cmd{i}), '#', txt);
    uicontrol(ph, 'Style', 'text', 'String', labls(i), 'BackGroundColor', clr, ...
        'Position', [loc(1)-12 posTxt], 'TooltipString', txt);
end

hs.xyz = uicontrol(ph, 'Style', 'text', 'Position', [190 38 siz(1)-190 20], ...
    'BackGroundColor', clr, 'UserData', false);

% Controls for each file
uipanel(ph, 'Units', 'normalized', 'Position', [0 2/64 1 0.53], ...
    'BorderType', 'etchedin', 'BorderWidth', 2);
hs.files = uicontrol(ph, 'Style', 'popup', 'BackgroundColor', 'w', ...
    'String', {niiName}, 'Position', [8 8 82 22], 'Value', 1, 'Callback', cb('files'), ...
    'TooltipString', 'All parameters at right are for the selected NIfTI');
hs.show = uicontrol(ph, 'Style', 'checkbox', 'Position', [90 8 22 20], 'Value', 1,  ...
    'Callback', cb('show'), 'BackGroundColor', clr, ...
    'TooltipString', 'Turn on/off selected image');

hs.lb = java_spinner([108 8 52 22], [rg(1) -inf inf 10], ph, cb('lb'), '#.##', 'min value (threshold)');
hs.ub = java_spinner([160 8 52 22], [rg(2) -inf inf 10], ph, cb('ub'), '#.##', 'max value (clipped)');
hs.lut = uicontrol(ph, 'Style', 'popup', 'Position', [214 8 74 22], ...
    'String', {'grayscale' 'red' 'green' 'blue' 'violet' 'yellow' 'cyan' ...
    'red-yellow' 'blue-green' 'two-sided' 'lines'}, ...
    'BackgroundColor', 'w', 'Callback', cb('lut'), 'Value', 1, ...
    'TooltipString', 'Lookup table options for non-RGB data');

hs.alpha = java_spinner([288 8 44 22], [1 0 1 0.1], ph, cb('alpha'), '#.#', ...
    'Alpha: 0 transparent, 1 opaque');

hs.smooth = uicontrol(ph, 'Style', 'checkbox', 'value', p.smooth, ...
    'FontSize', 8, 'Position', [332 8 60 22], 'String', 'smooth', ...
    'BackGroundColor', clr, 'Callback', cb('smooth'), ...
    'TooltipString', 'Smooth image in 3D');
hs.interp = uicontrol(ph, 'Style', 'popup', 'String', {'nearest' 'linear' 'cubic' 'spline'}, ...
    'Position', [392 8 68 22], 'value', p.interp, 'BackGroundColor', 'w', ...
    'Callback', cb('interp'), 'Enable', 'off', ... 
    'TooltipString', 'Interpolation method for overlay');
hs.volume = java_spinner([458 8 44 22], [1 1 nVol 1], ph, cb('volume'), '#', ...
    ['Volume number, 1:' num2str(nVol)]);
set(hs.volume, 'Enable', nVol>1);

%% finish gui
set(fh, 'ResizeFcn', cb('resize'), ...
    'WindowButtonMotionFcn', cb('mousemove'), ...    
    'WindowButtonDownFcn', cb('mousedown'), ...
    'WindowButtonUpFcn', cb('mouseup'), ...
    'WindowKeyPressFcn', @KeyPressFcn, ...
    'Interruptible', 'off', 'BusyAction', 'cancel', ...
    'PaperPositionMode', 'auto', 'HandleVisibility', 'Callback');
set(hs.files, 'UserData', p); % store file specific parameters
guidata(fh, hs); % store handles

% java dnd by Maarten van der Seijs: matlabcentral/fileexchange/53511 
try
    warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    jFrame = get(fh, 'JavaFrame');
    jObj = jFrame.getAxisComponent;
    dndcontrol.initJava();
    dndcontrol(jObj, @javaDropFcn);
catch me
    fprintf(2, '%s\n', me.message);
end
drawnow; % drawnow needed in case of multiple figures

if nargin>1, nii_viewer_cb('add', overlayName); end
nii_viewer_cb('update');
    
%% callbacks
function nii_viewer_cb(cmd, varargin)
h = gcbo;
try
    hs = guidata(h);
    fh = hs.fig;
catch % java objects won't work with guidata
    set(0, 'ShowHiddenHandles', 'on');
    fh = findobj('Type', 'figure', 'Tag', 'nii_viewer');
    set(0, 'ShowHiddenHandles', 'off');
    fh = fh(1); % top one if multiple
    hs = guidata(fh);
end
% if ~strcmp(cmd, 'mousemove'), disp(cmd); end

switch cmd
    case 'update'
        for ix = 1:3, set_cdata(ix, fh); end
    case 'spin_x'
        x = hs.spinner(1).getValue;
        set_cdata(1, fh, x);
        set(hs.cross([2 3]), 'XData', x+[-1 1 0 0]*hs.gap(1));
        set_xyz(hs);
    case 'spin_y'
        y = hs.spinner(2).getValue;
        set_cdata(2, fh, y);
        set(hs.cross(1), 'XData', y+[-1 1 0 0]*hs.gap(2));
        set(hs.cross(3), 'YData', y+[0 0 -1 1]*hs.gap(2));
        set_xyz(hs);
    case 'spin_z'
        z = hs.spinner(3).getValue;
        set_cdata(3, fh, z);
        set(hs.cross(1:2), 'YData', z+[0 0 -1 1]*hs.gap(3));
        set_xyz(hs);
    case {'lb' 'ub' 'lut' 'alpha' 'smooth' 'interp' 'volume'};
        i = get(hs.files, 'Value');
        p = get(hs.files, 'UserData');
        val = get(hs.(cmd), 'Value');
        
        if strcmp(cmd, 'lut') && val==11 % error check for vector lines
            set(hs.lut, 'UserData', p(i).lut); % remember old lut
            err = true;
            if i==1 % need to copy hsI from background when lut changes
                errordlg('"lines" applies only to overlay');
            elseif hs.q(i).interp % don't want interpolate normalized vector
                errordlg('"lines" needs identical-dim background image');
            elseif size(hs.q(i).nii.img,4)~=3
                errordlg('Not valid vector data: dim4 is not 3');
            else
                a = sum(hs.q(i).nii.img.^2, 4); a = a(a(:)>0);
                if any(abs(a-1)>0.01)
                    errordlg('Not valid vector data: squared sum is not 1');
                else err = false; % passed all checks
                end
            end
            if err, set(hs.lut, 'Value', p(i).lut); return; end
        end
        
        % hs.files.UserData(i).(cmd) = val; % for later matlab
        p(i).(cmd) = val;
        set(hs.files, 'UserData', p);
        if any(strcmp(cmd, {'lut' 'lb' 'ub'})) && ...
                strcmpi(get(hs.colorbar, 'Visible'), 'on')
            nii_viewer_cb('colorbar', 1);
        end
        nii_viewer_cb('update');
        if strcmp(cmd, 'volume'), set_xyz(hs); end
    case 'resize'
        if isempty(hs), return; end
        cb = get(fh, 'ResizeFcn');
        set(fh, 'ResizeFcn', ''); drawnow; % avoid weird effect
        clnObj = onCleanup(@() set(fh, 'ResizeFcn', cb)); % restore func
        
        posP = get(hs.panel, 'Position'); % get old height in pixels
        posF = get(fh, 'Position'); % asked position by user
        posI = get(hs.im_panel, 'Position'); % old size
        
        res = get(0, 'MonitorPositions');
        if size(res,1)<2 % single monitor
            res = res(1,3:4) - res(1,1:2) - min(res(:));
        elseif any(res(1,1:2) > res(1,1:2)) % dual
            res = res(1,1:2) + res(2,3:4) - min(res(:));
        else
            res = res(2,1:2) + res(1,3:4) - min(res(:));
        end
        
        oldF = round([posI(3) posI(4)+posP(4)]); % old fig size
        if isequal(oldF, posF(3:4)), return; end
        if all(posF(3:4) >= oldF) % enlarge
            a = max([posF(3) posF(4)-posP(4)] ./ hs.siz) * hs.siz;
            a(1) = min(a(1), res(1)-30); % leave space for MAC dock etc
            a(2) = min(a(2), res(2)-92-posP(4)); % leave space for title bar etc
            a = min(a ./ hs.siz) * hs.siz;
        elseif all(posF(3:4) <= oldF) % shrink
            a = min([posF(3) posF(4)-posP(4)] ./ hs.siz) * hs.siz;
        else % one side enlarge, another side shrink
            a = posI(3:4);
        end
        d = posF(1)+a(1)-res(1);
        if d>0, posF(1) = posF(1) - d; end
        d = posF(2)+a(2)+posP(4)+92-res(2);
        if d>0, posF(2) = posF(2) - d; end
        posF(1) = max(posF(1), 10);
        posF(2) = max(posF(2), 50);
        posF(3:4) = [a(1) a(2)+posP(4)]; % final figure size
        set(fh, 'Position', posF); % done for fig
        
        posP(2) = posF(4)-posP(4)+1; 
        posP(3) = posF(3);
        set(hs.panel, 'Position', posP); % done for control panel
        set(hs.im_panel, 'Position', [1 1 a]); % done for image panel
        
        pos = get(hs.xyz, 'Position');
        pos(3) = max(1, posP(3)-pos(1)-1);
        set(hs.xyz, 'Position', pos);
    case 'show' % turn on/off NIfTI
        i = get(hs.files, 'Value');
        p = get(hs.files, 'UserData');
        p(i).show = get(hs.show, 'Value');
        set(hs.files, 'UserData', p);
        states = {'off' 'on'};
        set(hs.q(i).hsI, 'Visible', states{p(i).show +1});
        if p(i).show, nii_viewer_cb('update'); end % skipped in set_cdata
    case 'files'
        i = get(hs.files, 'Value');
        p = get(hs.files, 'UserData');
        nam = {'lb' 'ub' 'alpha' 'volume' 'show' 'lut' 'smooth' 'interp' };
        cb = cell(1,4); % first 4 in nam are java objects
        for j = 1:4 % disable spinner callback, not needed for others
            cb{j} = get(hs.(nam{j}), 'StateChangedCallback');
            set(hs.(nam{j}), 'StateChangedCallback', '');
        end
       
        for j = 1:numel(nam)
            set(hs.(nam{j}), 'Value', p(i).(nam{j}));
        end
        if i==1, set(hs.interp, 'Enable', 'off');
        else set(hs.interp, 'Enable', 'on');
        end
        nVol = size(hs.q(i).nii.img, 4);
        set(hs.volume, 'Enable', nVol>1, ...
            'ToolTipText', ['Volume number, 1:' num2str(nVol)]);
        set(hs.volume.Model, 'Maximum', nVol);
        a = max(str2double(sprintf('%.1g', abs(p(i).lb)/10)), 0.01);
        set(hs.lb.Model, 'StepSize', a);
        a = str2double(sprintf('%.1g', p(i).ub/10));
        if a<0.01, a = str2double(sprintf('%.1g', abs(p(i).lb/10))); end
        set(hs.ub.Model, 'StepSize', a);
        for j = 1:4 % restore spinner callback
            set(hs.(nam{j}), 'StateChangedCallback', cb{j});
        end
        if strcmpi(get(hs.colorbar, 'Visible'), 'on')
            nii_viewer_cb('colorbar', 1);
        end
    case 'mousedown'
        set(hs.xyz, 'UserData', true);
        mouseClick(hs);
    case 'mouseup'
        set(hs.xyz, 'UserData', false);
    case 'mousemove'
        if get(hs.xyz, 'UserData'), mouseClick(hs); end
        c = get(fh, 'CurrentPoint'); 
        pos = get(hs.panel, 'Position'); height = pos(4);
        pos = get(fh, 'Position');
        sag = get(hs.ax(1), 'Position');
        cor = get(hs.ax(2), 'Position');
        tra = get(hs.ax(3), 'Position');

        c = c(1,1:2) ./ (pos(3:4) - [0 height]);
        if any(c<=0) || any(c>=1) % close to borders
            I = 0;
        elseif c(1)>sag(1) && c(2)>sag(2) % sag
            cp = get(hs.ax(1), 'CurrentPoint');
            I = [hs.spinner(1).getValue cp(1,1:2)];
        elseif c(1)<cor(3) && c(2)>cor(2) % cor
            cp = get(hs.ax(2), 'CurrentPoint');
            I = [cp(1,1) hs.spinner(2).getValue cp(1,2)];
        elseif c(1)<tra(3) && c(2)<tra(4) % tra
            cp = get(hs.ax(3), 'CurrentPoint');
            I = [cp(1,1:2) hs.spinner(3).getValue];
        else
            I = 0;
        end
        if any(I<1) || any(I>hs.dim(1:3)) % out of axis range
            I = cell2mat(get(hs.spinner, 'Value'))';
        end
        set_xyz(hs, I);
    case 'open'
        pf = get(hs.pref, 'UserData');
        [fname, pName] = uigetfile([pf.openPath '/*.nii; *.hdr;' ...
            '*.nii.gz; *.hdr.gz'], 'Select a NIfTI to view');
        if ~ischar(fname), return; end
        fname = fullfile(pName, fname);
        nam = get(fh, 'UserData');
        if strcmp(fname, nam{1}), return; end
        delete(fh);
        nii_viewer(fname);
    case 'add'
        R = [];
        if nargin<2
            pName = get(hs.add, 'UserData');
            [fname, pName] = uigetfile([pName '/*.nii; *.hdr;' ...
                '*.nii.gz; *.hdr.gz'], 'Select overlay NIfTI');
            if ~ischar(fname), return; end
            fname = fullfile(pName, fname);
            if ~isempty(strfind(get(h,'Label'), 'aligned'))
                [mtx, pName] = uigetfile([pName '/*.mat'], ['Select the ' ...
                    'text matrix file which aligns the nii to background']);
                if ~ischar(mtx), return; end
                fid = fopen(fullfile(pName, mtx));
                if fid<0, error('Transformation file not found.'); end
                R = str2num(fread(fid, '*char')');
                fclose(fid);
                if ~isequal(size(R), [4 4])
                    error('Invalid transformation matrix file.');
                end
            end
        else
            fname = varargin{1};
        end
        
        addOverlay(fname, hs, R);
        set(hs.overlay, 'Enable', 'on'); % enable Move/Close overlay
    case 'closeAll'
        for j = 2:numel(hs.q), delete(hs.q(j).hsI); end
        hs.q(2:end) = [];
        p = get(hs.files, 'UserData');
        p(2:end) = [];
        set(hs.files, 'UserData', p);
        
        nam = get(fh, 'UserData');
        nam(2:end) = [];
        set(fh, 'UserData', nam);

        str = get(hs.files, 'String');
        set(hs.files, 'Value', 1, 'String', str(1));
        guidata(fh, hs);
        nii_viewer_cb('files');
        set(hs.overlay, 'Enable', 'off'); % some menu items off
    case {'hdr' 'ext'}
        j = get(hs.files, 'Value');
        if strcmp(cmd, 'hdr')
            hdr = hs.q(j).nii.hdr;
        else
            try 
                hdr = hs.q(j).nii.ext.edata_decoded;
            catch
                errordlg('No known extension for the selected NIfTI'); 
                return; 
            end
        end
        nam = get(hs.files, 'String');
        nam = [nam{j} '_' cmd];
        assignin('base', nam, hdr);
        evalin('base', ['open ' nam]);
    case 'close'
        j = get(hs.files, 'Value');
        if j==1, return; end % no touch to background
        delete(hs.q(j).hsI); % 3 images
        hs.q(j) = [];
        p = get(hs.files, 'UserData');
        p(j) = [];
        set(hs.files, 'UserData', p);
        
        nam = get(fh, 'UserData');
        nam(j) = [];
        set(fh, 'UserData', nam);
        
        str = get(hs.files, 'String');
        str(j) = [];
        n = size(str ,1);
        j = min(n, j);
        set(hs.files, 'Value', j, 'String', str);
        guidata(fh, hs);
        nii_viewer_cb('files');
        on_off = 'on'; if n==1, on_off = 'off'; end
        set(hs.overlay, 'Enable', on_off);
    case 'cross' % show/hide crosshairs and RAS labels
        if strcmp(get(h, 'Checked'), 'on')
            set(h, 'Checked', 'off');
            set([hs.cross hs.ras], 'Visible', 'off');
        else
            set(h, 'Checked', 'on');
            set([hs.cross hs.ras], 'Visible', 'on');
        end
    case 'color'
        c = uisetcolor([], 'Pick crosshair color');
        if numel(c) ~= 3, return; end
        set([hs.cross hs.ras], 'Color', c);
    case 'thickness'
        c = strtok(get(h, 'Label'));
        set(hs.cross, 'LineWidth', str2double(c));
    case 'gap'
        c = str2double(strtok(get(h, 'Label')));
        hs.gap = min(hs.pixdim) ./ hs.pixdim * c / 2;
        guidata(fh, hs);
        nii_viewer_cb('spin_x');
        nii_viewer_cb('spin_y');
        nii_viewer_cb('spin_z');
    case 'copy'
        set(hs.panel, 'Visible', 'off');
        clnObj = onCleanup(@() set(hs.panel, 'Visible', 'on'));
        pf = get(hs.pref, 'UserData');
        print(fh, '-dbitmap', '-noui', ['-r' pf.dpi]);
        % print(fh, '-dmeta', '-painters');
    case 'save'
        ext = get(h, 'Label');
        fmt = ext;
        if strcmp(ext, 'jpg'), fmt = 'jpeg';
        elseif strcmp(ext, 'tif'), fmt = 'tiff';
        elseif strcmp(ext, 'eps'), fmt = 'epsc';
        elseif strcmp(ext, 'emf'), fmt = 'meta';
        end
        [fname, pName] = uiputfile(['*.' ext], 'Input file name to save figure');
        if ~ischar(fname), return; end
        fname = fullfile(pName, fname);
        if any(strcmp(ext, {'eps' 'pdf' 'emf'})), render = '-painters';
        else render = '-opengl';
        end
        pf = get(hs.pref, 'UserData');
        set(hs.panel, 'Visible', 'off');
        clnObj = onCleanup(@() set(hs.panel, 'Visible', 'on'));
        % print(fh, fname, render, ['-d' fmt]);
        print(fh, fname, render, '-noui', ['-d' fmt], ['-r' pf.dpi]);
    case 'colorbar'
        if nargin<2 && strcmpi(get(hs.colorbar, 'Visible'), 'on')
            set(hs.colorbar, 'Visible', 'off'); 
            if nargin<2, set(h, 'Checked', 'off'); end
            return; 
        end
        p = get(hs.files, 'UserData');
        i = get(hs.files, 'Value');
        if p(i).show==0, return; end
        lut = p(i).lut;
        if lut == 11, map = lut2map(get(hs.lut, 'UserData'));
        else map = lut2map(lut);
        end
        rg = sort([p(i).lb p(i).ub]);
        if lut~=10
            if rg(2)<0, rg = rg([2 1]); end
            mn = str2double(num2str(mean(rg), '%.4g'));
            labls = [rg(1) mn rg(2)];
        else
            rg = sort(abs(rg));
            labls = {num2str(-rg(2)) num2str(rg(1),'+/-%g') num2str(rg(2))};
        end
        % colormap in earlier matlab version changes values in colorbar.
        % So we have to turn on it first, and set those values each time.
        % set(hs.colorbar, 'Visible', 'on', 'YTickLabel', labls); % new matlab
        set(hs.colorbar, 'Visible', 'on');
        colormap(hs.ax(4), map);
        a = get(hs.colorbar, 'Children'); set(a, 'YData', [0 1]); % Tricky!
        set(hs.colorbar, 'YTickLabel', labls, 'YTick', [0 0.5 1], 'Ylim', [0 1]);
        if nargin<2, set(h, 'Checked', 'on'); end
    case 'about'
        str = sprintf(['nii_viewer.m by Xiangrui Li\n\n' ...
            'Feedback to: xiangrui.li@gmail.com\n\n' ...
            'Last updated on 20%s\n'], reviseDate);
        helpdlg(str, 'About nii_viewer')
    case 'stack'
        i = get(hs.files, 'Value');
        if i==1, return; end % no-op for background
        c = get(h, 'Label');
        n = numel(hs.q);
        ind = 1:n;
        if ~isempty(strfind(c, 'up')) % one level up
            if i==n, return; end
            for j = 1:3, uistack(hs.q(i).hsI(j)); end
            ind = ind([1:i-1 i+1 i i+2:n]);
        elseif ~isempty(strfind(c, 'down')) % one level down
            if i==2, return; end
            for j = 1:3, uistack(hs.q(i).hsI(j), 'down'); end
            ind = ind([1:i-2 i i-1 i+1:n]);
        elseif ~isempty(strfind(c, 'top')) % top within background + overlays
            step = n-i;
            if step==0, return; end
            for j = 1:3, uistack(hs.q(i).hsI(j), 'up', step); end
            ind = ind([1:i-1 i+1:n i]);
        elseif ~isempty(strfind(c, 'bottom')) % bottom for overlays
            step = i-2;
            if step==0, return; end
            for j = 1:3, uistack(hs.q(i).hsI(j), 'down', step); end
            ind = ind([1 i 2:i-1 i+1:n]);
        else error('Unknown stack level: %s', c);
        end
        
        str = get(hs.files, 'String');
        set(hs.files, 'String', str(ind), 'Value', find(ind==i));
        p = get(hs.files, 'UserData');
        p = p(ind);
        set(hs.files, 'UserData', p);
        
        nam = get(fh, 'UserData');
        nam = nam(ind);
        set(fh, 'UserData', nam);

        hs.q = hs.q(ind);
        guidata(fh, hs);
    case 'saveas' % this uses nii_tool, and not related to viewer
        i = get(hs.files, 'Value');
        c = get(h, 'Label');
        nam = get(fh, 'UserData');
        pName = fileparts(nam{i});
        
        if ~isempty(strfind(c, 'dim 4')) % fsl RGB
            nii = nii_tool('load', nam{i}); % re-load to be safe
            if any(size(nii.img,8) == 3:4)
                nii.img = permute(nii.img, [1:3 8 4:7]);
            elseif ~any(nii.hdr.dim(5) == 3:4)
                errordlg('Selected image is not RGB data.'); return;
            end
            [fname, pName] = uiputfile([pName '/*.nii'], ...
                'Input name for FSL RGB file');
            if ~ischar(fname), return; end
            fname = fullfile(pName, fname);
            nii_tool('save', nii, fname);
        elseif ~isempty(strfind(c, 'dim 3')) % old mricron RGB
            nii = nii_tool('load', nam{i});
            if any(nii.hdr.dim(5) == 3:4)
                nii.img = permute(nii.img, [1:3 5:7 4]);
            elseif ~any(size(nii.img,8) == 3:4)
                errordlg('Selected image is not RGB data'); return;
            end
            [fname, pName] = uiputfile([pName '/*.nii'], ...
                'Input name for old mricrom styte file');
            if ~ischar(fname), return; end
            fname = fullfile(pName, fname);
            old = nii_tool('RGBStyle', 'mricron');
            nii_tool('save', nii, fname);
            nii_tool('RGBStyle', old);
        elseif ~isempty(strfind(c, 'AFNI')) % NIfTI RGB
            nii = nii_tool('load', nam{i});
            if any(nii.hdr.dim(5) == 3:4)
                nii.img = permute(nii.img, [1:3 5:7 4]);
            elseif ~any(size(nii.img,8) == 3:4)
                errordlg('Selected image is not RGB data'); return;
            end
            [fname, pName] = uiputfile([pName '/*.nii'], ...
                'Input name for NIfTI standard RGB file');
            if ~ischar(fname), return; end
            fname = fullfile(pName, fname);
            old = nii_tool('RGBStyle', 'afni');
            nii_tool('save', nii, fname);
            nii_tool('RGBStyle', old);
        elseif ~isempty(strfind(c, '3D')) % SPM 3D
            nii = nii_tool('load', nam{i});
            if nii.hdr.dim(5)<2
                errordlg('Selected image is not multi-volume data'); return;
            end
            [fname, pName] = uiputfile([pName '/*.nii'], ...
                'Input base name for SPM 3D file');
            if ~ischar(fname), return; end
            fname = fullfile(pName, fname);
            nii_tool('save', nii, fname, 1); % force 3D
        elseif ~isempty(strfind(c, 'new resolution'))
            str = 'Resolution for three dimension in mm:'; 
            a = inputdlg(str, 'Input space resolution', 1, {'3 3 3'});
            if isempty(a), return; end
            res = sscanf(a{1}, '%g %g %g');
            if numel(res) ~= 3
                errordlg('Invalid space resolution');
                return;
            end
            if isequal(res, hs.q(i).nii.hdr.pixdim(2:4))
                warndlg('The input resolution is the same as current one');
                return;
            end
            [fname, pName] = uiputfile([pName '/*.nii;nii.gz'], ...
                'Input result name for the new resolution file');
            if ~ischar(fname), return; end
            fname = fullfile(pName, fname);
            pf = get(hs.pref, 'UserData');
            nii_xform(nam{i}, res, fname, pf.interp, pf.extraV)
        elseif ~isempty(strfind(c, 'matching background'))
            if i == 1
                errordlg('You seleted background image');
                return; 
            end
            [fname, pName] = uiputfile([pName '/*.nii;*.nii.gz'], ...
                'Input result file name');
            if ~ischar(fname), return; end
            fname = fullfile(pName, fname);
            pf = get(hs.pref, 'UserData');
            nii_xform(nam{i}, nam{1}, fname, pf.interp, pf.extraV)
        elseif ~isempty(strfind(c, 'aligned template'))
            [temp, pName] = uigetfile([pName '/*.nii;*.nii.gz'], ...
                'Select the aligned template file');
            if ~ischar(temp), return; end
            temp = fullfile(pName, temp);
            [mtx, pName] = uigetfile([pName '/*.mat'], ['Select the text ' ...
                'matrix file which aligns the nii to the template']);
            if ~ischar(mtx), return; end
            mtx = fullfile(pName, mtx);
            [fname, pName] = uiputfile([pName '/*.nii;*.nii.gz'], ...
                'Input result file name');
            if ~ischar(fname), return; end
            fname = fullfile(pName, fname);
            pf = get(hs.pref, 'UserData');
            nii_xform(nam{i}, {temp mtx}, fname, pf.interp, pf.extraV)
        else
            errordlg(sprintf('%s not implemented yet.', c));
        end
    case 'zoom'
        m = str2double(get(h, 'Label'));
        set_zoom(m, hs);
    case 'mask'
        i = get(hs.files, 'Value');
        nam = get(fh, 'UserData');
        pName = fileparts(nam{i});
        [fname, pName] = uigetfile([pName '/*.nii;*.hdr;*.nii.gz;*.hdr.gz'], ...
            'Select mask NIfTI');
        if ~ischar(fname), return; end
        nii = nii_tool('load', fullfile(pName, fname));
        if ~any([nii.hdr.sform_code nii.hdr.qform_code] == hs.form_code)
            str = ['The mask coordinate systems are inconsistent with the ' ...
                'selected image. Do you want to apply the mask anyway?'];
            btn = questdlg(str, 'Apply mask?', 'Cancel', 'Apply', 'Cancel');
            if isempty(btn) || strcmp(btn, 'Cancel'), return; end
        end
        R = nii_xform_mat(nii.hdr, hs.form_code);
        dim = single(size(hs.q(i).nii.img)); % may worth to use single
        dim(numel(dim)+1:3) = 1; dim = dim(1:3);
        [Y, X, Z] = meshgrid(1:dim(2), 1:dim(1), 1:dim(3));
        I = [X(:) Y(:) Z(:)]'-1; I(4,:) = 1;
        I = R \ (hs.q(i).R * I) + 1; % ijk+1 for mask
        I = round(I * 100) / 100;

        im = single(nii.img(:,:,:,1,1,1,1,1)); % first volume
        slope = nii.hdr.scl_slope;
        if slope==0, slope = 1; end
        im = im * slope + nii.hdr.scl_inter;
        im = interp3a(im, I, 'nearest');
        im1 = im(~isnan(im));
        im = reshape(im, dim);
        if numel(unique(im1(:)))<3
            thre = min(im(:));
        else
            rg = get_range(im);
            str = sprintf('Threshold for non-binary mask (%.3g to %.4g)', ...
                min(im(:)), max(im(:)));
            a = inputdlg(str,'Input mask threshold', 1, {num2str(rg(1),'%.3g')});
            if isempty(a), return; end
            thre = str2double(a{1});
        end
        hs.q(i).nii.mask = abs(im)>thre;
        guidata(fh, hs);
        nii_viewer_cb('update');
    case 'pref'
        pref_dialog(hs.pref);
    case 'background'
        if strcmp(get(h, 'Checked'), 'on');
            set(h, 'Checked', 'off');
            set(hs.im_panel, 'BackgroundColor', 'k');
            set(hs.colorbar, 'EdgeColor', 'w');
        else
            set(h, 'Checked', 'on');
            set(hs.im_panel, 'BackgroundColor', 'w');
            set(hs.colorbar, 'EdgeColor', 'k');
        end
        nii_viewer_cb('update');
    case 'flipLR'
        if strcmp(get(h, 'Checked'), 'on');
            set(h, 'Checked', 'off');
            set(hs.ax([2 3]), 'XDir', 'normal');
            set(hs.ras([3 5]), 'String', 'L');
        else
            set(h, 'Checked', 'on');
            set(hs.ax([2 3]), 'XDir', 'reverse');
            set(hs.ras([3 5]), 'String', 'R');
        end
    case 'keyHelp'
        str = sprintf(['Left or Right: Move crosshair left or right\n\n' ...
                       'Up or Down: Move crosshair superior or inferior\n\n' ...
                       '[ or ]: Move crosshair posterior or anterior\n\n' ...
                       '< or >: Decrease or increase volume number\n\n' ...
                       'Ctrl + or -: Zoom in or out by 10%%\n']);
        helpdlg(str, 'Key Shortcut');
    otherwise
        error('Unknown Callback: %s', cmd);
end

%% zoom in/out with a factor
function set_zoom(m, hs)
if m == 1 % restore, regardless of crosshair loc
    axis(hs.ax(1), [0 hs.dim(2) 0 hs.dim(3)]+0.5);
    axis(hs.ax(2), [0 hs.dim(1) 0 hs.dim(3)]+0.5);
    axis(hs.ax(3), [0 hs.dim(1) 0 hs.dim(2)]+0.5);
    return;
end

lim = zeros(3,2);
for ix = 1:3
    c = hs.spinner(ix).getValue;
    lim(ix,:) = c + 0.5 + [-1 1]/2*hs.dim(ix)/m;
end
axis(hs.ax(1), [lim(2,:) lim(3,:)]);
axis(hs.ax(2), [lim(1,:) lim(3,:)]);
axis(hs.ax(3), [lim(1,:) lim(2,:)]);

%% WindowKeyPressFcn
function KeyPressFcn(fh, evt)
key = evt.Character;
if isempty(key), return; end
hs = guidata(fh);
ctl = evt.Modifier;
if isempty(ctl)
    if key == 28 % left
        val = max(hs.spinner(1).getValue-1, 1);
        hs.spinner(1).setValue(val);
    elseif key == 29 % right
        val = min(hs.spinner(1).getValue+1, hs.dim(1));
        hs.spinner(1).setValue(val);
    elseif key == 30 % up
        val = min(hs.spinner(3).getValue+1, hs.dim(3));
        hs.spinner(3).setValue(val);
    elseif key == 31 % down
        val = max(hs.spinner(3).getValue-1, 1);
        hs.spinner(3).setValue(val);
    elseif key == ']' 
        val = min(hs.spinner(2).getValue+1, hs.dim(2));
        hs.spinner(2).setValue(val);
    elseif key == '['
        val = max(hs.spinner(2).getValue-1, 1);
        hs.spinner(2).setValue(val);
    elseif key == '.' % >
        val = min(hs.volume.getValue+1, get(hs.volume.Model,'Maximum'));
        hs.volume.setValue(val);
    elseif key == ',' % <
        val = max(hs.volume.getValue-1, 1);
        hs.volume.setValue(val);
    end
elseif strcmpi(ctl, 'control') || strcmpi(ctl, 'command')
    lim = get(hs.ax(1), 'XLim');
    m = hs.dim(2) / abs(lim(2)-lim(1));
    if key=='+' || key=='='
        if m>40, return; end % ignore
        m = m * 1.1;
    elseif key == '-'
        if m<=1, return; end
        m = m / 1.1;
        if m<1.01, m = 1; end
    end
    set_zoom(m, hs);
end

%% Drop callback: to open as background
function javaDropFcn(~, evt)
try
    fname = evt.Data{1};
    nii_tool('hdr', fname); % make sure valid nii
    set(0, 'ShowHiddenHandles', 'on');
    delete(gcf);
    nii_viewer(fname);
catch me
    errordlg(me.message);
end

%% update CData/AlphaData for one of the sag/cor/tra view
function set_cdata(ix, fh, ind)
hs = guidata(fh); % read it here to catch change within previous set_cdata call
if nargin<3, ind = hs.spinner(ix).getValue; end
ind = round(ind);
if ind<1 || ind>hs.dim(ix), return; end
interStr = get(hs.interp, 'String');
p = get(hs.files, 'UserData');
n = numel(p);
for i = 1:n
    if ~p(i).show, continue; end % save time, but need to update when enabled
    lut = p(i).lut;
    if lut == 11
        vector_lines(hs, i, ix); continue; 
    elseif ~strcmpi(get(hs.q(i).hsI(1), 'Type'), 'image')
        delete(hs.q(i).hsI); % delete quiver
        for j = 1:3, hs.q(i).hsI(j) = copyobj(hs.q(1).hsI(j), hs.ax(j)); end
        guidata(fh, hs);
        crossFront(hs);
        if i<n, for j=1:3; uistack(hs.q(i).hsI(j), 'down', n-i); end; end
    end
    t = round(p(i).volume);
    img = permute(hs.q(i).nii.img(:,:,:,t,:,:,:,:), [1:3 8 4:7]);
    img = single(img);
    dim4 = size(img, 4); % in case of RGB
    if hs.q(i).scl_slope~=1 || hs.q(i).scl_inter~=0
        img = img * hs.q(i).scl_slope + hs.q(i).scl_inter;
    end
    
    if isfield(hs.q(i).nii, 'mask')
        img = bsxfun(@times, img, hs.q(i).nii.mask);
    end
    
    i0 = ind; % background or the same R as background
    if hs.q(i).interp % skip if false: save time and avoid numeric inaccuracy
        dim = hs.dim;
        dim(ix) = 1; % 1 slice at dim ix
        [Y, X, Z] = meshgrid(1:dim(2), 1:dim(1), 1:dim(3));
        I = [X(:) Y(:) Z(:)]' - 1; % ijk grids in figure
        I(ix,:) = ind-1; I(4,:) = 1;
        I = hs.q(i).R \ (hs.q(1).R * I) + 1; % ijk+1 for overlay with fraction
 
        im = zeros([dim dim4], 'single');
        for j = 1:dim4
            if p(i).smooth
                nv = 3; % number of voxels (odd) used for smooth
                d3 = dim; d3(ix) = nv; 
                b = zeros(d3, 'single');
                vec = {':' ':' ':'};
                I0 = I(1:3,:);
                for k = 1:nv
                    I0(ix,:) = I(ix,:) - (nv+1)/2 + k;
                    a = interp3a(img(:,:,:,j), I0, interStr{p(i).interp});
                    vec{ix} = k; b(vec{:}) = reshape(a, dim);
                end
                b = smooth23(b, 'gaussian', nv);
                vec{ix} = (nv+1)/2;
                im(:,:,:,j) = b(vec{:}); % middle one
            else
                a = interp3a(img(:,:,:,j), I, interStr{p(i).interp});
                im(:,:,:,j) = reshape(a, dim);
            end
        end
        img = im; i0 = 1; % overwrite img with interpolated single slice
    elseif p(i).smooth
        nv = 3; % odd number
        vec = {':' ':' ':'};
        a = ind - (nv+1)/2 + (1:nv);
        if any(a<1 | a>hs.dim(ix)), a = ind; end
        vec{ix} = a;
        a = smooth23(img(vec{:}), 'gaussian', nv);
        vec{ix} = mean(1:size(a,ix)); % middle slice
        img = a(vec{:}); i0 = 1;
    end

    if     ix == 1, im = permute(img(i0,:,:,:), [3 2 4 1]);
    elseif ix == 2, im = permute(img(:,i0,:,:), [3 1 4 2]);
    elseif ix == 3, im = permute(img(:,:,i0,:), [2 1 4 3]);
    end
    
    if dim4 == 1 % not RGB
        rg = sort([p(i).lb p(i).ub]);
        if rg(2)<0 % asking for negative data
            rg = -rg([2 1]);
            if lut~=10, im = -im; end
        end
        if lut == 10 % two-sided, store negative value
            rg = sort(abs(rg));
            im_neg = -single(im) .* (im<0);
            im_neg = (im_neg-rg(1)) / (rg(2)-rg(1));
            im_neg(im_neg>1) = 1; im_neg(im_neg<0) = 0; 
            im_neg = repmat(im_neg, [1 1 3]); % gray now
        end
        im = (single(im)-rg(1)) / (rg(2)-rg(1));
        im(im>1) = 1; im(im<0) = 0;
        im = repmat(im, [1 1 3]); % gray now

        if     lut == 1, % gray do nothing
        elseif lut == 2, im(:,:,2:3) = 0; % red
        elseif lut == 3, im(:,:,[1 3]) = 0; % green
        elseif lut == 4, im(:,:,1:2) = 0; % blue
        elseif lut == 5, im(:,:,2) = 0; % violet
        elseif lut == 6, im(:,:,3) = 0; % yellow
        elseif lut == 7, im(:,:,1) = 0; % cyan
        elseif lut == 8 % red_yellow
            im(:,:,3) = 0;
            a = im(:,:,1); a(a>0) = 1; im(:,:,1) = a;
        elseif lut == 9 % blue_green
            im(:,:,1) = 0;
            a = im(:,:,3); a(a==0) = 1; a = 1 - a; im(:,:,3) = a;
        elseif lut == 10 % two-sided: combine red_yellow & blue_green
            im(:,:,3) = 0;
            a = im(:,:,1); a(a>0) = 1; im(:,:,1) = a;
            im_neg(:,:,1) = 0;
            a = im_neg(:,:,3); a(a==0) = 1; a = 1 - a; im_neg(:,:,3) = a;
            im = im + im_neg;
        end
    else % RGB
        if max(im(:))>2, im = single(im) / 255; end % guess uint8
        im(im>1) = 1; im(im<0) = 0;
    end
    set(hs.q(i).hsI(ix), 'CData', im);
    
    alfa = mean(im,3);
    % This simple 2D extraction avoids seeing white background
    if i==1 && isequal(get(hs.im_panel,'BackgroundColor'), [1 1 1])
        m = mean(alfa(alfa(:)>0));
        alfa = smooth23(alfa/m, 'box', 5)> m/2;
        a = alfa;
        for j = 1:size(alfa,1) % left-right boundary
            i1 = find(alfa(j,:), 1);
            i2 = find(alfa(j,:), 1, 'last'); 
            alfa(j, i1:i2) = 1;
        end
        for j = 1:size(a,2) % top-bottom boundary
            i1 = find(a(:,j), 1);
            i2 = find(a(:,j), 1, 'last'); 
            a(i1:i2, j) = 1;
        end
        alfa = alfa .* a;
    end
    alfa = p(i).alpha * alfa>0;
    set(hs.q(i).hsI(ix), 'AlphaData', alfa);
end
 
%% Get mouse location, and set it to spinner
function mouseClick(hs)
ax = gca;
if isequal(ax, hs.ax(4)), return; end
c = get(ax, 'CurrentPoint');
c = round(c(1, 1:2));
x = get(ax, 'XLim'); y = get(ax, 'YLim'); 
if c(1)<x(1) || c(1)>x(2) || c(2)<y(1) || c(2)>y(2), return; end

i = 1:3;
i(ax==hs.ax(1:3)) = [];
hs.spinner(i(1)).setValue(c(1)); % evoke update by spinner
hs.spinner(i(2)).setValue(c(2));

%% Add an overlay
function addOverlay(fname, hs, mtx)
hdr = nii_tool('hdr', fname);
codes = [hs.q(1).nii.hdr.sform_code hs.q(1).nii.hdr.qform_code];
frm = 0;
addAligned = nargin>2 && ~isempty(mtx);
if hdr.sform_code>0 && any(hdr.sform_code == codes)
    frm = hdr.sform_code;
elseif hdr.qform_code>0 && any(hdr.qform_code == codes)
    frm = hdr.qform_code;
elseif ~addAligned
    warndlg(['The coordinate systems are inconsistent for background image' ...
        ' and the overlay image. The overlay is likely meaningless.'], ...
        'Transform Inconsistent');
end
i = numel(hs.q) + 1;

if addAligned % aligned mtx: do it in special way
    [hs.q(i), frm, rg, dim, pixdim] = read_nii(fname, frm, 0); % no re-orient
    R0 = nii_xform_mat(hs.q(1).nii.hdr, hs.form_code); % original background R
    
    % see nii_xform for more comment on following method
    R = R0 / diag([hs.q(1).nii.hdr.pixdim(2:4) 1]) * mtx * diag([pixdim 1]);
    [~, i1] = max(abs(hs.q(i).R(1:3,1:3)));
    [~, i0] = max(abs(R(1:3,1:3)));
    flp = sign(R(i0+[0 4 8])) ~= sign(hs.q(i).R(i1+[0 4 8]));
    if any(flp)
        rotM = diag([1-flp*2 1]);
        rotM(1:3,4) = (dim-1).* flp;
        R = R / rotM;
    end
    hs.q(i).R = R; % no re-orient or flip: no problem since not for DTI lines
else
    [hs.q(i), frm, rg] = read_nii(fname, frm);
end
if frm>0 && frm ~= hs.form_code % update background form_code
    hs.form_code = frm;
    hs.q(1).R = nii_xform_mat(hs.q(1).nii.hdr, frm); % img re-orient done
end
hs.q(i).interp = any(abs(hs.q(1).R(:) - hs.q(i).R(:))>1e-3);

% duplicate image obj for overlay
for j = 1:3, hs.q(i).hsI(j) = copyobj(hs.q(1).hsI(j), hs.ax(j)); end
crossFront(hs);

% set default for overlay
p = get(hs.files, 'UserData');
p(i).show = 1;
p(i).lb = rg(1);
p(i).ub = rg(2);
p(i).lut = mod(p(i-1).lut, 7) + 1; % 1:7, use next lut of previous
p(i).alpha = 1;
p(i).smooth = 0;
p(i).interp = 1;
p(i).volume = 1;
set(hs.files, 'UserData', p);

[pName, niiName, ext] = fileparts(fname);
if strcmpi(ext, '.gz'), [~, niiName] = fileparts(niiName); end
str = get(hs.files, 'String');
str{i} = niiName; % add the bottom of the list
set(hs.files, 'String', str, 'Value', i);
set(hs.add, 'UserData', pName);

nam = get(hs.fig, 'UserData');
nam{i} = fname;
set(hs.fig, 'UserData', nam);

guidata(hs.fig, hs); % update hs.q for nii
nii_viewer_cb('files'); % show para for new overlay
nii_viewer_cb('update');

% save addPath to file
pf = get(hs.pref, 'UserData');
pf.addPath = pName;
para = load(pf.pref_file); para = para.para;
para.nii_viewer = pf;
try save(pf.pref_file, 'para'); catch, end

%% Load, re-orient nii, return essential nii stuff
% nii.img may be re-oriented, but nii.hdr is not touched
function [q, frm, rg, dim, pixdim] = read_nii(fname, ask_code, reOri)
q.hsI = []; % assign the field for the struct
q.interp = false; % default of background
q.nii = nii_tool('load', fname);
q.scl_slope = q.nii.hdr.scl_slope;
if q.scl_slope==0, q.scl_slope = 1; end
q.scl_inter = q.nii.hdr.scl_inter;
ndim = q.nii.hdr.dim(1);
dim = q.nii.hdr.dim(2:8);
dim(dim<1 | dim>32767 | mod(dim,1)>0) = 1;
if ndim>4 % 4+ dim, put all into dim4
    if sum(dim(4:7)>1)>1
        warndlg([fname ' has 5 or more dimension. Dimension above 4 are ' ...
            'all treated as volumes for visualization']);        
    end
    dim(4) = prod(dim(4:7)); dim(5:7) = 1;
    q.nii.img = reshape(q.nii.img, [dim size(q.nii.img, 8)]);
end

if nargin<2, ask_code = []; end
[q.R, frm] = nii_xform_mat(q.nii.hdr, ask_code);
dim = dim(1:3);
pixdim = q.nii.hdr.pixdim(2:4);
q.flip = false(1,3);
if nargin<3 || reOri
    [~, ixyz] = max(abs(q.R(1:3,1:3)));
    [~, perm] = sort(ixyz);
    if ~isequal(perm, 1:3) && isequal(sort(perm), 1:3)
        dim = dim(perm);
        pixdim = pixdim(perm);
        q.R(:,1:3) = q.R(:,perm); 
        q.nii.img = permute(q.nii.img, [perm 4:8]);
    end

    q.flip = q.R([1 6 11]) < 0;
    rotM = diag([1-q.flip*2 1]); % 1 or -1 on diagnal
    rotM(1:3, 4) = (dim-1) .* q.flip; % 0 or dim-1
    q.R = q.R / rotM; % xform matrix after flip
    for i = 1:3, if q.flip(i), q.nii.img = flipdim(q.nii.img, i); end; end %#ok
end
if size(q.nii.img,4)<4 && ~isfloat(q.nii.img)
    q.nii.img = single(q.nii.img); 
end
if q.scl_slope~=1 || q.scl_inter~=0 && isfloat(q.nii.img)
    q.nii.img = q.nii.img  * q.scl_slope + q.scl_inter;
    q.scl_slope = 1; q.scl_inter = 0;
end

img = single(q.nii.img(:,:,:,1,1,1,1,1)) * q.scl_slope + q.scl_inter;
rg = get_range(img);

%% Return xform mat
function [R, frm] = nii_xform_mat(hdr, ask_code)
if nargin<2 || isempty(ask_code), ask_code = hdr.sform_code; end
if hdr.sform_code == ask_code
    frm = hdr.sform_code;
    R = [hdr.srow_x; hdr.srow_y; hdr.srow_z; 0 0 0 1];
elseif hdr.qform_code == ask_code
    frm = hdr.qform_code;
    R = quat2R(hdr);
elseif hdr.sform_code>0
    frm = hdr.sform_code;
    R = [hdr.srow_x; hdr.srow_y; hdr.srow_z; 0 0 0 1];
elseif hdr.qform_code>0
    frm = hdr.qform_code;
    R = quat2R(hdr);
else
    frm = 0;
    R = diag(hdr.pixdim(2:4));
    R(:,4) = -hdr.pixdim(2:4) .* hdr.dim(2:4) / 2;
    R = [R; 0 0 0 1];
end

%% quatenion to xform_mat
function R = quat2R(hdr)
b = hdr.quatern_b;
c = hdr.quatern_c;
d = hdr.quatern_d;
a = sqrt(1-b^2-c^2-d^2);
R = [1-2*(c^2+d^2)  2*(b*c-d*a)     2*(b*d+c*a);
     2*(b*c+d*a)    1-2*(b^2+d^2)   2*(c*d-b*a);
     2*(b*d-c*a )   2*(c*d+b*a)     1-2*(b^2+c^2)];
R = R * diag(hdr.pixdim(2:4));
if hdr.pixdim(1)<0, R(:,3)= -R(:,3); end
R = [R [hdr.qoffset_x hdr.qoffset_y hdr.qoffset_z]'; 0 0 0 1];

%% Create java SpinnerNumber
function h = java_spinner(pos, val, parent, callback, fmt, helpTxt)
mdl = javax.swing.SpinnerNumberModel(val(1), val(2), val(3), val(4));
jSpinner = com.mathworks.mwswing.MJSpinner(mdl);
h = javacomponent(jSpinner, pos, parent);
set(h, 'StateChangedCallback', callback, 'ToolTipText', helpTxt);
jEditor = javaObject('javax.swing.JSpinner$NumberEditor', h, fmt);
h.setEditor(jEditor);

%% estimate intensity for lower and upper bound of display
function rg = get_range(img)
if size(img,8)>2
    if max(img(:))>2, rg = [0 255];
    else rg = [0 1];
    end
    return;
end
ind = abs(img(:))>50;
if sum(ind)<numel(img)/10, ind = abs(img(:))>std(img(:))/2; end
im = img(ind);
mn = mean(im);
sd = std(im);
rg = mn + [-2 2]*sd;
mi = min(img(:)); ma = max(img(:));
if rg(1)<=0 && mn-sd>0, rg(1) = sd/5; end
if rg(1)<mi || isnan(rg(1)), rg(1) = mi; end
if rg(2)>ma || isnan(rg(2)), rg(2) = ma; end
if rg(1)==rg(2), rg(1) = mi; end
rg = str2num(sprintf('%.2g ', rg)); %#ok<*ST2NM>

%% Get the last date string in history
function dStr = reviseDate
dStr = '151107?';
fid = fopen(which(mfilename));
if fid<1, return; end
str = fread(fid, '*char')';
fclose(fid);
ind = strfind(str, '% End of history. Don''t edit this line!');
if isempty(ind), return; end
ind = ind(1);
ret = str(ind-1); % new line char: \r or \n
str = str(max(1, ind-500):ind+2); % go back several lines
ind = strfind(str, [ret '% ']); % new line with % and space
for i = 1:numel(ind)-1
    ln = str(ind(i)+3 : ind(i+1)-1);
    if numel(ln)>5 && all(isstrprop(ln(1:6), 'digit'))
        dStr = ln(1:6);
    end 
end

%% Draw vector lines
function vector_lines(hs, i, ix)
if strcmpi(get(hs.q(i).hsI(1), 'Type'), 'image')
    delete(hs.q(i).hsI);
    lut = get(hs.lut, 'UserData'); % last lut
    if isempty(lut), lut = 2; end % default red
    clr = lut2map(lut); clr = clr(end,:);
    for j = 1:3
        hs.q(i).hsI(j) = quiver(hs.ax(j), 1, 1, 0, 0, 'Color', clr, ...
            'ShowArrowHead', 'off', 'AutoScale', 'off');
    end
    crossFront(hs); % to be safe before next
    iDown = numel(hs.q) - i;
    if iDown>0, for j = 1:3, uistack(hs.q(i).hsI(j), 'down', iDown); end; end
    guidata(hs.fig, hs);
end

img = hs.q(i).nii.img;
% This is needed since vec is in image ref, at least for fsl
for j = 1:3, if hs.q(i).flip(j), img(:,:,:,j) = -img(:,:,:,j); end; end
if isfield(hs.q(i).nii, 'mask')
    img = bsxfun(@times, img, hs.q(i).nii.mask);
end

I = hs.spinner(ix).getValue;
if     ix==1, im = img(I,:,:,[2 3]); im = permute(im, [3 2 4 1]);
elseif ix==2, im = img(:,I,:,[1 3]); im = permute(im, [3 1 4 2]);
elseif ix==3, im = img(:,:,I,[1 2]); im = permute(im, [2 1 4 3]);
end
im(im==0) = nan; % avoid dots in emf and eps
dim = single(size(im));
[X, Y] = meshgrid(1:dim(2), 1:dim(1));
X = X - im(:,:,1)/2;
Y = Y - im(:,:,2)/2;
set(hs.q(i).hsI(ix), 'XData', X, 'YData', Y, 'UData', im(:,:,1), 'VData', im(:,:,2));

%% Bring cross and label to front
function crossFront(hs)
for i = 1:3
    txt = allchild(hs.ax(i));
    ind = strcmp(get(txt, 'Type'), 'text');
    txt = txt(ind); % two letters, plus junk text with matlab 2010b
    uistack([txt' hs.cross(i)], 'top');
end

%% Compute color map for LUT
function map = lut2map(lut)
map = linspace(0,1,64)'*[1 1 1]; % gray
if     lut == 1, return; % gray
elseif lut == 2, map(:,2:3) = 0; % red
elseif lut == 3, map(:,[1 3]) = 0; % green
elseif lut == 4, map(:,1:2) = 0; % blue
elseif lut == 5, map(:,2) = 0; % violet
elseif lut == 6, map(:,3) = 0; % yellow
elseif lut == 7, map(:,1) = 0; % cyan
elseif lut == 8, map(:,3) = 0; map(:,1) = 1; % red_yellow
elseif lut == 9, map(:,1) = 0; map(:,3) = map(end:-1:1,3); % blue_green
elseif lut == 10 % two-sided
    map = map(1:2:end,:); % half
    map_neg = map;
    map(:,3) = 0; map(:,1) = 1; % red_yellow
    map_neg(:,1) = 0; map_neg(:,3) = map_neg(end:-1:1,3); % blue_green
    map = [map_neg(end:-1:1,:); map];
elseif lut == 11, map(:,2:3) = 0; % vector lines
else error('Unknown LUT: lut = %g', lut);
end

%% Preference dialog
function pref_dialog(pref)
pf = get(pref, 'UserData');
if ~isfield(pf, 'dpi'), pf.dpi = '0'; end
    
d = dialog('Name', 'Preferences', 'Visible', 'off');
pos = get(d, 'Position');
pos(3:4) = [396 230];
set(d, 'Position', pos, 'Visible', 'on');

uicontrol('Parent', d, 'Position', [300 10 70 24], 'Tag', 'OK', ...
    'String', 'OK', 'Callback', {@pref_dialog_cb, pref});
uicontrol('Parent', d, 'Position',[200 10 70 24], ...
    'String', 'Cancel', 'Callback', 'delete(gcf)');

uicontrol('Parent', d, 'Style', 'text', 'Position', [8 200 300 22], ...
    'String', 'Background (template) image folder:', 'HorizontalAlignment', 'left');
h.openPath = uicontrol(d, 'Style', 'edit', 'String', pf.openPath, ...
    'Position', [8 180 350 22], 'BackGroundColor', 'w', 'HorizontalAlignment', 'left', ...
    'TooltipString', 'nii_viewer will point to this folder when you "Open" image');
uicontrol('Parent', d, 'Position', [358 181 30 22], 'Tag', 'browse', ...
    'String', '...', 'Callback', @pref_dialog_cb);

uipanel(d, 'Units', 'Pixels', 'Position', [4 110 390 56], 'BorderType', 'etchedin', ...
    'BorderWidth', 2, 'Title', 'For "Save NIfTI as" if interpolation is applicable');
str = {'nearest' 'linear' 'cubic' 'spline'};
val = find(strcmp(str, pf.interp));
uicontrol('Parent', d, 'Style', 'text', 'Position', [8 116 140 22], ...
    'String', 'Interpolation method:', 'HorizontalAlignment', 'right');
h.interp = uicontrol(d, 'Style', 'popup', 'String', str, ...
    'Position', [150 120 68 22], 'Value', val, 'BackGroundColor', 'w');

uicontrol('Parent', d, 'Style', 'text', 'Position', [230 116 90 22], ...
    'String', 'Missing value:', 'HorizontalAlignment', 'right');
h.extraV = uicontrol(d, 'Style', 'edit', 'String', num2str(pf.extraV), ...
    'Position', [324 120 60 22], 'BackGroundColor', 'w', ...
    'TooltipString', 'NaN or 0 is typical, but can be any number');

str = strtrim(cellstr(num2str([0 120 150 200 300 600 1200]')));
val = find(strcmp(str, pf.dpi));
uipanel(d, 'Units', 'Pixels', 'Position', [4 40 390 56], 'BorderType', 'etchedin', ...
    'BorderWidth', 2, 'Title', 'For "Save figure as" and "Copy figure"');
uicontrol('Parent', d, 'Style', 'text', 'Position', [8 46 90 22], ...
    'String', 'Resolution:', 'HorizontalAlignment', 'right');
h.dpi = uicontrol(d, 'Style', 'popup', 'String', str, ...
    'Position', [110 50 50 22], 'Value', val, 'BackGroundColor', 'w', ...
    'TooltipString', 'in DPI (0 means screen resolution)');
guidata(d, h);

%% Preference dialog callback
function pref_dialog_cb(h, ~, pref)
hs = guidata(h);
if strcmp(get(h, 'Tag'), 'OK') % done
    pf = get(pref, 'UserData');
    
    i = get(hs.interp, 'Value');
    str = get(hs.interp, 'String');
    pf.interp = str{i};
    
    pf.extraV = str2double(get(hs.extraV, 'String'));
    pf.openPath = get(hs.openPath, 'String');

    i = get(hs.dpi, 'Value');
    str = get(hs.dpi, 'String');
    pf.dpi = str{i}; 
    
    set(pref, 'UserData', pf);
    para = load(pf.pref_file); para = para.para;
    para.nii_viewer = pf;
    try save(pf.pref_file, 'para'), catch, end
    delete(get(h, 'Parent'));
elseif strcmp(get(h, 'Tag'), 'browse') % set openPath
    pth = uigetdir(pwd, 'Select folder for background image');
    if ~ischar(pth), return; end
    set(hs.openPath, 'String', pth);
end

%% Simple version of interp3
function V = interp3a(V, I, method)
if exist('griddedInterpolant', 'class')
    d = size(V); d(numel(d)+1:3) = 1;
    F = griddedInterpolant({1:d(1), 1:d(2), 1:d(3)}, V, method, 'none');
    V = F(I(1,:), I(2,:), I(3,:)); % interpolate
else % earlier matlab
    V = interp3(V, I(2,:), I(1,:), I(3,:), method, nan);
end

%% 2D/3D smooth wrapper: no input check for 2D
function out = smooth23(in, method, n, varargin)
if size(in,3)>1, out = smooth3(in, method, n, varargin{:}); return; end
if nargin<3 || isempty(n), n = 3; end
if numel(n)==1, n = [n n]; end
k = floor(n/2);
if k<1, out = in; return; end
n = k * 2 + 1; % odd number
if strcmp(method, 'box')
    kernal = ones(n) / n(1)/n(2);
elseif strcmp(method, 'gaussian')
    if nargin<4 || isempty(varargin{1}), sd = 0.65; 
    else sd = varargin{1}; 
    end
    [x, y] = meshgrid(-k(2):k(2), -k(1):k(1));
    kernal = exp(-(x.*x  +  y.*y) / (2*sd*sd));
    kernal = kernal / sum(kernal(:));
end
in = [repmat(in(1,:), [k(1) 1]); in; repmat(in(end,:), [k(1) 1])];
in = [repmat(in(:,1), [1 k(2)])  in  repmat(in(:,end), [1 k(2)])];
out = conv2(in, kernal, 'valid');

%% Show ijk/xyz and value
function set_xyz(hs, I)
if nargin<2, I = cell2mat(get(hs.spinner, 'Value'))'; end

I = round(I);
val = round(hs.q(1).R * [I-1 1]');
str = sprintf('(%g,%g,%g)=(%g,%g,%g): ', I, val(1:3));

p = get(hs.files, 'UserData');
for i = 1:numel(hs.q)
    if p(i).show == 0, continue; end
    t = round(p(i).volume);
    if hs.q(i).interp % need interpolation
        I0 = hs.q(i).R \ (hs.q(1).R * [I-1 1]'); % overlay ijk
        I0 = round(I0+1);
    else I0 = I;
    end
    try
        val = hs.q(i).nii.img(I0(1),I0(2),I0(3),t,:);
        val = val * hs.q(i).scl_slope + hs.q(i).scl_inter;
    catch
        val = nan; % out of range
    end
    fmtstr = '%.4g ';
    if numel(val)>1
        fmtstr = repmat(fmtstr, 1, numel(val));
        fmtstr = ['[' fmtstr]; fmtstr(end) = ']'; %#ok
    end
    str = sprintf(['%s ' fmtstr], str, val);
end
set(hs.xyz, 'String', str);
%%