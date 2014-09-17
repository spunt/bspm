function bspm_checkreg_gui(in, titles)
% BSPM_CHECKREG_BATCH Wrapper for Checking Registration
%
% USAGE: bspm_checkreg_batch(in)
%
% ARGUMENTS
%   in: an array of cells, with each cell containing paths for images
%   to loop over
%

% ------------------------ Copyright (C) 2014 ------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

% CHECK ARGS
if nargin<1, error('USAGE: bspm_checkreg_batch(in)'); end
if ~iscell(in), in = cellstr(in); end
% GUI
S.fig = spm_figure('Create','Graphics', 'Check Reg GUI', 'off');
set(S.fig, 'CloseRequestFcn', @cb_close, 'MenuBar', 'none'); 
if nargin<2 
    titles = strcat({'Group'}, cellstr(num2str((1:length(in))'))); 
end
checkreg(S.fig, in{1}, [0 0 0]', titles{1});
data = guihandles(S.fig);
data.titles = titles; 
data.in = in; 
data.count = 1;
data.xyz = [0 0 0]'; 
S.count = uimenu('Parent', S.fig, 'Label', sprintf('Group %d of %d', data.count, length(data.in))); 
guidata(S.fig, data);
pos.next        = [.900 .965 .09 .03];
pos.prev        = [.010 .965 .09 .03];
pos.save        = [.405 .965 .09 .03];
pos.which       = [.505 .965 .09 .03];

S.prev = uicontrol('Parent', S.fig, 'Units', 'Normal', 'FontUnits', 'Normal', ...
        'Style', 'Push', ...
        'Position', pos.prev, ...
        'FontUnits', 'norm', ...
        'FontName', 'Arial', ...
        'FontSize', .90, ...
        'FontWeight', 'Bold', ...  
        'String', '<<', ...
        'Enable', 'off', ...
        'Callback', {@cb_prev, S}); 
S.next = uicontrol('Parent', S.fig, 'Units', 'Normal', 'FontUnits', 'Normal', ...
        'Style', 'Push', ...
        'Position', pos.next, ...
        'FontUnits', 'norm', ...
        'FontName', 'Arial', ...
        'FontSize', .90, ...
        'FontWeight', 'Bold', ...  
        'String', '>>', ...
        'Enable', 'on', ...
        'Callback', {@cb_next, S});
S.which = uicontrol('Parent', S.fig, 'Units', 'Normal', 'FontUnits', 'Normal', ...
        'Style', 'Push', ...
        'Position', pos.which, ...
        'FontUnits', 'norm', ...
        'FontName', 'Arial', ...
        'FontSize', .55, ...
        'FontWeight', 'norm', ...  
        'String', 'Which', ...
        'Enable', 'on', ...
        'Callback', {@cb_which, S}); 
S.save = uicontrol('Parent', S.fig, 'Units', 'Normal', 'FontUnits', 'Normal', ...
        'Style', 'Push', ...
        'Position', pos.save, ...
        'FontUnits', 'norm', ...
        'FontName', 'Arial', ...
        'FontSize', .55, ...
        'FontWeight', 'norm', ...  
        'String', 'Save', ...
        'Enable', 'on', ...
        'Callback', {@cb_save, S}); 
set(S.fig, 'Visible', 'on');
end
%% CALLBACKS
function cb_next(varargin)
    action = varargin{1}; 
    parent = varargin{3};
    data = guidata(parent.fig);
    data.xyz = spm_orthviews('pos'); 
    data.count = data.count + 1; 
    set(parent.count, 'Label', sprintf('Group %d of %d', data.count, length(data.in))); 
    if data.count==length(data.in)
        set(action, 'Enable', 'off'); 
    end
    set(parent.fig, 'Visible', 'off');
    if data.count==2
        set(parent.prev, 'Enable', 'on');
    end
    checkreg(parent.fig, data.in{data.count}, data.xyz, data.titles{data.count});
    set(parent.fig, 'Visible', 'on');
    guidata(parent.fig, data);
end
function cb_prev(varargin)
    action = varargin{1}; 
    parent = varargin{3};
    data = guidata(parent.fig);
    data.xyz = spm_orthviews('pos'); 
    data.count = data.count - 1; 
    set(parent.count, 'Label', sprintf('Group %d of %d', data.count, length(data.in))); 
    if data.count==1
        set(action, 'Enable', 'off'); 
    end
    set(parent.fig, 'Visible', 'off');
    checkreg(parent.fig, data.in{data.count}, data.xyz, data.titles{data.count});
    set(parent.fig, 'Visible', 'on');
    guidata(parent.fig, data);
end
function cb_save(varargin)
    parent = varargin{3};
    data = guidata(parent.fig);
    outname = fullfile(pwd, sprintf('CheckReg_%s_x=%d_y=%d_z=%d.pdf', regexprep(data.titles{data.count}, ' ', '_'), round(data.xyz(:)))); 
    saveas(parent.fig, outname, 'pdf');
    disp(['Saved to: ' outname]); 
end
function cb_which(varargin)
    parent = varargin{3};
    data = guidata(parent.fig);
    gui_cellshow(data.in{data.count}, data.titles{data.count});
end
function cb_close(varargin)
    if length(varargin)==3
        parent = varargin{3};
        h = parent.fig;
    else
        h = varargin{1};
    end
    delete(h); % Bye-bye figure
end
%% SUBFUNCTIONS
function gui_cellshow(in, title)
% CELLSHOW Show cell array
%
%       USAGE: h = cellshow(in, [title], [keeppath], [hideoutput])
%           
% -------------------------------------------------------------------
if nargin<2, title = 'File List'; end;
if nargin<1, error('USAGE: h = cellshow(in, [title])'); end
if ischar(in), in = cellstr(in); end
if iscell(title), title = char(title); end
N = cell(length(in),1);
N(:,1) = num2cell(1:length(in));
str = sprintf('%%0%dd', max(cellfun('length', cellfun(@num2str, N, 'Unif', false))));
for i = 1:length(in), N{i} = sprintf(str, i); end
in = strcat({' - '}, in);
in = strcat(cellstr(N), in);
%% GUI Figure
h.fh = figure('Units','Normal',...
        'Color',[0 0 0],...
        'MenuBar','none',...
        'Name', title,...
        'NumberTitle','off',...
        'Position',[.25 .15 .50 .25],... % left bottom width height
        'Resize','on');
%% List
h.list = uicontrol('Parent',h.fh,...
        'FontUnits','Normal',...
        'Style','list',...
        'Units','Normal',...
        'Enable', 'on', ...
        'Position',[0 0 1 1],... % left bottom width height
        'String', in,...
        'FontName','Arial',...
        'BackgroundColor', [1 1 1], ...
        'ForegroundColor',[0 0 0],...
        'HorizontalAlignment','left',...
        'FontUnits', 'points', ...
        'FontWeight', 'bold', ...
        'FontSize', 20); 
end
function checkreg(F, images, xyz, title)
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
    if ischar(images), images = cellstr(images); end
    captions = cellfun(@fileparts, images, 'Unif', false); 
    [~,captions] = cellfun(@fileparts, captions, 'Unif', false);
    if iscell(images{1}), images = spm_vol(vertcat(images{:})); end 
    chkeep  = [findobj(F, 'type', 'uicontrol'); findobj(F, 'type', 'uimenu')]; 
    chall   = get(F, 'children'); 
    delete(chall(~ismember(chall, chkeep))); 
    spm_orthviews('Reset');
    mn = length(images);
    n  = round(mn^0.4);
    m  = ceil(mn/n);
    w  = 1/n;
    h  = 1/m;
    ds = (w+h)*0.02;
    if nargin==4
        if mn==2, 
            set(get(gca,'Title'), 'String', title, 'FontSize', 24, 'FontName', 'Arial', 'FontWeight', 'Bold', 'Position', [0.6 0.5]); 
        elseif mn>2,
            set(get(gca,'Title'), 'String', title, 'FontSize', 24, 'FontName', 'Arial', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center'); 
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
                spm_orthviews('Caption', ij, captions{ij}, 'FontSize', 11, 'FontName', 'Arial');
            end
        end
    end
    spm_orthviews('Reposition', xyz); 
    set(findobj(gcf, 'type', 'line'), 'color',[0.7020    0.8039    0.8902], 'linewidth', 1);  
end

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
