function varargout=lui_spm(varargin)
% lui_spm: Statistical Parametric Mapping (startup function)
%_______________________________________________________________________
%  ___  ____  __  __
% / __)(  _ \(  \/  )  
% \__ \ )___/ )    (   Statistical Parametric Mapping
% (___/(__)  (_/\/\_)  lui_spm - http://www.fil.ion.ucl.ac.uk/lui_spm/
%_______________________________________________________________________
%
% lui_spm (Statistical Parametric Mapping) is a package for the analysis
% functional brain mapping experiments. It is the in-house package of
% the Wellcome Department of Cognitive Neurology, and is available to
% the scientific community as copyright freeware under the terms of the
% GNU General Public Licence.
% 
% Theoretical, computational and other details of the package are
% available in lui_spm's "Help" facility. This can be launched from the
% main lui_spm Menu window using the "Help" button, or directly from the
% command line using the command `lui_spm_help`.
%
% Details of this release are available via the "About lui_spm" help topic
% (file lui_spm.man), accessible from the lui_spm splash screen.  (Or type
% `lui_spm_help lui_spm.man` in the MatLab command window)
% 
% This lui_spm function initialises the default parameters, and displays a
% splash screen with buttons leading to the PET(SPECT) & fMRI
% modalities Alternatively, `lui_spm('pet')` and `lui_spm('fmri')`
% (equivalently `lui_spm pet` and `lui_spm mri`) lead directly to the respective
% modality interfaces.
%
% Once the modality is chosen, (and it can be toggled mid-session) the
% lui_spm user interface is displayed. This provides a constant visual
% environment in which data analysis is implemented. The layout has
% been designed to be simple and at the same time show all the
% facilities that are available. The interface consists of three
% windows: A menu window with pushbuttons for the lui_spm routines (each
% button has a 'CallBack' string which launches the appropriate
% function/script); A blank panel used for interaction with the user;
% And a graphics figure with various editing and print facilities (see
% lui_spm_figure.m). (These windows are 'Tag'ged 'Menu', 'Interactive', and
% 'Graphics' respectively, and should be referred to by their tags
% rather than their figure numbers.)
%
% Further interaction with the user is (mainly) via questioning in the
% 'Interactive' window (managed by lui_spm_input), and file selection
% (managed by lui_spm_select). See the help on lui_spm_input.m and lui_spm_select.m for
% details on using these functions.
%
% If a "message of the day" file named lui_spm_motd.man exists in the lui_spm
% directory (alongside lui_spm.m) then it is displayed in the Graphics
% window on startup.
%
% Arguments to this routine (lui_spm.m) lead to various setup facilities,
% mainly of use to lui_spm power users and programmers. See programmers
% FORMAT & help in the main body of lui_spm.m
%
%_______________________________________________________________________
% lui_spm is developed by members and collaborators of the
% Wellcome Department of Imaging Neuroscience

%-SVN ID and authorship of this program...
%-----------------------------------------------------------------------
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Andrew Holmes
% $Id: lui_spm.m 1019 2007-12-04 17:46:25Z guillaume $


%=======================================================================
% - FORMAT specifications for embedded CallBack functions
%=======================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take. Recall   )
%( MatLab's command-function duality: `lui_spm Welcome` is equivalent to   )
%( `lui_spm('Welcome')`.                                                   )
%
% FORMAT lui_spm
% Defaults to lui_spm('Welcome')
%
% FORMAT lui_spm('Welcome')
% Clears command window, deletes all figures, prints welcome banner and
% splash screen, sets window defaults.
%
% FORMAT lui_spm('AsciiWelcome')
% Prints ASCII welcome banner in MatLab command window.
%
% FORMAT lui_spm('PET') lui_spm('FMRI')
% Closes all windows and draws new Menu, Interactive, and Graphics
% windows for an lui_spm session. The buttons in the Menu window launch the
% main analysis routines.
%
% FORMAT Fmenu = lui_spm('CreateMenuWin',Vis)
% Creates lui_spm menu window, 'Tag'ged 'Menu'
% F   - handle of figure created
% Vis - Visibility, 'on' or 'off'
%
% Finter = FORMAT lui_spm('CreateIntWin',Vis)
% Creates an lui_spm Interactive window, 'Tag'ged 'Interactive'
% F   - handle of figure created
% Vis - Visibility, 'on' or 'off'
%
% FORMAT lui_spm('ChMod',Modality)
% Changes modality of lui_spm: Currently lui_spm supports PET & MRI modalities,
% each of which have a slightly different Menu window and different
% defaults. This function switches to the specified modality, setting
% defaults and displaying the relevant buttons.
%
% FORMAT lui_spm('defaults',Modality)
% Sets default global variables for the specified modality.
%
% FORMAT [Modality,ModNum]=lui_spm('CheckModality',Modality)
% Checks the specified modality against those supported, returns
% upper(Modality) and the Modality number, it's position in the list of
% supported Modalities.
%
% FORMAT WS=lui_spm('WinScale')
% Returns ratios of current display dimensions to that of a 1152 x 900
% Sun display. WS=[Xratio,Yratio,Xratio,Yratio]. Used for scaling other
% GUI elements.
% (Function duplicated in lui_spm_figure.m, repeated to reduce inter-dependencies.)
%
% FORMAT [FS,sf] = lui_spm('FontSize',FS)
% FORMAT [FS,sf] = lui_spm('FontSizes',FS)
% Returns fontsizes FS scaled for the current display.
% FORMAT sf = lui_spm('FontScale')
% Returns font scaling factor
% FS     - (vector of) Font sizes to scale [default [1:36]]
% sf     - font scaling factor (FS(out) = floor(FS(in)*sf)
%
% Rect = lui_spm('WinSize',Win,raw)
% Returns sizes and positions for lui_spm windows.
% Win  - 'Menu', 'Interactive', 'Graphics', or '0'
%      -  Window whose position is required. Only first character is
%         examined. '0' returns size of root workspace.
% raw  - If specified, then positions are for a 1152 x 900 Sun display.
%        Otherwise the positions are scaled for the current display.
%
% FORMAT lui_spmdir=lui_spm('Dir',Mfile)
% Returns the directory containing the version of lui_spm in use,
% identified as the first in MATLABPATH containing the Mfile lui_spm (this
% file) (or Mfile if specified).
%
% FORMAT [v,c]=lui_spm('Ver',Mfile,ReDo,Cache,Con)
% Returns the current version (v) & copyright notice, extracted from
% the top line of the Contents.m file in the directory containing the
% currently used file Mfile (defaults on omission or empty to 'lui_spm').
%
%-The version and copyright information are saved in a global
% variable called [upper(lui_spm_str_manip(Mfile,'rt')),'_VER'], as a
% structure with fields 'v' and 'c'. This enables repeat use without
% recomputation.
%
%-If Con [default (missing or empty) 1] is false, then the version
% information is extracted from Mfile itself, rather than the
% Contents.m file in the same directory. When using a Contents.m file,
% the first line is read. For other files, the second line (the H1 help
% line) is used. This is for consistency with MatLab's ver and help
% commands respectively. (This functionality enables toolboxes to be
% identified by a function rather than a Contents.m file, allowing
% installation in a directory which already has a Contents.m file.)
%
%-If Cache [default (missing or empty) 1] is true, then the version and
% copyright information cached in the global variable
% [upper(Mfile),'_VER'], as a structure with fields 'v' and 'c'. This
% enables repeat use without recomputation.
%
%-If ReDo [default (missing or empty) 0] is true, then the version and
% copyright information are recomputed (regardless of any stored global
% data).
%
% FORMAT xTB = lui_spm('TBs')
% Identifies installed lui_spm toolboxes: lui_spm toolboxes are defined as the
% contents of sub-directories of fullfile(lui_spm('Dir'),'toolbox') - the
% lui_spm toolbox installation directory. For lui_spm to pick a toolbox up,
% there must be a single mfile in the directory whose name ends with
% the toolbox directory name. (I.e. A toolbox called "test" would be in
% the "test" subdirectory of lui_spm('Dir'), with a single file named
% *test.m.) This M-file is regarded as the launch file for the
% toolbox.
% xTB      - structure array containing toolbox definitions
% xTB.name - name of toolbox (taken as toolbox directory name)
% xTB.prog - launch program for toolbox
% xTB.dir  - toolbox directory
%
% FORMAT lui_spm('TBlaunch',xTB,i)
% Launch a toolbox, prepending TBdir to path if necessary
% xTB      - toolbox definition structure (i.e. from lui_spm('TBs')
% xTB.name - name of toolbox
% xTB.prog - name of program to launch toolbox
% xTB.dir  - toolbox directory (prepended to path if not on path)
%
% FORMAT [c,cName] = lui_spm('Colour')
% Returns the RGB triple and a description for the current en-vogue lui_spm
% colour, the background colour for the Menu and Help windows.
%
% FORMAT [v1,v2,...] = lui_spm('GetGlobal',name1,name2,...)
% Returns values of global variables (without declaring them global)
% name1, name2,... - name strings of desired globals
% a1, a2,...       - corresponding values of global variables with given names
%                    ([] is returned as value if global variable doesn't exist)
%
% FORMAT CmdLine = lui_spm('CmdLine',CmdLine)
% Command line lui_spm usage?
% CmdLine (input)  - CmdLine preference
%                    [defaults (missing or empty) to global defaults.cmdline,]
%                    [if it exists, or 0 (GUI) otherwise.                    ]
% CmdLine (output) - true if global CmdLine if true,
%                    or if on a terminal with no support for graphics windows.
%
% FORMAT v = lui_spm('MLver')
% Returns MatLab version, truncated to major & minor revision numbers
%
% FORMAT lui_spm('SetCmdWinLabel',WinStripe,IconLabel)
% Sets the names on the headers and icons of Sun command tools.
% WinStripe defaults to a summary line identifying the user, host and
% MatLab version; IconLabel to 'MatLab'.
%
% FORMAT lui_spm('PopUpCB',h)
% Callback handler for PopUp UI menus with multiple callbacks as cellstr UserData
%
% FORMAT str = lui_spm('GetUser',fmt)
% Returns current users login name, extracted from the hosting environment
% fmt   - format string: If USER is defined then sprintf(fmt,USER) is returned
%
% FORMAT lui_spm('Beep')
% plays the keyboard beep!
%
% FORMAT lui_spm('time')
% Returns the current time and date as hh:mm dd/mm/yyyy
%
% FORMAT lui_spm('Pointer',Pointer)
% Changes pointer on all lui_spm (HandleVisible) windows to type Pointer
% Pointer defaults to 'Arrow'. Robust to absence of windows
%
% FORMAT h = lui_spm('alert',Message,Title,CmdLine,wait)
% FORMAT h = lui_spm('alert"',Message,Title,CmdLine,wait)
% FORMAT h = lui_spm('alert*',Message,Title,CmdLine,wait)
% FORMAT h = lui_spm('alert!',Message,Title,CmdLine,wait)
% Displays an alert, either in a GUI msgbox, or as text in the command window.
%  ( 'alert"' uses the 'help' msgbox icon, 'alert*' the )
%  ( 'error' icon, 'alert!' the 'warn' icon             )
% Message - string (or cellstr) containing message to print
% Title   - title string for alert
% CmdLine - CmdLine preference [default lui_spm('CmdLine')]
%         - If CmdLine is complex, then a CmdLine alert is always used,
%           possibly in addition to a msgbox (the latter according
%           to lui_spm('CmdLine').)
% wait    - if true, waits until user dismisses GUI / confirms text alert
%           [default 0] (if doing both GUI & text, waits on GUI alert)
% h       - handle of msgbox created, empty if CmdLine used
%
% FORMAT lui_spmid = lui_spm('FnBanner', Fn,FnV)
% Prints a function start banner, for version FnV of function Fn, & datestamps
% FORMAT lui_spmid = lui_spm('SFnBanner',Fn,FnV)
% Prints a sub-function start banner
% FORMAT lui_spmid = lui_spm('SSFnBanner',Fn,FnV)
% Prints a sub-sub-function start banner
% Fn    - Function name (string)
% FnV   - Function version (string)
% lui_spmid - ID string: [lui_spmver: Fn (FnV)] 
%
% FORMAT [Finter,Fgraph,CmdLine] = lui_spm('FnUIsetup',Iname,bGX,CmdLine)
% Robust UIsetup procedure for functions:
%   Returns handles of 'Interactive' and 'Graphics' figures.
%   Creates 'Interactive' figure if ~CmdLine, creates 'Graphics' figure if bGX.
% Iname   - Name for 'Interactive' window
% bGX     - Need a Graphics window? [default 1]
% CmdLine - CommandLine usage? [default lui_spm('CmdLine')]
% Finter  - handle of 'Interactive' figure
% Fgraph  - handle of 'Graphics' figure
% CmdLine - CommandLine usage?
%
% FORMAT F = lui_spm('FigName',Iname,F,CmdLine)
% Set name of figure F to "lui_spmver (User): Iname" if ~CmdLine
% Robust to absence of figure.
% Iname      - Name for figure
% F (input)  - Handle (or 'Tag') of figure to name [default 'Interactive']
% CmdLine    - CommandLine usage? [default lui_spm('CmdLine')]
% F (output) - Handle of figure named
%
% FORMAT lui_spm('GUI_FileDelete')
% CallBack for GUI for file deletion, using lui_spm_select and confirmation dialogs
%
% FORMAT Fs = lui_spm('Show')
% Opens all lui_spm figure windows (with HandleVisibility) using `figure`.
%   Maintains current figure.
% Fs - vector containing all HandleVisible figures (i.e. get(0,'Children'))
%
% FORMAT lui_spm('Clear',Finter, Fgraph)
% Clears and resets lui_spm-GUI, clears and timestamps MatLab command window
% Finter  - handle or 'Tag' of 'Interactive' figure [default 'Interactive']
% Fgraph  - handle or 'Tag' of 'Graphics' figure [default 'Graphics']
%
% FORMAT lui_spm('Help',varargin)
% Merely a gateway to lui_spm_help(varargin) - so you can type "lui_spm help"
% 
%_______________________________________________________________________

%-Disable warning messages due to dll files still existing
%-----------------------------------------------------------------------

%-Disable Java if necessary
%-----------------------------------------------------------------------
try
    feature('JavaFigures',0);
end

%-Parameters
%-----------------------------------------------------------------------
Modalities = {'PET','FMRI','EEG'};

%-Format arguments
%-----------------------------------------------------------------------
if nargin == 0, Action='Welcome'; else, Action = varargin{1}; end


%=======================================================================
switch lower(Action), case 'welcome'             %-Welcome splash screen
%=======================================================================
check_installation
lui_spm_defaults;
global defaults
if isfield(defaults,'modality'), lui_spm(defaults.modality); return; end;

%-Open startup window, set window defaults
%-----------------------------------------------------------------------
openfig('lui_spm_Welcome');


%=======================================================================
case 'asciiwelcome'                           %-ASCII lui_spm banner welcome
%=======================================================================
disp( ' ___  ____  __  __                                            ');
disp( '/ __)(  _ \(  \/  )                                           ');
disp( '\__ \ )___/ )    (   Statistical Parametric Mapping           ');
disp(['(___/(__)  (_/\/\_)  ',lui_spm('Ver'),' - http://www.fil.ion.ucl.ac.uk/lui_spm/']);
fprintf('\n');


%=======================================================================
case lower(Modalities)            %-Initialise lui_spm in PET, fMRI modality
%=======================================================================
% lui_spm(Modality)
check_installation

%-Initialisation and workspace canonicalisation
%-----------------------------------------------------------------------
local_clc, lui_spm('SetCmdWinLabel')
lui_spm('AsciiWelcome');                    fprintf('\n\nInitialising lui_spm');
Modality = upper(Action);                                  fprintf('.');
delete(get(0,'Children'));                                 fprintf('.');

%-Draw lui_spm windows
%-----------------------------------------------------------------------
Fmenu  = lui_spm('CreateMenuWin','off');fprintf('.');
Finter = lui_spm('CreateIntWin','off');	                       fprintf('.');
Fgraph = lui_spm_figure('Create','Graphics','Graphics','off'); fprintf('.');
   
lui_spm_figure('WaterMark',Finter,lui_spm('Ver'),'',45);           fprintf('.');

Fmotd  = fullfile(lui_spm('Dir'),'lui_spm_motd.man');
if exist(Fmotd), lui_spm_help('!Disp',Fmotd,'',Fgraph,lui_spm('Ver')); end
                                                           fprintf('.');

%-Load startup global defaults
%-----------------------------------------------------------------------
lui_spm_defaults;                                              fprintf('.');

%-Setup for current modality
%-----------------------------------------------------------------------
lui_spm('ChMod',Modality);                                     fprintf('.');

%-Reveal windows
%-----------------------------------------------------------------------
set([Fmenu,Finter,Fgraph],'Visible','on');          fprintf('done\n\n');

%-Print present working directory
%-----------------------------------------------------------------------
fprintf('lui_spm present working directory:\n\t%s\n',pwd)


%=======================================================================
case 'createmenuwin'                              %-Draw lui_spm menu window
%=======================================================================
% Fmenu = lui_spm('CreateMenuWin',Vis)

%-Close any existing 'Menu' 'Tag'ged windows
%-----------------------------------------------------------------------
delete(lui_spm_figure('FindWin','Menu'))
Fmenu     = openfig('lui_spm_Menu');

%-Patch up the callback to realign_unwarp
tmp = findobj(Fmenu,'ToolTipString','realignment','Tag','FMRI');
set(tmp,'Callback',@realign_unwarp);

%-Set lui_spm colour
%-----------------------------------------------------------------------
set(findobj(Fmenu,'Tag', 'frame'),'backgroundColor',lui_spm('colour'));

%-Set toolbox
%-----------------------------------------------------------------------
xTB       = lui_spm('tbs');
set(findobj(Fmenu,'Tag', 'Toolbox'),'String',{'Toolbox:' xTB.name });
set(findobj(Fmenu,'Tag', 'Toolbox'),'UserData',xTB);
varargout = {Fmenu};
return

%=======================================================================
case 'createintwin'                      %-Create lui_spm interactive window
%=======================================================================
% Finter = lui_spm('CreateIntWin',Vis)
%-----------------------------------------------------------------------
delete(lui_spm_figure('FindWin','Interactive'))
Finter    = openfig('lui_spm_Interactive');
varargout = {Finter};
return


%=======================================================================
case 'chmod'                            %-Change lui_spm modality PET<->fMRI
%=======================================================================
% lui_spm('ChMod',Modality)
%-----------------------------------------------------------------------

%-Sort out arguments
%-----------------------------------------------------------------------
if nargin<2, Modality = ''; else, Modality = varargin{2}; end
[Modality,ModNum] = lui_spm('CheckModality',Modality);

%-Sort out global defaults
%-----------------------------------------------------------------------
lui_spm('defaults',Modality);

%-Sort out visability of appropriate controls on Menu window
%-----------------------------------------------------------------------
Fmenu = lui_spm_figure('FindWin','Menu');
if isempty(Fmenu), error('lui_spm Menu window not found'), end

if strcmpi(Modality,'PET')
	set(findobj(Fmenu, 'Tag', 'FMRI'),    'Visible', 'off');
	set(findobj(Fmenu, 'Tag', 'EEG'),     'Visible', 'off');
	set(findobj(Fmenu, 'Tag', 'PETFMRI'), 'Visible', 'on' );
	set(findobj(Fmenu, 'Tag', 'PET'),     'Visible', 'on' );
elseif strcmpi(Modality,'FMRI')
	set(findobj(Fmenu, 'Tag', 'EEG'),     'Visible', 'off');
	set(findobj(Fmenu, 'Tag', 'PET'),     'Visible', 'off');
	set(findobj(Fmenu, 'Tag', 'PETFMRI'), 'Visible', 'on' );
	set(findobj(Fmenu, 'Tag', 'FMRI'),    'Visible', 'on' );
else
	set(findobj(Fmenu, 'Tag', 'PETFMRI'), 'Visible', 'off');
	set(findobj(Fmenu, 'Tag', 'PET'),     'Visible', 'off');
	set(findobj(Fmenu, 'Tag', 'FMRI'),    'Visible', 'off');
	set(findobj(Fmenu, 'Tag', 'EEG'),     'Visible', 'on' );
end
set(findobj(Fmenu,'Tag','Modality'),'Value',ModNum,'UserData',ModNum);
lui_spm_jobman('chmod',Modality);

%=======================================================================
case 'defaults'                 %-Set lui_spm defaults (as global variables)
%=======================================================================
% lui_spm('defaults',Modality)
%-----------------------------------------------------------------------
global defaults
if isempty(defaults), lui_spm_defaults; end;

%-Sort out arguments
%-----------------------------------------------------------------------
if nargin<2, Modality=''; else, Modality=varargin{2}; end
Modality          = lui_spm('CheckModality',Modality);
defaults.modality = Modality;
defaults.SWD      = lui_spm('Dir');              % lui_spm directory
defaults.TWD      = lui_spm_platform('tempdir'); % Temp directory
	
%-Set Modality specific default (global variables)
%-----------------------------------------------------------------------
global UFp
if strcmpi(defaults.modality,'pet')
	UFp	= defaults.stats.pet.ufp;		% Upper tail F-prob
elseif strcmpi(defaults.modality,'fmri')
	UFp	= defaults.stats.fmri.ufp;		% Upper tail F-prob
elseif strcmpi(defaults.modality,'eeg')
    ;
elseif strcmpi(defaults.modality,'unknown')
else
	error('Illegal Modality');
end


%=======================================================================
case 'quit'                                      %-Quit lui_spm and clean up
%=======================================================================
% lui_spm('Quit')
%-----------------------------------------------------------------------
delete(get(0,'Children'));
local_clc;
fprintf('Bye for now...\n\n');


%=======================================================================
case 'checkmodality'              %-Check & canonicalise modality string
%=======================================================================
% [Modality,ModNum] = lui_spm('CheckModality',Modality)
%-----------------------------------------------------------------------
if nargin<2, Modality=''; else, Modality=upper(varargin{2}); end
if isempty(Modality)
	global defaults
	if isfield(defaults,'modality'), Modality = defaults.modality;
	else, Modality = 'UNKNOWN'; end
end
if ischar(Modality)
	ModNum = find(ismember(Modalities,Modality));
else
	if ~any(Modality == [1:length(Modalities)])
		Modality = 'ERROR';
		ModNum   = [];
	else
		ModNum   = Modality;
		Modality = Modalities{ModNum};
	end
end

if isempty(ModNum), error('Unknown Modality'), end
varargout = {upper(Modality),ModNum};


%=======================================================================
case {'winscale','getwinscale'}  %-Window scale factors (to fit display)
%=======================================================================
% WS = lui_spm('WinScale')
%-----------------------------------------------------------------------
if strcmp(lower(Action),'getwinscale')
	warning('lui_spm(''GetWinScale'' GrandFathered, use ''WinScale''')
end
if lui_spm_matlab_version_chk('7') >=0
	S0 = get(0, 'MonitorPosition');
	S0 = S0(1,:);
else
	S0   = get(0,'ScreenSize');
end;
if all(S0==1), error('Can''t open any graphics windows...'), end

tmp = [S0(3)/1152 (S0(4)-50)/900 S0(3)/1152 (S0(4)-50)/900];
varargout = {min(tmp)*[1 1 1 1]};

% Make sure that aspect ratio is about right - for funny shaped screens
% varargout = {[S0(3)/1152 (S0(4)-50)/900 S0(3)/1152 (S0(4)-50)/900]};


%=======================================================================
case {'fontsize','fontsizes','fontscale'}                 %-Font scaling
%=======================================================================
% [FS,sf] = lui_spm('FontSize',FS)
% [FS,sf] = lui_spm('FontSizes',FS)
% sf = lui_spm('FontScale')
%-----------------------------------------------------------------------
if nargin<3, c=0; else, c=1; end
if nargin<2, FS=[1:36]; else, FS=varargin{2}; end

sf  = 1 + 0.85*(min(lui_spm('WinScale'))-1);

if strcmp(lower(Action),'fontscale')
	varargout = {sf};
else
	varargout = {ceil(FS*sf),sf};
end


%=======================================================================
case 'winsize'                 %-Standard lui_spm window locations and sizes
%=======================================================================
% Rect = lui_spm('WinSize',Win,raw)
%-----------------------------------------------------------------------
if nargin<3, raw=0; else, raw=1; end
if nargin<2, Win=''; else, Win=varargin{2}; end

Rect = [	[108 466 400 445];...
		[108 045 400 395];...
		[515 015 600 865] ];

WS = lui_spm('WinScale');

if isempty(Win)
	WS = ones(3,1)*WS;
elseif upper(Win(1))=='M'
	%-Menu window
	Rect = Rect(1,:);
elseif upper(Win(1))=='I'
	%-Interactive window
	Rect = Rect(2,:);
elseif upper(Win(1))=='G'
	%-Graphics window
	Rect = Rect(3,:);
elseif Win(1)=='0'
	%-Root workspace
if lui_spm_matlab_version_chk('7') >=0
		Rect = get(0, 'MonitorPosition');
		Rect = Rect(1,:);
	else
		Rect = get(0,'ScreenSize');
	end;
else
	error('Unknown Win type');
end

if ~raw, Rect = Rect.*WS; end
varargout = {Rect};


%=======================================================================
case 'dir'                           %-Identify specific (lui_spm) directory
%=======================================================================
% lui_spm('Dir',Mfile)
%-----------------------------------------------------------------------
if nargin<2, Mfile='lui_spm'; else, Mfile=varargin{2}; end
lui_spmdir = which(Mfile);
if isempty(lui_spmdir)			%-Not found or full pathname given
	if exist(Mfile,'file')==2	%-Full pathname
		lui_spmdir = Mfile;
	else
		error(['Can''t find ',Mfile,' on MATLABPATH']);
	end
end
[lui_spmdir,junk] = fileparts(lui_spmdir);

if exist('isdeployed') && isdeployed,
    ind = findstr(lui_spmdir,'_mcr')-1;
    [lui_spmdir,junk] = fileparts(lui_spmdir(1:ind(1)));
end;
varargout = {lui_spmdir};


%=======================================================================
case 'ver'                                                 %-lui_spm version
%=======================================================================
% lui_spmver = lui_spm('Ver',Mfile,ReDo,Cache,Con)
%-----------------------------------------------------------------------
if nargin<5, Con=[]; else, Con=varargin{5}; end
if isempty(Con), Con=1; end
if nargin<4, Cache=[]; else, Cache=varargin{4}; end
if isempty(Cache), Cache=1; end
if nargin<3, ReDo=[]; else, ReDo=varargin{3}; end
if isempty(ReDo), ReDo=0; end
if nargin<2, Mfile=''; else, Mfile=varargin{2}; end
if isempty(Mfile), Mfile='lui_spm'; end

xVname = [upper(lui_spm_str_manip(Mfile,'rt')),'_VER'];

%-See if version info exists in global variable
%-----------------------------------------------------------------------
xV = lui_spm('GetGlobal',xVname);
if ~ReDo & ~isempty(xV)
	if isstruct(xV) & isfield(xV,'v') & isfield(xV,'c')
		varargout = {xV.v,xV.c};
		return
	end
end

%-Work version out from file
%-----------------------------------------------------------------------
if Con
	Vfile = fullfile(lui_spm('Dir',Mfile),'Contents.m');
	skip = 0;	%-Don't skip first line
else
	Vfile = which(Mfile);
	if isempty(Vfile), error(['Can''t find ',Mfile,' on MATLABPATH']); end
	skip = 1;	%-Skip first line
end
if exist(Vfile)
	fid = fopen(Vfile,'r');
	str = fgets(fid);
	if skip, str=fgets(fid); end
	fclose(fid);
	str(1:max(1,min(find(str~='%' & str~=' '))-1))=[];
	tmp = min(find(str==10|str==32));
	v = str(1:tmp-1);
	if str(tmp)==32
		c = str(tmp+1:tmp+min(find(str(tmp+1:end)==10))-1);
	else
		c = '(c) Copyright reserved';
	end
else
	v = 'lui_spm';
	c = '(c) Copyright reserved';
end

%-Store version info in global variable
%-----------------------------------------------------------------------
if Cache
	eval(['global ',xVname])
	eval([xVname,' = struct(''v'',v,''c'',c);'])
end

varargout = {v,c};


%=======================================================================
case 'tbs'                                %-Identify installed toolboxes
%=======================================================================
% xTB = lui_spm('TBs')
%-----------------------------------------------------------------------

% Toolbox directory
%-----------------------------------------------------------------------
Tdir  = fullfile(lui_spm('Dir'),'toolbox');

%-List of potential installed toolboxes directories
%-----------------------------------------------------------------------
if exist(Tdir,'dir')
	d = dir(Tdir); 
	d = {d([d.isdir]).name};
	d = {d{cellfun('isempty',regexp(d,'^\.'))}};
else
	d = {};
end


%-Look for a "main" M-file in each potential directory
%-----------------------------------------------------------------------
xTB = [];
for i = 1:length(d)
    tdir = fullfile(Tdir,d{i});
    fn   = cellstr(lui_spm_select('List',tdir,['^.*' d{i} '\.m$']));

    if ~isempty(fn{1}),
        xTB(end+1).name = strrep(d{i},'_','');
        xTB(end).prog   = lui_spm_str_manip(fn{1},'r');
        xTB(end).dir    = tdir;
    end;

end

varargout{1} = xTB;


%=======================================================================
case 'tblaunch'                                  %-Launch an lui_spm toolbox
%=======================================================================
% xTB = lui_spm('TBlaunch',xTB,i)
%-----------------------------------------------------------------------
if nargin < 3, i   = 1;          else i   = varargin{3}; end
if nargin < 2, xTB = lui_spm('TBs'); else xTB = varargin{2}; end

if i > 0
	%-Addpath (& report)
	%-------------------------------------------------------------------
	if isempty(findstr(xTB(i).dir,path))
		addpath(xTB(i).dir,'-begin');
		lui_spm('alert"',{'Toolbox directory prepended to Matlab path:',...
			xTB(i).dir},...
			[xTB(i).name,' toolbox'],1);
	end

	%-Launch
	%-------------------------------------------------------------------
	evalin('base',xTB(i).prog);
end


%=======================================================================
case 'colour'                                     %-lui_spm interface colour
%=======================================================================
% lui_spm('Colour')
%-----------------------------------------------------------------------
%-Pre-developmental livery
% varargout = {[1.0,0.2,0.3],'fightening red'};
%-Developmental livery
% varargout = {[0.7,1.0,0.7],'flourescent green'};
%-Alpha release livery
% varargout = {[0.9,0.9,0.5],'over-ripe banana'};
%-Beta release livery
  varargout = {[0.9 0.8 0.9],'blackcurrant purple'};
%-Distribution livery
% varargout = {[0.8 0.8 1.0],'vile violet'};;

global defaults
if isempty(defaults), lui_spm_defaults; end;
if isfield(defaults,'ui') && isfield(defaults.ui,'colour2'),
	varargout{1} = defaults.ui.colour2;
end;

%=======================================================================
case 'getglobal'                           %-Get global variable cleanly
%=======================================================================
% varargout = lui_spm('GetGlobal',varargin)
%-----------------------------------------------------------------------
wg = who('global');
for i=1:nargin-1
	if any(strcmp(wg,varargin{i+1}))
		eval(['global ',varargin{i+1},', tmp=',varargin{i+1},';'])
		varargout{i} = tmp;
	else
		varargout{i} = [];
	end
end

%=======================================================================
case {'cmdline','isgcmdline'}                   %-lui_spm command line mode?
%=======================================================================
% CmdLine = lui_spm('CmdLine',CmdLine)
% isGCmdLine usage is Grandfathered
%-----------------------------------------------------------------------
if nargin<2, CmdLine=[]; else, CmdLine = varargin{2}; end
if isempty(CmdLine),
	global defaults
	if ~isempty(defaults) & isfield(defaults,'cmdline'),
		CmdLine = defaults.cmdline;
	else,
		CmdLine = 0;
	end;
end
varargout = {CmdLine * (get(0,'ScreenDepth')>0)};

%=======================================================================
case 'mlver'                       %-MatLab major & point version number
%=======================================================================
% v = lui_spm('MLver')
%-----------------------------------------------------------------------
v = version; tmp = find(v=='.');
if length(tmp)>1, varargout={v(1:tmp(2)-1)}; end

%=======================================================================
case 'setcmdwinlabel'      %-Set command window label (Sun OpenWin only)
%=======================================================================
% lui_spm('SetCmdWinLabel',WinStripe,IconLabel)
%-----------------------------------------------------------------------

%-Only label Sun command tools
%-----------------------------------------------------------------------
Term        = getenv('TERM');
if ~strcmp(Term,'sun-cmd'), return, end

%-Work out label text
%-----------------------------------------------------------------------
User        = lui_spm('GetUser');
[null,Host] = unix('echo `hostname` | sed -e ''s/\..*$//''');
Host        = Host(1:length(Host)-1); 
v           = lui_spm('MLver');

if nargin<3, IconLabel = ['MatLab',v(1)]; end
if nargin<2, WinStripe = [User,' - ',Host,' : MatLab ',v]; end

%-Set window stripe
%-----------------------------------------------------------------------
disp([']l' WinStripe '\]L' IconLabel '\'])


%=======================================================================
case 'popupcb'               %-Callback handling utility for PopUp menus
%=======================================================================
% lui_spm('PopUpCB',h)
%-----------------------------------------------------------------------
if nargin<2, h=gcbo; else, h=varargin{2}; end
v   = get(h,'Value');
if v==1, return, end
set(h,'Value',1)
CBs = get(h,'UserData');
evalin('base',CBs{v-1})


%=======================================================================
case 'getuser'                                           %-Get user name
%=======================================================================
% str = lui_spm('GetUser',fmt)
%-----------------------------------------------------------------------
str = lui_spm_platform('user');
if ~isempty(str) & nargin>1, str = sprintf(varargin{2},str); end
varargout = {str};


%=======================================================================
case 'beep'                                %-Emit a keyboard "bell" beep
%=======================================================================
% lui_spm('Beep')
fprintf('%c',7)


%=======================================================================
case 'time'                          %-Return formatted date/time string
%=======================================================================
% [timestr, date_vec] = lui_spm('Time')
%-----------------------------------------------------------------------
tmp = clock;
varargout = {sprintf('%02d:%02d:%02d - %02d/%02d/%4d',...
			tmp(4),tmp(5),floor(tmp(6)),tmp(3),tmp(2),tmp(1)),...
		tmp};


%=======================================================================
case 'pointer'                 %-Set mouse pointer in all MatLab windows
%=======================================================================
% lui_spm('Pointer',Pointer)
%-----------------------------------------------------------------------
if nargin<2, Pointer='Arrow'; else, Pointer=varargin{2}; end
set(get(0,'Children'),'Pointer',Pointer)


%=======================================================================
case {'alert','alert"','alert*','alert!'}                %-Alert dialogs
%=======================================================================
% h = lui_spm('alert',Message,Title,CmdLine,wait)
%-----------------------------------------------------------------------

%- Globals 
%-----------------------------------------------------------------------
if nargin<5, wait    = 0;  else, wait    = varargin{5}; end
if nargin<4, CmdLine = []; else, CmdLine = varargin{4}; end
if nargin<3, Title   = ''; else, Title   = varargin{3}; end
if nargin<2, Message = ''; else, Message = varargin{2}; end
Message = cellstr(Message);

if isreal(CmdLine)
	CmdLine  = lui_spm('CmdLine',CmdLine);
	CmdLine2 = 0;
else
	CmdLine  = lui_spm('CmdLine');
	CmdLine2 = 1;
end
timestr = lui_spm('Time');
lui_spmv    = lui_spm('ver');

switch(lower(Action))
case 'alert',	icon = 'none';	str = '--- ';
case 'alert"',	icon = 'help';	str = '~ - ';
case 'alert*',	icon = 'error'; str = '* - ';
case 'alert!',	icon = 'warn';	str = '! - ';
end

if CmdLine | CmdLine2
	Message(strcmp(Message,'')) = {' '};
	tmp = sprintf('%s: %s',lui_spmv,Title);
	fprintf('\n    %s%s  %s\n\n',str,tmp,repmat('-',1,62-length(tmp)))
	fprintf('        %s\n',Message{:})
	fprintf('\n        %s  %s\n\n',repmat('-',1,62-length(timestr)),timestr)
	h = [];
end

if ~CmdLine
	tmp = max(size(char(Message),2),42) - length(lui_spmv) - length(timestr);
	str = sprintf('%s  %s  %s',lui_spmv,repmat(' ',1,tmp-4),timestr);
	h   = msgbox([{''};Message(:);{''};{''};{str}],...
		sprintf('%s%s: %s',lui_spmv,lui_spm('GetUser',' (%s)'),Title),...
		icon,'non-modal');
	drawnow
	set(h,'windowstyle','modal');
end

if wait
	if isempty(h)
		input('        press ENTER to continue...');
	else
		uiwait(h)
		h = [];
	end
end

if nargout, varargout = {h}; end


%=======================================================================
case {'fnbanner','sfnbanner','ssfnbanner'}  %-Text banners for functions
%=======================================================================
% lui_spmid = lui_spm('FnBanner', Fn,FnV)
% lui_spmid = lui_spm('SFnBanner',Fn,FnV)
% lui_spmid = lui_spm('SSFnBanner',Fn,FnV)
%-----------------------------------------------------------------------
time = lui_spm('time');
str  = lui_spm('ver');
if nargin>=2, str = [str,': ',varargin{2}]; end
if nargin>=3, str = [str,' (v',varargin{3},')']; end

switch lower(Action)
case 'fnbanner'
	tab = '';
	wid = 72;
	lch = '=';
case 'sfnbanner'
	tab = sprintf('\t');
	wid = 72-8;
	lch = '-';
case 'ssfnbanner'
	tab = sprintf('\t\t');
	wid = 72-2*8;
	lch = '-';
end

fprintf('\n%s%s',tab,str)
fprintf('%c',repmat(' ',1,wid-length([str,time])))
fprintf('%s\n%s',time,tab)
fprintf('%c',repmat(lch,1,wid)),fprintf('\n')
varargout = {str};


%=======================================================================
case 'fnuisetup'                %-Robust UI setup for main lui_spm functions
%=======================================================================
% [Finter,Fgraph,CmdLine] = lui_spm('FnUIsetup',Iname,bGX,CmdLine)
%-----------------------------------------------------------------------
if nargin<4, CmdLine=lui_spm('CmdLine'); else, CmdLine=varargin{4}; end
if nargin<3, bGX=1; else, bGX=varargin{3}; end
if nargin<2, Iname=''; else, Iname=varargin{2}; end
if CmdLine
	Finter = lui_spm_figure('FindWin','Interactive');
	if ~isempty(Finter), lui_spm_figure('Clear',Finter), end
	%if ~isempty(Iname), fprintf('%s:\n',Iname), end
else
	Finter = lui_spm_figure('GetWin','Interactive');
	lui_spm_figure('Clear',Finter)
	if ~isempty(Iname)
		str = sprintf('%s (%s): %s',lui_spm('ver'),lui_spm('GetUser'),Iname);
	else
		str = '';
	end
	set(Finter,'Name',str)
end

if bGX
	Fgraph = lui_spm_figure('GetWin','Graphics');
	lui_spm_figure('Clear',Fgraph)
else
	Fgraph = lui_spm_figure('FindWin','Graphics');
end
varargout = {Finter,Fgraph,CmdLine};	


%=======================================================================
case 'figname'                                %-Robust lui_spm figure naming
%=======================================================================
% F = lui_spm('FigName',Iname,F,CmdLine)
%-----------------------------------------------------------------------
if nargin<4, CmdLine=lui_spm('CmdLine'); else, CmdLine=varargin{4}; end
if nargin<3, F='Interactive'; else, F=varargin{3}; end
if nargin<2, Iname=''; else, Iname=varargin{2}; end

%if ~isempty(Iname), fprintf('\t%s\n',Iname), end
if CmdLine, varargout={[]}; return, end
F = lui_spm_figure('FindWin',F);
if ~isempty(F) & ~isempty(Iname)
	set(F,'Name',sprintf('%s (%s): %s',lui_spm('ver'),lui_spm('GetUser'),Iname))
end
varargout={F};


%=======================================================================
case 'gui_filedelete'                                %-GUI file deletion
%=======================================================================
% lui_spm('GUI_FileDelete')
%-----------------------------------------------------------------------
P = cellstr(lui_spm_select(Inf,'.*','Select file(s) to delete'));
n = numel(P);
if n==0
	lui_spm('alert"','Nothing selected to delete!','file delete',0);
	return
elseif n<4
	str=[{' '};P];
elseif n<11
	str=[{' '};P;{' ';sprintf('(%d files)',n)}];
else
	str=[{' '};P(1:min(n,10));{'...';' ';sprintf('(%d files)',n)}];
end
if lui_spm_input(str,-1,'bd','delete|cancel',[1,0],[],'confirm file delete')
	lui_spm_unlink(P{:})
	lui_spm('alert"',P,'file delete',1);
end


%=======================================================================
case 'show'                   %-Bring visible MatLab windows to the fore
%=======================================================================
% Fs = lui_spm('Show')
%-----------------------------------------------------------------------
cF = get(0,'CurrentFigure');
Fs = get(0,'Children');
Fs = findobj(Fs,'flat','Visible','on');
for F=Fs', figure(F), end
set(0,'CurrentFigure',cF)
lui_spm('FnBanner','GUI show');
varargout={Fs};


%=======================================================================
case 'clear'                                             %-Clear lui_spm GUI
%=======================================================================
% lui_spm('Clear',Finter, Fgraph)
%-----------------------------------------------------------------------
if nargin<3, Fgraph='Graphics'; else, Fgraph=varargin{3}; end
if nargin<2, Finter='Interactive'; else, Finter=varargin{2}; end
lui_spm_figure('Clear',Fgraph)
lui_spm_figure('Clear',Finter)
lui_spm('Pointer','Arrow')
lui_spm_select('clearvfiles');
lui_spm_conman('Initialise','reset');
local_clc, lui_spm('FnBanner','GUI cleared');
fprintf('\n');
%evalin('Base','clear')


%=======================================================================
case 'help'                                  %-Pass through for lui_spm_help
%=======================================================================
% lui_spm('Help',varargin)
%-----------------------------------------------------------------------
if nargin>1, lui_spm_help(varargin{2:end}), else, lui_spm_help, end


%=======================================================================
otherwise                                        %-Unknown action string
%=======================================================================
error('Unknown action string')

%=======================================================================
end


%=======================================================================
function realign_unwarp(ob,varargin)
% Choose either realign or unwarp
%=======================================================================
if get(ob,'Value')==1,
	lui_spm_jobman('interactive','','jobs.spatial.realign');
else
	lui_spm_jobman('interactive','','jobs.spatial.realignunwarp');
end


%=======================================================================
function local_clc
%=======================================================================
if ~(exist('isdeployed') && isdeployed),
    clc
end

%=======================================================================
function check_installation
%=======================================================================
d = lui_spm('Dir');

% check the search path
if ~ismember(lower(d),lower(strread(path,'%s','delimiter',pathsep)))
    error(sprintf([...
        'You do not appear to have the MATLAB search path\n'...
        'set up to include your lui_spm5 distribution.  This\n'...
        'means that you can start lui_spm in this directory,\n'...
        'but if you change to another directory then MATLAB\n'...
        'will be unable to find the lui_spm functions.  You\n'...
        'can use the editpath command in MATLAB to set it up.\n'...
        'For more information, try typing the following:\n'...
        '    help path\n    help editpath\n']));
end

% Ensure that the original release - as well as the updates - was installed.
if ~exist(fullfile(d,'lui_spm_showdoc.m')), % This is a file that should not have changed
    if isunix,
        error(sprintf([...
            'There appears to be some problem with the installation.\n'...
            'The original lui_spm5.tar.gz distribution should be installed\n'...
            'and the updates installed on top of this.  Unix commands\n'...
            'to do this are:\n'...
            '   gunzip < lui_spm5.tar.gz | tar xvf -\n'...
            '   cd lui_spm5\n'...
            '   gunzip < Updates_????.tar.gz | tar xvf -\n'...
            'You may need help from your local network administrator.']));
    else
        error(sprintf([...
            'There appears to be some problem with the installation.\n'...
            'The original lui_spm5.tar.gz distribution should be installed\n'...
            'and the updates installed on top of this.  If in doubt,\n'...
            'consult your local network administrator.']));
    end
end

% Ensure that the files were unpacked correctly
if ispc
    try
        t = load(fullfile(d,'Split.mat'));
    catch
        error(sprintf([...
            'There appears to be some problem reading the MATLAB\n'...
            '.mat files from the lui_spm distribution.  This is probably\n'...
            'something to do with the way that the distribution was\n'...
            'unpacked.  If you used WinZip, then ensure that\n'...
            'TAR file smart CR/LF conversion is disabled\n'...
            '(under the Miscellaneous Configuration Options).\n']));
    end
    if ~exist(fullfile(d,'toolbox','DARTEL','diffeo3d.c'),'file'),
        error(sprintf([...
            'There appears to be some problem with the installation.\n'...
            'This is probably something to do with the way that the\n'...
            'distribution was unbundled from the original .tar.gz files.'...
            'Please ensure that the files are unpacked so that the\n'...
            'directory structure is retained.\n']));
    end
end

% Check the mex files
try
    lui_spm_atranspa(1);
catch
    error(sprintf([...
        'lui_spm uses a number of "mex" files, which are compiled functions.\n'...
        'These need to be compiled for the various platforms on which lui_spm\n'...
        'is run.  At the FIL, where lui_spm is developed, the number of\n'...
        'computer platforms is limited.  It is therefore not possible to\n'...
        'release a version of lui_spm that will run on all computers. See\n'...
        '%s%csrc%cMakefile for information about how \n'...
        'to compile mex files for %s in MATLAB %s.\n'],...
        d,filesep,filesep,computer,version));
end

if lui_spm_matlab_version_chk('7') >=0,
    warning('off','MATLAB:dispatcher:ShadowedMEXExtension');
end

