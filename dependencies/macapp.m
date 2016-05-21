function macapp(appname, openquit, varargin)
% MACAPP Open/Close OSX Application
%
%  USAGE: macapp(appname, openquit, varargin)
%   

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-08-04
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
def = { ... 
	'forcequit',		0	...
	};
vals = setargs(def, varargin);
if nargin<2, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
if iscell(appname), appname =char(appname); end
if forcequit
    [err, str] = system(sprintf('killall %s', appname));
    if err, fprintf('Error: %s\n', str); end; 
    return;
end
if isnumeric(openquit)
    if openquit, openquit = 'open';
    else openquit = 'quit'; end;e
else
    openquit = strtrim(lower(openquit));
    if strcmpi(openquit, 'close'), openquit = 'quit'; end
end
if ~ismember(openquit, {'open' 'quit'})
    fprintf('\n - | INPUT %s IS INVALID. VALID OPTIONS ARE: open, quit\n\n', openquit); 
end
if strcmpi(openquit, 'open')
    cmd = sprintf('open -a "%s"', appname); 
else
    cmd = sprintf('osascript -e ''quit app "%s"''', appname);
end
[err, str] = system(cmd); 
if err, fprintf('Error: %s\n', str); end; 