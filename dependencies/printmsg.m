function varargout = printmsg(msg, varargin)
% PRINTMSG Create and print a formatted message with title
%
%	USAGE: varargout = printmsg(msg, varargin)
%

% --------------------------- Copyright (C) 2014 ---------------------------
%	Author: Bob Spunt
%	Email: bobspunt@gmail.com
% 
%	$Created: 2014_09_27
% _________________________________________________________________________
def = { ...
        'msgtop',          '_',        ...
        'msgbottom',        '_',          ...
        'msgwidth',     75, ...
        'msgtitle',     '', ...
        'nodisplay',      0,      ...
        };
vals = setargs(def, varargin);
if nargin < 1, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
if iscell(msg), msg = char(msg); end
if iscell(msgtitle), msgtitle = char(msgtitle); end
msgwidth        = max([msgwidth length(msgtitle)+10]); 
msgtopborder    = repmat(msgtop,1,msgwidth);
msgbottomborder = repmat(msgbottom,1,msgwidth);
if ~isempty(msgtitle), msgtitle = sprintf('%s %s %s', msgtop, strtrim(msgtitle), msgtop); end
titleln         = length(msgtitle);
msgln           = length(msg); 
msgtopborder(floor(.5*msgwidth-.5*titleln):floor(.5*msgwidth-.5*titleln) + titleln-1) = msgtitle;
fmtmessage      = repmat(' ', 1, msgwidth);
fmtmessage(floor(.5*msgwidth-.5*msgln):floor(.5*msgwidth-.5*msgln) + msgln-1) = msg;
fmtmessage      = sprintf('%s\n\n%s\n%s', msgtopborder, fmtmessage, msgbottomborder);
if ~nodisplay, disp(fmtmessage); end
if nargout > 0, varargout{1} = fmtmessage; end
end