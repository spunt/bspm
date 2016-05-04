function bspm_display_message(message, border, titlestr)
% BOB_DISPLAY_MESSAGE
%
% Prints "message" to screen with border character defined by "border"
% (default = '='), and with optional title defined by "titlestr"
%
% USAGE: bspm_display_message(message, border, titlestr)
%
if nargin<1, mfile_showhelp; return; end
if nargin<2, border = '='; end
if nargin<3, titlestr = []; end
mln = length(message);
boxTop(1:mln)=border;
if isempty(titlestr)
    fprintf('%s\n%s\n%s\n', boxTop, message, boxTop)
    return
end
bsp      = 5; % extra space to add to title
titlestr = strtrim(titlestr);
titlestr = [' ' titlestr ' '];
tln = length(titlestr);
if (tln + bsp*2) > mln
    spln = (tln + bsp*2 - mln)/2; 
    lrsp = [floor(spln) ceil(spln)]; 
    message = sprintf('%s%s%s', repmat(' ',1,lrsp(1)), message, repmat(' ',1,lrsp(2))); 
    mln = length(message);
    boxTop(1:mln)=border;
end
flankln = (mln-tln)/2;
if ~mod(flankln, 2), leftrightln = [flankln flankln]; 
else leftrightln = [floor(flankln) ceil(flankln)]; end
fprintf('\n%s%s%s\n%s\n%s\n', repmat(border,1,leftrightln(1)), titlestr, ...
  repmat(border,1,leftrightln(2)), message, boxTop);
end
 
 
 
 
