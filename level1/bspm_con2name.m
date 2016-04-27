function name = bspm_con2name(con, numberit)
% BSPM_CON2NAME
% USAGE: name = bspm_con2name(con)
%
if nargin<2, numberit = 0; end
if nargin==0, mfile_showhelp; return; end
if iscell(con), con = char(con); end
h       = spm_vol(con);
name    = {h.descrip}';
rmpat   = {
    '.*\d+:\s' 
    '(\s\W\sAll\sS\w*)$'
    };
for i = 1:length(rmpat), name = regexprep(name, rmpat{i}, ''); end
name    = regexprep(name, '\s', '_'); 
if numberit
    ncon = length(name);
    fmt = ['%0' num2str(length(num2str(ncon))) 'd']; 
    for i = 1:length(name)
        name{i} = sprintf(['C' fmt '_%s'], i, name{i});
    end
end

