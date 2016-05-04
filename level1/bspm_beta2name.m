function name = bspm_beta2name(beta)
% BSPM_BETA2NAME
% USAGE: name = bspm_beta2name(beta)
%
if nargin<1, mfile_showhelp; return; end
if iscell(beta), beta = char(beta); end
h       = spm_vol(beta);
name    = {h.descrip}'; 
name    = regexprep(name, '^.*Sn\(\d\)\s', ''); 
name    = regexprep(name, '\*bf\(\d\)', '');
end
function mfile_showhelp(varargin)
    % MFILE_SHOWHELP
    ST = dbstack('-completenames');
    if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
    eval(sprintf('help %s', ST(2).file));  
end
