function out = bspm_expand4Dfilename(in)    
% USAGE: out = bspm_expand4Dfilename(in)  
%
if nargin<1, disp('USAGE: out = bspm_expand4Dfilename(in)'); return; end
if iscell(in), in = char(in); end
[pcim, ncim, ecim] = fileparts(char(in));
if strcmp(ecim, '.gz')
    gunzip(in); 
    in = fullfile(pcim, ncim); 
end
nvol    = length(spm_vol(in)); 
append  = cellfun(@num2str, num2cell(1:nvol)', 'Unif', false);
out     = strcat(in, {','}, append);
end