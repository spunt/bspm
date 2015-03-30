function fn = bspm_expand4D(fn4D)
% BSPM_EXPAND4D
%   USAGE: fn = bspm_expand4D(fn4D)
%
if nargin<1, error('USAGE: fn = bspm_expand4D(fn4D)'); end
if iscell(fn4D), fn4D = char(fn4D); end
[pcim, ncim, ecim] = fileparts(fn4D);
if strcmp(ecim, '.gz')
    gunzip(fn4D); 
    fn4D = fullfile(pcim, ncim); 
end
nvol    = length(spm_vol(fn4D)); 
append  = cellfun(@num2str, num2cell(1:nvol)', 'Unif', false);
fn      = strcat(fn4D, {','}, append);
end
 
 
 
 
