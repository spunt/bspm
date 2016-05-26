function fn = bspm_expand4D(fn4D)
% BSPM_EXPAND4D
%   USAGE: fn = bspm_expand4D(fn4D)
%
if nargin<1, mfile_showhelp; return; end
if ischar(fn4D), fn4D = cellstr(fn4D); end
fn      = [];
fn4D    = regexprep(fn4D, ',\d+$', '');
[pcim, ncim, ecim] = cellfun(@fileparts, fn4D, 'unif', false); 
for i = 1:length(fn4D)
    
    % | If compressed, try to decompress
    if strcmp(ecim{i}, '.gz')
        try
            pigz(fn4D{i})
        catch
            gunzip(fn4D{i});
        end
        fn4D = cellstr(fullfile(pcim{i}, ncim{i}));  
    end
    nii     = nifti(fn4D{i}); 
    nvol    = size(nii.dat, 4); 
    append  = cellfun(@num2str, num2cell(1:nvol)', 'Unif', false);
    fn      = [fn; strcat(fn4D{i}, {','}, append)]; 
end
end
 
 
 
 
