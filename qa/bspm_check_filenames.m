function fn = bspm_check_filenames(inname)
% BSPM_CHECK_FILENAMES
% 
%   USAGE outname = bspm_check_orientations(inname)
%       imname can be char or cell 
%       output is a cell array
%
if nargin==0, mfile_showhelp; return; end
if ischar(inname), inname = cellstr(inname); end
inname              = regexprep(inname, ',\d+$', '');
[pcim, ncim, ecim]  = cellfun(@fileparts, inname, 'unif', false);
fn = []; 
for i = 1:length(inname)
    
    if ~ismember(ecim{i}, {'.gz' '.nii' '.img'})
        continue
    end
    
    % | If compressed, try to decompress
    if strcmp(ecim{i}, '.gz')
        try
            pigz(inname{i})
        catch
            gunzip(inname{i});
        end
        inname{i} = cellstr(fullfile(pcim{i}, ncim{i}));  
    end
    
    % | Check to see if 4D
    nii     = nifti(inname{i}); 
    nvol    = size(nii.dat, 4); 
    append  = cellfun(@num2str, num2cell(1:nvol)', 'Unif', false);
    fn      = [fn; strcat(inname{i}, {','}, append)];
end
end
%  
% 
% 
% 
% if ~isstruct(inname)
%     if iscell(inname), inname = char(inname); end
%     try inname = spm_vol(inname); catch err, rethrow err; end
% end
% flag = 0;
% volinfo.fname   = {inname.fname}'; 
% volinfo.dims    = cat(1,inname.dim);
% volinfo.mats    = reshape(cat(3,inname.mat),[16,numel(inname)]); 
% if any(any(diff(volinfo.dims,1,1),1))
%     flag = flag + 1;
%     volinfo.uniquedims = unique(volinfo.dims, 'rows'); 
%     if verbose, disp('The images do not all have the same dimensions.'); end 
% end
% if any(any(abs(diff(volinfo.mats,1,2))>1e-4))
%     flag = flag + 2; 
%     volinfo.uniquemats = unique(volinfo.mats', 'rows');
%     if verbose, disp('The images do not all have same orientation and/or voxel sizes.'); end
% end
% end
