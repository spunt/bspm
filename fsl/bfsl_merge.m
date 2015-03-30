function bfsl_merge(images, outname)
% BFSL_MERGE Call fslview from MATLAB
%
% USAGE: bfsl_merge(images, outname)
%
% Usage: fslmerge <-x/y/z/t/a/tr> <output> <file1 file2 .......> [tr value in seconds]
%      -t : concatenate images in time
%      -x : concatenate images in the x direction
%      -y : concatenate images in the y direction
%      -z : concatenate images in the z direction
%      -a : auto-choose: single slices -> volume, volumes -> 4D (time series)
%      -tr : concatenate images in time and set the output image tr to the final option value
% --------------------------------------------------------------------------------------------
if nargin < 1, error('USAGE: bfsl_merge(images, outname)'); end
if ischar(images), images = cellstr(images); end
[p1,n1] = fileparts(images{1});
if nargin < 2
    charidx = intersect(regexp(n1,'\D'), regexp(n1, '\w')); 
    outname = fullfile(p1, sprintf('%s-4D-T%d.nii', n1(charidx), length(images))); 
elseif isempty(fileparts(outname))
    outname = fullfile(p1, outname);
end
cmd = sprintf(['fslmerge -t %s' repmat(' %s', 1, length(images))], outname, images{:}); 
system(cmd);