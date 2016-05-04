function bfsl_reslice(images)
% BFSL_RESLICE Skull strip with FSL BET
%
%   USAGE: bfsl_reslice
% 
%   images = volumes to filter
%   f = fractional intensity threshold (0:1); default=.5; smaller values
%   give larger brain outline estimates
%   outprefix = prefix for output files
%
% ------------------------------------------------
if nargin==0, mfile_showhelp; return; end
if ischar(images), images = cellstr(images); end
nim = length(images);
for i = 1:nim
    input = images{i};
    [p n e] = fileparts(input);
    output = sprintf('%s/ds%s.nii.gz', p, n);
    command = sprintf('system(''fslmaths %s -subsamp2 %s'')',input, output);
    eval(command);
    gunzip(output);
    delete(output);
end