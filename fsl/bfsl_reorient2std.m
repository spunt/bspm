function bfsl_reorient2std(images)
% BFSL_REORIENT2STD
%
% USAGE: bfsl_reorient2std(images)
%
% -----------------------------------------------------
if nargin<1, mfile_showhelp; return; end
if ischar(images), images = cellstr(images); end
nim = length(images);
for i = 1:nim
    input = images{i};
    [p, n] = fileparts(input);
    output = [p filesep 'o' n '.nii.gz'];
    command = sprintf('system(''fslreorient2std %s %s'')',input, output);
    eval(command);
    gunzip(output);
    delete(output);
end
end

