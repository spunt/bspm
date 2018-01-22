function bfsl_reorient2std(images)
% BFSL_REORIENT2STD
%
% Usage: fslreorient2std <input_image> [output_image]
%
% fslreorient2std is a tool for reorienting the image to match the
% approximate orientation of the standard template images (MNI152).
% It only applies 0, 90, 180 or 270 degree rotations.
% It is not a registration tool.
% It requires NIfTI images with valid orientation information
% in them (seen by valid labels in FSLView).  This tool
% assumes the labels are correct - if not, fix that before using this.
% If the output name is not specified the equivalent transformation
%  matrix is written to the standard output
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

