function bfsl_optibet(images)
% BFSL_OPTIBET 
%
%   USAGE: bfsl_optibet(images)
% 
%     How To Use:   sh optiBET.sh -i <input_image> -options
% 
%     * if option is -f script uses FSL for initial extraction (default)
%     * if option is -a script uses AFNI for initial extraction
%     * if option is -o script uses MNI152_T1_1mm_brain_mask.nii.gz for mask (default)
%     * if option is -t script uses MNI152_T1_2mm_brain_mask.nii.gz for mask
%     * if option is -g script uses avg152T1_brain.nii.gz for mask
%     * if option is -d use debug mode (will NOT delete intermediate files)
%     * script requires proper installation of FSL and AFNI
%     * input image should be in standard orientation
%     * use .nii.gz image for input
%     * outputs binarized brain-extraction mask, saved as:  <input_image>_optiBET_brain_mask.nii.gz
%     * and full-intensity brain-extraction, saved as: <input_image>_optiBET_brain.nii.gz
%
%
% in .bash_profile:
% SHELLDIR=/Users/bobspunt/Github/osx-sync/shell
% PATH=${PATH}:${SHELLDIR}
% export SHELLDIR PATH
%
% ------------------------------------------------
if nargin<1, mfile_showhelp; return; end
if ischar(images), images = cellstr(images); end
nim = length(images); 
for i = 1:nim
    input       = images{i};
    [p, n, e]   = fileparts(input);
    if ~strcmp(e, '.gz')
       if strcmp(e, '.img'), bspm_img2nii(input); input = fullfile(p, [n '.nii']); end
       gzip(input);
       old      = input; 
       input    = strcat(input, '.gz'); 
       delete(old);
    end
    fprintf('\n | - Working on Image %d of %d', i, nim);
    command     =   sprintf('system(''sh optiBET.sh -i %s'')', input);
    eval(command);
end
fprintf('\n | - Done\n');  
 
 
 
 
