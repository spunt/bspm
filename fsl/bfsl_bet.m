function bfsl_bet(images, f, outprefix, nowaitbar, g, dorobust)
% BFSL_BET Skull strip with FSL BET
%
%   USAGE: bfsl_bet(images, f, outprefix, nowaitbar, g, dorobust)
% 
%   images = volumes to filter
%   f = fractional intensity threshold (0:1); default=.5; smaller values
%   give larger brain outline estimates
%   outprefix = prefix for output files
%   nowaitbar = option to suppress waitbar
%   -g <g>      vertical gradient in fractional intensity threshold (-1->1); default=0; positive values give larger brain outline at bottom, smaller at top
%   dorobust    = -R          robust brain centre estimation (iterates BET several times)
%
% ------------------------------------------------
if nargin<6, dorobust = 0; end
if nargin<5, g = 0; end
if nargin<4, nowaitbar = 0; end
if nargin<3, outprefix = 'b'; end
if nargin<2, f = 0.5; end
if nargin<1, mfile_showhelp; return; end
if ischar(images), images = cellstr(images); end
nim = length(images);
if ~nowaitbar
    h = waitbar(0,sprintf('Completed 0 of %d Volumes',nim));
end
for i = 1:nim
    input = images{i};
    [p n e] = fileparts(input);
    output = [p filesep outprefix n '.nii.gz'];
    if dorobust
        command = sprintf('system(''bet %s %s -f %2.2f -g %2.2f -R'')',input, output, f, g);
    else
        command = sprintf('system(''bet %s %s -f %2.2f -g %2.2f'')',input, output, f, g);
    end
    eval(command);
    gunzip(output);
    delete(output);
    if ~nowaitbar
        waitbar(i/nim,h,sprintf('Completed %d of %d Volumes',i,nim));
    end
end
if ~nowaitbar
close(h);
end

% Usage:    bet <input> <output> [options]
% 
% Main bet2 options:
%   -o          generate brain surface outline overlaid onto original image
%   -m          generate binary brain mask
%   -s          generate approximate skull image
%   -n          don't generate segmented brain image output
%   -f <f>      fractional intensity threshold (0->1); default=0.5; smaller values give larger brain outline estimates
%   -g <g>      vertical gradient in fractional intensity threshold (-1->1); default=0; positive values give larger brain outline at bottom, smaller at top
%   -r <r>      head radius (mm not voxels); initial surface sphere is set to half of this
%   -c <x y z>  centre-of-gravity (voxels not mm) of initial mesh surface.
%   -t          apply thresholding to segmented brain image and mask
%   -e          generates brain surface as mesh in .vtk format
% 
% Variations on default bet2 functionality (mutually exclusive options):
%   (default)   just run bet2
%   -R          robust brain centre estimation (iterates BET several times)
%   -S          eye & optic nerve cleanup (can be useful in SIENA)
%   -B          bias field & neck cleanup (can be useful in SIENA)
%   -Z          improve BET if FOV is very small in Z (by temporarily padding end slices)
%   -F          apply to 4D FMRI data (uses -f 0.3 and dilates brain mask slightly)
%   -A          run bet2 and then betsurf to get additional skull and scalp surfaces (includes registrations)
%   -A2 <T2>    as with -A, when also feeding in non-brain-extracted T2 (includes registrations)
% 
% Miscellaneous options:
%   -v          verbose (switch on diagnostic messages)
%   -h          display this help, then exits
%   -d          debug (don't delete temporary intermediate images)