function bfsl_fslview(image)
% BFSL_FSLVIEW Call fslview from MATLAB
%
% USAGE: bfsl_fslview(image)
%
% fslview [-m 3d|ortho|lightbox] <baseimage> [-l lutname] [-b low,hi]
% 	[ <overlay> [-l lutname] [-b low,hi] ] ...
% fslview -m ortho,lightbox filtered_func_data thresh_zstat1 -t 0.5 thresh_zstat2 -l "Cool" -t 0.5
% 
% Optional arguments (You may optionally specify one or more of):
% 	-V,--verbose	switch on diagnostic messages
% 	-h,--help	display this message
% 	-m,--mode	Initial viewer mode. Comma separated list of: 3d; single, ortho; lightbox
% 
% 
% Per-image options
% 
% Usage:
% image [-l GreyScale] [-t 0.1] [-b 2.3,6]
% 	-l,--lut	Lookup table name. As per GUI, one of: Greyscale;
% 			"Red-Yellow"; "Blue-Lightblue"; Red; Green;
% 			Blue; Yellow; Pink; Hot; Cool; Copper, etc.
% 	-b,--bricon	Initial bricon range, e.g., 2.3,6
% 	-t,--trans	Initial transparency, e.g., 0.2
% -----------------------------------------------------
if nargin==0, system('fslview &'), end;
if iscell(image), image = char(image); end
system(['fslview ' image ' &']);