function bspm_display(imname)
% Image and header display
% USAGE: bspm_display(imname)
%__________________________________________________________________________
%
% spm_image is an interactive facility that allows orthogonal sections
% from an image volume to be displayed.  Clicking the cursor on either
% of the three images moves the point around which the orthogonal
% sections are viewed.  The co-ordinates of the cursor are shown both
% in voxel co-ordinates and millimeters within some fixed framework.
% The intensity at that point in the image (sampled using the current
% interpolation scheme) is also given. The position of the crosshairs
% can also be moved by specifying the co-ordinates in millimeters to
% which they should be moved.  Clicking on the horizontal bar above
% these boxes will move the cursor back to the origin  (analogous to
% setting the crosshair position (in mm) to [0 0 0]).
%
% The images can be re-oriented by entering appropriate translations,
% rotations and zooms into the panel on the left.  The transformations
% can then be saved by hitting the "Reorient images..." button.  The
% transformations that were applied to the image are saved to the header
% information of the selected images.  The transformations are considered
% to be relative to any existing transformations that may be stored.
% Note that the order that the transformations are applied in is the
% same as in spm_matrix.m.
%
% The ``Reset...'' button next to it is for setting the orientation of
% images back to transverse.  It retains the current voxel sizes,
% but sets the origin of the images to be the centre of the volumes
% and all rotations back to zero.
%
% The right panel shows miscellaneous information about the image.
% This includes:
%   Dimensions - the x, y and z dimensions of the image.
%   Datatype   - the computer representation of each voxel.
%   Intensity  - scalefactors and possibly a DC offset.
%   Miscellaneous other information about the image.
%   Vox size   - the distance (in mm) between the centres of
%                neighbouring voxels.
%   Origin     - the voxel at the origin of the co-ordinate system
%   Dir Cos    - Direction cosines.  This is a widely used
%                representation of the orientation of an image.
%
% There are also a few options for different resampling modes, zooms
% etc.  You can also flip between voxel space or world space.  If you
% are re-orienting the images, make sure that world space is specified.
% Blobs (from activation studies) can be superimposed on the images and 
% the intensity windowing can also be changed.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_image.m 4205 2011-02-21 15:39:08Z guillaume $
if nargin<1
    P = uigetvol('Select an image'); 
    if isempty(P), return; end
else
    P = imname;
end
if iscell(P), P = char(P); end
spm_image('Display', P); 