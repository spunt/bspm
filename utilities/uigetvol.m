function vol        = uigetvol(message, multitag, defaultdir)
% UIGETVOL Dialogue for selecting image volume file
%
%   USAGE: vol = uigetvol(message, multitag)
%       
%       message = to display to user
%       multitag = (default = 0) tag to allow selecting multiple images
%
% EX: img = uigetvol('Select Image to Process'); 
%
if nargin < 3, defaultdir = pwd; end
if nargin < 2, multitag = 0; end
if nargin < 1, message = 'Select Image File'; end
if ~multitag
    [imname, pname] = uigetfile({fullfile(defaultdir, '*.img;*.nii;*.nii.gz'), 'Image File (*.img, *.nii, *.nii.gz)'; '*.*', 'All Files (*.*)'}, message);
else
    [imname, pname] = uigetfile({fullfile(defaultdir, '*.img;*.nii;*.nii.gz'), 'Image File (*.img, *.nii, *.nii.gz)'; '*.*', 'All Files (*.*)'}, message, 'MultiSelect', 'on');
end
if isequal(imname,0) || isequal(pname,0)
    vol = [];
else
    vol = fullfile(pname, strcat(imname));
end