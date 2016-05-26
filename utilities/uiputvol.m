function vol    = uiputvol(defname, prompt)
% UIPUTVOL Dialogue for selecting output directory/name for volume
%
%   USAGE: vol    = uiputvol(defname, prompt)
%       
%       defname     = default outputname (can include path)
%       prompt      = message to display to user
%
if nargin < 1, defname = 'myimage.nii'; end
if nargin < 2, prompt = 'Save image as'; end
[imname, pname] = uiputfile({'*.img; *.nii', 'Image File'; '*.*', 'All Files (*.*)'}, prompt, defname);
if isequal(imname,0) || isequal(pname,0)
    vol = [];
else
    vol = fullfile(pname, imname); 
end