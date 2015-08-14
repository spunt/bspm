function nii_smooth(P, FWHM)
% Creates smoothed image with prefix 's'
%   Returns output filename
%  P : filename(s)
%  FWHM : full-width half-maximum of smoothing kernel
% Example
%   nii_smooth ('C:\ct\script\xwsctemplate_final.nii',3);

% if nargin <1 %no files
%  P = spm_select(inf,'image','Select images to smooth');
% end;
% if nargin < 2 %no FWHM specified
%  FWHM = 8;
% end;
if nargin < 2, error('USAGE: nii_smooth(P, FWHM)'); end
if iscell(P), P = char(P); end
for i=1:size(P,1)
  ref = deblank(P(i,:));
  [pth,nam,ext] = spm_fileparts(ref);
  src = fullfile(pth,[nam ext]);
  smth = fullfile(pth, ['s' nam ext]);
  spm_smooth(src,smth,FWHM,0);  
end