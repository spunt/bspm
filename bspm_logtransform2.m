function bspm_logtransform2(input)
% BSPM_LOGTRANSFORM
% USAGE: bspm_logtransform(input)
%
%
%       
% Based on code provided by David Langers

% ---------- Copyright (C) 2014 ----------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, error('No input!'); end
fn = {'epipat'};
[status, msg] = checkfields(input, fn);
if ~status, error(msg); end
prefix = 'l';
epipat = input.epipat;
if iscell(epipat), epipat = char(epipat); end
in = files(epipat);

nv = length(in);
fprintf('\nReading data for %d volumes', nv);
[d h] = bspm_read_vol(in, 'reshape'); % read data
level = nanmean(d(d(:,1) > nanmean(d(:,1))/10)); % determine mean tissue signal
level = level/exp(0.0001*32768/2); % 16-bit integer range at 0.01% signal change resolution
if ~isempty(prefix)
    fn = {h.fname}';
    [p, n, e] = cellfun(@fileparts, fn, 'Unif', false);
    fn = strcat(p,filesep,{'l'},n, e);
    [h.fname] = fn{:};
end
fprintf('\nApplying log transform to ');
str1 = sprintf('%%0%dd', length(num2str(nv)));
str2 = repmat('\b', 1, 1+(2*length(num2str(nv)))); 
for f = 1:nv
  hdr = h(f);
  hdr.pinfo = [0.01; 0; 0];
  img = reshape(d(:,f),h(f).dim); 
  img = 100.0*log(max(img/level, 1.0));
  spm_write_vol(hdr, img);
  fprintf([str1 '/%d' str2], f, nv); 
end
fprintf('\nDONE\n');

% 
% % read representative (i.e. first) image to assess scale
% img = spm_read_vols(spm_vol(in{1}));
% 
% % determine mean tissue signal
% level = nanmean(img(img > nanmean(img(:))/10));
% 
% % use 16-bit integer range at 0.01% signal change resolution
% level = level/exp(0.0001*32768/2);
% 
% % apply transform to all images
% for f = 1:length(in)
%     
%   [path name ext] = fileparts(in{f});
%   if strcmp(ext,'.img')
%       hdr = spm_vol([path filesep name '.hdr']);
%   else
%       hdr = spm_vol(in{f});
%   end
%   img = spm_read_vols(hdr);
%   hdr.pinfo = [0.01; 0; 0];
%   hdr.fname = [path filesep prefix name '.nii'];
%   img = 100.0*log(max(img/level, 1.0));
%   spm_write_vol(hdr, img);
%   
% end






 
 
 
 
