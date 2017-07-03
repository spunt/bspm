%_______________________________________________________________________
%
% Compute variability for image data
%_______________________________________________________________________
%
% Input
%
%    img | image data  (2D matrix, scans x voxel)
% metric | variability metric (string)
%
%_______________________________________________________________________
%
% Output
%
% res_img | variability of image data (2D matrix, 1 x voxel)
%
%_______________________________________________________________________
%
% This file is part of the Variability Toolbox for SPM
% Published by the Lifespan Neural Dynamics Group
% Provided under the GNU General Public License
% See LICENSE for terms and conditions

function res_img = shared_var(img, metric)

  if size(img, 1) < 3
    error('Need at least three trials per condition for variability computation.')
  end

  switch metric
    case 'var'
      res_img = squeeze(var(img));
    case 'sd'
      res_img = squeeze(std(img));
    case 'mssd'
      res_img = squeeze(sum(diff(img).^2)/(size(img,1)-1));
    case 'sqrt_mssd'
      res_img = sqrt(squeeze(sum(diff(img).^2)/(size(img,1)-1)));
    otherwise
      error('Invalid variability metric.');
  end

end
