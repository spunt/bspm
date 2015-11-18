%_______________________________________________________________________
%
% Save image data as nifti file
%_______________________________________________________________________
%
% Input
%
%  input_file | nifti file used as template for header (string)
% output_file | nifti file where image data is written (string)
%         img | image data (2D double matrix)
%
%_______________________________________________________________________
%
% This file is part of the Variability Toolbox for SPM
% Published by the Lifespan Neural Dynamics Group
% Provided under the GNU General Public License
% See LICENSE for terms and conditions

function shared_save_nifti(input_file, output_file, img)

  cfg = shared_config;

  %%
  %% prepare file header and image data
  %%
	input_hdr = spm_vol([input_file ',1']);
  hdr.dim = input_hdr.dim;
  hdr.mat = input_hdr.mat;
  hdr.fname = output_file;
	hdr.descrip = [cfg.name ' ' cfg.version];
  hdr.dt = [spm_type('float32') spm_platform('bigend')];

  dim = numel(size(img));

  if dim == 4

    for n = 1:size(img,1)
      hdr.n = [n 1];
      spm_write_vol(hdr, squeeze(img(n,:,:,:)));
    end

    %% 
    %% remove automatically generated .mat file
    %%
    [pth, name, ext] = fileparts(output_file);
    result_mat = fullfile(pth,[name '.mat']);
    if exist(result_mat)
      delete(result_mat);
    end
    
  else

    spm_write_vol(hdr, img);
  end

end
