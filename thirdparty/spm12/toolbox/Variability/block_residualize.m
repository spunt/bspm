%_______________________________________________________________________
%
% Regress out nuisance variables such as motion parameters
%_______________________________________________________________________
%
% Input
%
%     img | image data of a session (2D double matrix, scans x voxel)
%     reg | regressor data of a session (double matrix, scans x regressors)
%    sess | number of session (integer)
%
%_______________________________________________________________________
%
% Output
%
% res_img | residualized image data (2D double matrix, scans x voxel)
%
%_______________________________________________________________________
%
% This file is part of the Variability Toolbox for SPM
% Published by the Lifespan Neural Dynamics Group
% Provided under the GNU General Public License
% See LICENSE for terms and conditions

function res_img = block_residualize(img, reg, sess)

  cfg = shared_config;

  %%
  %%  compute the temporal mean for all voxels
  %%
	avg = mean(img);
	step_size = 5000;
	num_vox = size(img,2);
  num_step = ceil(num_vox/step_size);

  if cfg.gui
	  label_text = sprintf('Computing Session %i', sess);
	  label_format = '\fontsize{16}';
	  label = [label_format label_text];
	  spm_progress_bar('Init', 100, label, '', 't');
  end

  %%
  %% residualize step-wise to allow progress indicator
  %%
	for step = 1:num_step

    vox_begin = step_size * (step - 1) + 1;

    if step < num_step
      vox_end = vox_begin + step_size - 1;
    else
      vox_end = num_vox;
    end

    vox = vox_begin:vox_end;
    res_img(:,vox) = residualize(img, avg, reg, vox);

		if cfg.gui
      spm_progress_bar('Set', round(100*step/num_step));
    end
	end

end

%%
%% wrapper for residualization
%%
function result_img = residualize(img, avg, reg, vox)
  img = img(:,vox);
  avg = avg(vox);
  num_scan = size(img,1);
  resid_img = tbx_regress(img, reg);
  mean_img = repmat(avg, num_scan, 1);
  result_img = resid_img + mean_img;
end

%%
%% residualize using linear regression
%%
function result_img = tbx_regress(img, reg)

	%% perform linear regression on regressors and image data
	num_scan = size(img, 1);
	intercept = ones(num_scan, 1);
	reg_coef = [intercept reg] \ img;
	
	%% predict scans with linear model and remove correlations
	predict = [intercept reg] * reg_coef;
	result_img = img - predict;
end
