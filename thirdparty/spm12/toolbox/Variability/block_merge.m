%_______________________________________________________________________
%
% Merge image data of all sessions condition-wise
%_______________________________________________________________________
%
% Input
%
%  session | session data of first-level design (struct)
%    index | condition-wise scan indices (3D cell array)
%    coord | mask for image data (double vector)
%
%_______________________________________________________________________
%
% Output
%
% img | condition-wise arranged image data (2D cell array)
%
%_______________________________________________________________________
%
% This file is part of the Variability Toolbox for SPM
% Published by the Lifespan Neural Dynamics Group
% Provided under the GNU General Public License
% See LICENSE for terms and conditions

function img = block_merge(session, index, coord)

  %%
  %% retrieve total number of blocks per condition
  %%
  cnd_blk_total = zeros(numel(index),1);
  for sess = 1:numel(session)
    for cnd = 1:numel(index)
      cnd_blk_total(cnd) = cnd_blk_total(cnd) + numel(index{cnd}{sess});
    end
  end

  for cnd = 1:numel(index)
    img{cnd} = cell(cnd_blk_total(cnd),1);
    blk_count{cnd} = 0;
  end

  %%
  %% session-wise loading avoids high memory peaks
  %%
  for sess = 1:numel(session)

    scan_img = block_load_data(session, sess, coord);

    %%
  	%% residualize if regressors have been specified
    %%
    if isfield(session(sess),'regress') ...
    && isstruct(session(sess).regress) ...
    && numel(session(sess).regress) >= 1

      for r = 1:numel(session(sess).regress)
        if ~isnumeric(session(sess).regress(r).val)
          error('Regressor field contains non-numerical data')
        end
        size_reg = size(session(sess).regress(r).val,1);
        size_img = size(scan_img,1);
        if size_reg ~= size_img
          msg = 'Differing number of scans in image data and regressor field ';
          msg = [msg sprintf('(%i vs. %i)',size_img, size_reg)];
          error(msg)          
        end
      end

      regressor = [session(sess).regress.val];
      scan_img = block_residualize(scan_img, regressor, sess);
    end

    %%
    %% merge blocks condition-wise
    %%
  	for cnd = 1:numel(index)
  		for blk = 1:numel(index{cnd}{sess})

        %%
        %% image data of current condition, session and block
        %%
  			block_img = scan_img(index{cnd}{sess}{blk},:);

        %%
  			%% normalize to block global mean
        %%
  			block_img = 100 * block_img / mean(mean(block_img));

        %%
        %% using cell array to delay detrending
        %%
        cur_blk = blk_count{cnd} + blk;
        num_scan = size(block_img, 2);
        img{cnd}{cur_blk} = block_img;
        img{cnd};

  		end % block
      blk_count{cnd} = blk_count{cnd} + numel(index{cnd}{sess});
  	end % condition
  end % session
end
