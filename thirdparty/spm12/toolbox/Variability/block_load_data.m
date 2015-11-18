%_______________________________________________________________________
%
% Read nifti data of a session from disk
%
% Reading the nifti files from disk is rather expensive,
% so we load all image data only once and keep it in memory.
%
% As this can be quite heavy in memory usage the volumes
% are flattened to enable applying a mask if supplied.
%
% e.g. 91x109x91 -> 1x902629 -> 1x151420
%
% 902629 vs. 151420 voxels per volume
% 2.6 GB vs. 440 MB per subject (3 runs)
%
%_______________________________________________________________________
%
% Input
%
% session | session data of first-level design (struct)
%    sess | number of session to load (numeric)
%   coord | mask for image data (numeric vector)
%
%_______________________________________________________________________
%
% Output
%
% img | image data of session (2D numeric matrix)
% 
%_______________________________________________________________________
%
% This file is part of the Variability Toolbox for SPM
% Published by the Lifespan Neural Dynamics Group
% Provided under the GNU General Public License
% See LICENSE for terms and conditions

function img = block_load_data(session, sess, coord)

  cfg = shared_config;

	scan_file = session(sess).scans;
	num_scan = numel(scan_file);
	hdr = cell2mat(spm_vol(scan_file));

  if numel(session) == 1
    label_text = sprintf('Reading Session');
  else
    label_text = sprintf('Reading Session %i', sess);
  end
  
  if cfg.gui
	  label = [cfg.format label_text];
	  spm_progress_bar('Init', 100, label, '', 't');
  else
    fprintf('%+*s', cfg.pad, label_text)
    prog_bar = shared_progbar(num_scan);
  end

	img = zeros(num_scan, prod(hdr(1).dim));
  
	for v = 1:num_scan
		tmp = spm_read_vols(hdr(v));
    img(v,:) = tmp(:);
    if cfg.gui
		  spm_progress_bar('Set', round(100*v/num_scan));
    else
      prog_bar.update(v)
    end
	end
  
  if not(cfg.gui)
    prog_bar.stop;
  end

	img = img(:, coord);

end
