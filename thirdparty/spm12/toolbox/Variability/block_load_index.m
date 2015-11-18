%_______________________________________________________________________
%
% Compute condition-wise index for image data
%_______________________________________________________________________
%
% Input
%
% fmri_spec | first-level design (struct)
%
%_______________________________________________________________________
%
% Output
%
% index | condition-wise scan index (3D cell array)
%
%_______________________________________________________________________
%
% This file is part of the Variability Toolbox for SPM
% Published by the Lifespan Neural Dynamics Group
% Provided under the GNU General Public License
% See LICENSE for terms and conditions

function index = block_load_index(fmri_spec)

  %%
  %% assuming that all sessions have the same conditions
  %%
  session = fmri_spec.sess;
	condition = {session(1).cond.name};

	for cnd = 1:numel(condition)
		num_scan_cond = 0;

		for sess = 1:numel(session)
			onset = session(sess).cond(cnd).onset;
			duration = session(sess).cond(cnd).duration;

			if strcmp(fmri_spec.timing.units, 'secs')
				[onset, duration] = secs_to_scans(onset, duration, fmri_spec.timing.RT);
			else
        %% expect the user to specify the first scan as zero
				onset = onset + 1;
				if not(positive_ints(onset) && positive_ints(duration))
					error('Onsets and duration must be positive integers.')
				end
			end

			%% if duration is defined only once extend it to a vector
			if numel(duration) == 1
				duration = ones(1, numel(onset)) * duration;
			end

			for blk = 1:numel(onset)
				index{cnd}{sess}{blk} = onset(blk) - 1 + [1:duration(blk)];
				num_scan_cond = num_scan_cond + numel(index{cnd}{sess}{blk});
			end
		end
	end

end

%%
%% convert onsets/durations from seconds to scans
%%
function [onset, duration] = secs_to_scans(secs_onset, secs_duration, RT)
	onset_valid = sum(mod(secs_onset, RT)) == 0;
	duration_valid = sum(mod(secs_duration, RT)) == 0;

	if not(onset_valid && duration_valid)
		message = 'Onsets and/or durations not a multiple of RT,\n';
		message = 'rounding result to nearest scan number.';
		warning('%s', message);
	end

	onset = round((secs_onset / RT) + 1);
	duration = round(secs_duration / RT);
end

%%
%% check if all entries of a vector are positive integers
%%
function valid = positive_ints(vector)
	valid = (sum(vector >= 1) == length(vector)) && (sum(mod(vector,1)) == 0);
end
