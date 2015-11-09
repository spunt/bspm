function tbx_gui_Masking
% tbx_gui_Masking - start SPM5 or SPM8 batch GUI with Masking toolbox.
% Called when Masking is selected from the toolbox pull-down menu.
%
% See also: opt_thresh, make_majority_mask

% spm.m puts the first (spm_select 'list'ed) file matching the pattern
% toolbox/masking/*_masking.m (this file) as the call-back for the toolbox
% name in the toolboxes pulldown menu, calling it with no output arguments.
% The following code handles this callback for SPM5 or SPM8(b)

if strcmp(spm('Ver'), 'SPM5')
    spm_jobman('interactive', '', 'jobs.tools.masking');
elseif strcmp(spm('Ver'), 'SPM8b') || strcmp(spm('Ver'), 'SPM8')
    % can't add non-exec branch jobs.tools.masking, hence add children
    job{1}.spm.tools.masking{1}.makeavg = struct;
    job{2}.spm.tools.masking{1}.optthr  = struct;
    spm_jobman('interactive', job);
end

% Note that "masking{1}." could be just "masking." if I change cfg_repeat
% to cfg_choice, but then it would also be necessary to remove the ".num="
% field from the auto-generated tbx_cfg_masking
