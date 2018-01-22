function cmd = bfsl_realign_and_unwarp(epi4d, se_neg, se_pos, paramfile)
    %  USAGE: cmd = bfsl_realign_and_unwarp(epi4d, se_neg, se_pos, paramfile)
    %
    if nargin < 4, mfile_showhelp; return; end
    if iscell(epi4d), epi4d = char(epi4d); end
    if iscell(se_neg), se_neg = char(se_neg); end
    if iscell(se_pos), se_pos = char(se_pos); end
    if iscell(paramfile), paramfile = char(paramfile); end
    [epidir,epibasename,e] = fileparts(epi4d);

    % | In & Out Filenames (Full Path)
    se_negpos = fullfile(epidir, 'se_negpos');
    topup_output = fullfile(epidir, 'topup_output');
    uref = fullfile(epidir, 'u_se_negpos')
    repi4d = fullfile(epidir, strcat('r_', epibasename))
    urepi4d = fullfile(epidir, strcat('ur_', epibasename))
    burepi4d = fullfile(epidir, strcat('bur_', epibasename))
    mean_burepi4d = fullfile(epidir, strcat('mean_bur_', epibasename));

    cmd = cell(6, 1);
    
    % merge negative & positive spin-echo volumes
    cmd{1} = sprintf('fslmerge -t %s %s %s', se_negpos, se_neg, se_pos);

    % run topup
    cmd{2} = sprintf('topup --imain=%s --datain=%s --config=b02b0.cnf --out=%s --iout=%s', se_negpos, paramfile, topup_output, uref);

    % motion correct epi, then register to unwarped reference
    cmd{3} = sprintf('mcflirt -in %s -out %s -plots -reffile %s', epi4d, repi4d, uref);

    % apply topup to epi
    cmd{4} = sprintf('applytopup --imain=%s --datain=%s --inindex=2 --method=jac --topup=%s --out=%s', repi4d, paramfile, topup_output, urepi4d);

    % extract brain from unwarped epi
    cmd{5} = sprintf('bet %s %s -F', urepi4d, burepi4d);

    % compute mean of realigned & unwarped
    cmd{6} = sprintf('fslmaths %s -Tmean %s', burepi4d, mean_burepi4d);

end



