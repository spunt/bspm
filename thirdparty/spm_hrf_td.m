function td = spm_hrf_td(TR)
% returns the temporal derivative of the hemodynamic response function
% FORMAT td = spm_hrf_td(TR);
% TR = scan repeat time
%
% Pulled from code in spm5
%_______________________________________________________________________

% get canonical hemodynaic response function
%---------------------------------------------------------------
[bf p]         = spm_hrf(TR);

% compute time derivative
%---------------------------------------------------------------
dp     = 1;
p(6)   = p(6) + dp;
D      = (bf(:,1) - spm_hrf(TR,p))/dp;
bf     = [bf D(:)];
p(6)   = p(6) - dp;
td      = D(:);