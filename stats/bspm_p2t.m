function t = bspm_p2t(alpha, df)
% BOB_P2T Get t-value from p-value + df
%
%   USAGE: t = bspm_p2t(alpha, df)
%       
%   OUTPUT
%       t = crtical t-value
%
%   ARGUMENTS
%       alpha = p-value
%       df = degrees of freedom
%
% =========================================
if nargin<2, disp('USAGE: bspm_p2t(p, df)'); return, end
t = tinv(1-alpha, df);


 
 
 
 
