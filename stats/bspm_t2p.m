function p = bspm_t2p(t, df)
% BOB_T2P Get p-value from t-value + df
%
%   USAGE: p = bspm_t2p(t, df)
%       
%   OUTPUT
%       p = p-value
%
%   ARGUMENTS
%       t = t-value
%       df = degrees of freedom
%
% =========================================
if nargin<2, disp('USAGE: bspm_t2p(p, df)'); return, end
p = tcdf(t, df);
p = 1 - p;


 
 
 
 
