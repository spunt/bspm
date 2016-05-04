function z = bspm_t2z(t, df)
% BOB_T2P Get p-value from t-value + df
%
%   USAGE: p = bspm_t2z(t, df)
%       
%   OUTPUT
%       z = z-value
%
%   ARGUMENTS
%       t = t-value
%       df = degrees of freedom
%
% =========================================
if nargin<1, mfile_showhelp; return; end
z = -sqrt(2) * erfcinv((tcdf(t, df))*2);

 
 
 
 
