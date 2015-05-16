function z = bspm_p2z(p)
% BOB_P2Z Get z-stat form p-value (t-dist)
%
%   USAGE: z = bspm_p2z(p)
%       
% =========================================
if nargin<1, error('USAGE: z = bspm_p2z(p)'); end
z = -sqrt(2) * erfcinv((1-p)*2);

 
 
 
 
