function zx = oneoutzscore(x, returnas)
% ONEOUTZSCORE Perform columnwise leave-one-out zscoring
% 
% USAGE: zx = oneoutzscore(x, returnas)
% 
%   returnas: 0, signed values (default); 1, absolute values
%
if nargin<1, disp('USAGE: zx = oneoutzscore(x, returnas)'); return; end
if nargin<2, returnas = 0; end
if size(x,1)==1, x=x'; end
zx              = x; 
[nrow, ncol]    = size(x);
for c = 1:ncol
    cin         = repmat(x(:,c), 1, nrow);
    theoneout   = cin(logical(eye(nrow)))';
    theleftin   = reshape(cin(logical(~eye(nrow))),nrow-1,nrow);
    cz          = (theoneout-nanmean(theleftin))./nanstd(theleftin);
    zx(:,c)     = cz';
end
if returnas, zx = abs(zx); end
end