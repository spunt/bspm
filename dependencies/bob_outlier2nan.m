function y      = bob_outlier2nan_new(x, sd, rmcase)
% BOB_OUTLIER2NAN
% 
% Find outliers recusviely using leave-one-out method and convert to NaN.
% If input X is a matrix, outliers computed columnwise. 
% 
%   USAGE: y = bob_outlier2nan(x, sd, rmcase)
% _________________________________________________________________________
% INPUTS
%   x       = input matrix
%   sd      = number of standard deviations to truncate
%   rmcase  = option to remove case/row 
%
    if nargin < 1, mfile_showhelp; return; end
    if nargin < 2, sd = 3; end
    if nargin < 3, rmcase = 0; end
    
    y       = zeros(size(x));
    ncol    = size(x,2);
    idx2y   = 1:ncol; 
    
    while ~isempty(x)

        nnan1                      = sum(isnan(x));
        x(oneoutzscore(x, 1) > sd) = NaN;
        nancol                     = (sum(isnan(x)) - nnan1)==0;
        y(:, idx2y(nancol))        = x(:, nancol);
        x(:, nancol)               = [];
        idx2y(nancol)              = [];
        
    end

    if rmcase, y(isnan(mean(y,2)),:) = []; end

end
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