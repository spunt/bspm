function [a, R2] = ipctb_ica_regress(y, X)
% Inputs are the observed data and the model matrix (one or more models)
% Assuming the rows are greater in number than columns.

a = pinv(X)*y;

% calculating R-square statistic
if nargout > 1
    residual = y - X*a;
    R2 = 1 - (norm(residual, 2) / norm(y - mean(y), 2))^2;
end