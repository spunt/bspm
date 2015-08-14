function [p_value, t_value, df] = icatb_ttest(x, tail)
%% Compute one sample t-test of sample x
%
% Inputs:
% 1. x - Sample
% 2. tail - Options are 0 (two_tailed), 1 (right tailed) or -1 (left
% tailed)
%
% Outputs:
% 1. p_value - P value
% 2. t_value - T value
% 3. df - Degrees of freedom

% Convert x to column vector
x = x(:);

%% Remove NaN's
x(isnan(x)) = [];

if isempty(x)
    error('Data is empty or it has all NaN''s');
end

if ~exist('tail', 'var')
    tail = 0;
end

% Degrees of freedom
df = length(x) - 1;

%% T-value
t_value = mean(x) / (std(x) / sqrt(length(x)));

%% P value
p_value = icatb_get_pvalue(t_value, df, tail);

