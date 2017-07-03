function [Xi2 p]= bramila_chisquare(x)
% usage
%	[Xi2 p] = bramila_chisquare(x)
% where:
%	x = a vector with ids of answers (usually integers e.g. -1 = no and 1 = yes)

%
% The chi-square compares frequencies obtained  in the sample to those expected according to the null hypothesis (i.e., no difference in the population). 
% see http://se.mathworks.com/matlabcentral/answers/96572-how-can-i-perform-a-chi-square-test-to-determine-how-statistically-different-two-proportions-are-in
% see http://www.upa.pdx.edu/IOA/newsom/da1/ho_z-test.pdf

uu=unique(x);
observed=histc(x,uu); % counts each item
expected = length(x)/length(uu);
Xi2=sum((observed-expected).^2 ./ expected);

p = 1 - chi2cdf(Xi2,1);
			
