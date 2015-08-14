function dTVal = ipctb_P2T(dPVal, nDegFree)
% Converts two sided pval threshold to corresponding T-score threshold
% This function demands global values of pval and nSubjects to run
% pval is the confidence interval of the tails in the T-distribution
% nSubjects are used to calculate the degrees of freedom
global gipctb_dPVal gipctb_nDegFree;

gipctb_dPVal = dPVal; 
gipctb_nDegFree = nDegFree; 
dTVal = fminsearch('ipctb_P2TSub',1);
