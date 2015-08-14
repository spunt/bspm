function dRet = ipctb_P2TSub(dT)
% Converts two sided pval threshold to corresponding T-score threshold
% This function demands global values of pval and nSubjects to run
% pval is the confidence interval of the tails in the T-distribution
% nSubjects are used to calculate the degrees of freedom
global gipctb_dPVal gipctb_nDegFree;

dP = 2*(1-ipctb_spm_Tcdf(dT, gipctb_nDegFree));
dRet = abs(dP-gipctb_dPVal);