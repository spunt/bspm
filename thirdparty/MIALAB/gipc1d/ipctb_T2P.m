function dPVal = ipctb_T2P(dTVal, nDegFree)
% Converts T-score threshold to corresponding two sided pval threshold
% This function requires global values to run
% pval is the confidence interval of the tails in the T-distribution
global gipctb_dPVal gipctb_nDegFree;

dPVal = 2*(1-ipctb_spm_Tcdf(dTVal, nDegFree));