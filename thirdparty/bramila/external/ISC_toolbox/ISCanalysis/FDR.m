function [pID,pN] = FDR(p,q)
  % FORMAT pt = FDR(p,q)
  %
  % p   - vector of p-values
  % q   - False Discovery Rate level
  %
  % pID - p-value threshold based on independence or positive dependence
  % pN  - Nonparametric p-value threshold
  %______________________________________________________________________________
  % @(#)FDR.m     1.3 Tom Nichols 02/01/18
  
  
  p = sort(p(:));
  V = length(p);
  I = (1:V)';
  %I = (1:length(p))';
  
  cVID = 1;
  cVN = sum(1./(1:V));
  %cVN = sum(1./(1:length(p)));
  
  pID = p(max(find(p<=I/V*q/cVID)));
  pN = p(max(find(p<=I/V*q/cVN)));
  
  if isempty(pID)
    pID = NaN;
  end
  if isempty(pN)
    pN = NaN;
  end
  
