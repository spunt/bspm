function evaluatePvals2(Params,vals,sessionNr,winOn,nrBand)

% This function obtains the following thresholds for the mean correlation maps (GISC maps):
%
% 1.  p < 0.05,  no multiple comparisons correction
% 2.  p < 0.05,  FDR corrected using independence or positive dependence assumption
% 3.  p < 0.05,  FDR corrected (no assumptions)
% 4.  p < 0.05,  Bonferroni corrected
% 5.  p < 0.005,  no multiple comparisons correction
% 6.  p < 0.005, FDR corrected using independence or positive dependence assumption
% 7.  p < 0.005, FDR corrected (no assumptions)
% 8.  p < 0.005, Bonferroni corrected
% 9.  p < 0.001,  no multiple comparisons correction
% 10. p < 0.001, FDR corrected using independence or positive dependence assumption
% 11. p < 0.001, FDR corrected (no assumptions)
% 12. p < 0.001, Bonferroni corrected
%
% inputs:
% Params - struct containing all necessary parameters  
% vals - vector of data values for which thresholds are calculated
% sessionNr - session index
% winOn - perform calculations for windowed (winOn=1) or across-session data (winOn=0)
% nrBand - frequency subband index
%
% See also:
% CALCULATETHRESHOLDS
% INITPARAMS
% RUNanalysis
% 
% 
% last modified 26.10.2009
% Jukka-Pekka Kauppi
% Tampere University of Technology

plotFigs = 0;

Pub = Params.PublicParams;
Priv = Params.PrivateParams;

% load general look-up table parameters:
if winOn
  load([Priv.statsDestination Priv.prefixLUT 'generalwin'])
else
  load([Priv.statsDestination Priv.prefixLUT 'general' Priv.prefixSession num2str(sessionNr) Priv.prefixFreqBand num2str(nrBand)])
end

% remove nan-values and redundant values from the data of interest:
valsTmp = vals;
M = length(vals);
vals = vals(:)';
crapVals = isnan(vals) | vals < critVal;
vals(crapVals) = [];

% sort values into descending order:
[vals sinds] = sort(vals);
[trash sinds] = sort(sinds); % sinds keeps track of original indexing
vals = fliplr(vals);

SS = 0;
% calculate p-values from look-up tables using linear interpolation:
p = zeros(size(vals));
for k = 1:length(limits)-1
    idx = find(vals >= limits(k+1) & vals < limits(k));
    if ~isempty(idx)
      % load appropriate look-up table:
      if winOn
        load([Priv.statsDestination Priv.prefixLUT num2str(k) 'win'])
      else
        load([Priv.statsDestination Priv.prefixLUT num2str(k) Priv.prefixSession num2str(sessionNr) Priv.prefixFreqBand num2str(nrBand)])
      end      
        if k == 1
            % to avoid nans in interpolation, set last point in the highest-end
            % look-up table at least as large as the highest value in true distribution:
            limits(1) = max(vals)+1e-6;
            X(end+1) = max(max(vals),max(X))+1e-6;
            F(end+1) = F(end);    
        end
        p(idx) = interp1(X(2:end),F(2:end),vals(idx));
        p(idx) = intVals(k) + (intVals(k+1) - p(idx));    
        SS = SS + length(idx);
    end
end

%SS
%sum(~crapVals)
%sum(crapVals)

% return p-values using original indexing:
pvals = NaN*ones(M,1);
p = fliplr(p);
p = p(sinds);

pvals(~crapVals) = p;
pvals = pvals(:)';

size(pvals)
size(valsTmp)

if plotFigs
    figure,plot(pvals(:),valsTmp(:),'.');grid on;xlabel('p-value');
    ylabel('observation');hold on;
end

% perform multiple comparison corrections:
[pID,pN] = FDR(pvals,Priv.q);
[pID2,pN2] = FDR(pvals,(Priv.q)/10);
[pID3,pN3] = FDR(pvals,(Priv.q)/(50));

% corrected thresholds:
pvals_Th = [pID pN Priv.q/length(pvals) Priv.q/10 pID2 pN2 ((Priv.q)/10)/length(pvals) Priv.q/50 pID3 pN3 ((Priv.q)/50)/length(pvals)];

Th = zeros(1,11);

% use look-up tables to find threshold values:
for k = 1:length(limits)-1
  idx = find(pvals_Th < intVals(k+1) & pvals_Th >= intVals(k));
  if ~isempty(idx)
    % load appropriate look-up table:
    if winOn
      load([Priv.statsDestination Priv.prefixLUT num2str(k) 'win'])
    else    
      load([Priv.statsDestination Priv.prefixLUT num2str(k) Priv.prefixSession num2str(sessionNr)])
    end
    if k == 1
      % to avoid nans in interpolation, set last point in the highest-end
      % look-up table at least as large as the highest value in true distribution:
      X(end+1) = max(max(vals),max(X))+1e-6;
      F(end+1) = 0;
    end
    pvals_Th_flipped = intVals(k) + (intVals(k+1) - pvals_Th(idx));
    Th(idx) = interp1(F(2:end),X(2:end),pvals_Th_flipped);
  end
end

Th = [critVal Th]
pvals_Th = [Priv.q pvals_Th]

Th_info{1} = ['none_' num2str(Priv.q)];
Th_info{2} = ['pID_' num2str(Priv.q)];
Th_info{3} = ['pN_' num2str(Priv.q)];
Th_info{4} = ['bonf_' num2str(Priv.q)];
Th_info{5} = ['none_' num2str(Priv.q/10)];
Th_info{6} = ['pID_' num2str(Priv.q/10)];
Th_info{7} = ['pN_' num2str(Priv.q/10)];
Th_info{8} = ['bonf_' num2str(Priv.q/10)];
Th_info{9} = ['none_' num2str(Priv.q/50)];
Th_info{10} = ['pID_' num2str(Priv.q/50)];
Th_info{11} = ['pN_' num2str(Priv.q/50)];
Th_info{12} = ['bonf_' num2str(Priv.q/50)];

Th_info

% save thresholds:
save([Priv.statsDestination 'Th' Priv.prefixFreqBand num2str(nrBand) Priv.prefixSession... 
num2str(sessionNr) 'win' num2str(winOn)],'Th','Th_info','pvals_Th')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction FDR by Tom Nichols:

function [pID,pN] = FDR(p,q)
% FORMAT pt = FDR(p,q)
% 
% p   - vector of p-values
% q   - False Discovery Rate level
%
% pID - p-value threshold based on independence or positive dependence
% pN  - Nonparametric p-value threshold
%______________________________________________________________________________
% @(#)FDR.m	1.3 Tom Nichols 02/01/18

p = sort(p(:));
V = length(p); 
I = (1:V)';
cVID = 1;
cVN = sum(1./(1:V));

pID = p(max(find(p<=I/V*q/cVID)));
pN = p(max(find(p<=I/V*q/cVN)));

if isempty(pID)
  pID = NaN;
end
if isempty(pN)
  pN = NaN;
end
