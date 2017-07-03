function evaluatePvals(Params,vals,sessionNr,winOn,nrBand,mapType)

% This function obtains the following thresholds for the mean correlation maps (ISC maps):
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
% mapType: 
% 1 = mean ISC + permutation 
% 2 = t-stats 
% 3 = phase synhcronization + permutation
%
% See also:
% CALCULATETHRESHOLDS
% INITPARAMS
% RUNANALYSIS
%
%
% last modified 31.10.2013
% Jukka-Pekka Kauppi
% Tampere University of Technology/University of Helsinki

if nargin == 5
    mapType = 1;
end

plotFigs = 0;

Pub = Params.PublicParams;
Priv = Params.PrivateParams;

switch mapType
    % load general look-up table parameters:
    case 1
        if winOn
            load([Priv.statsDestination Priv.prefixLUT 'generalwin' Priv.prefixSession num2str(sessionNr)])
        else
            load([Priv.statsDestination Priv.prefixLUT 'general' Priv.prefixSession num2str(sessionNr)])
        end
    case 2
%        load([Priv.statsDestination Priv.prefixLUT 'general' Priv.prefixSession num2str(sessionNr)])
%        critVal = critValFis;
%        N = NFis;
%        limits = limitsFis;
    case 3
        load([Priv.statsDestination Priv.prefixLUT 'general' Priv.prefixSession num2str(sessionNr)])
        critVal = critValPh;
        N = NPh;
        limits = limitsPh;
end

Th_info{1} = ['none_' num2str(Priv.q)];
Th_info{2} = ['pID_' num2str(Priv.q)];
Th_info{3} = ['pN_' num2str(Priv.q)];
Th_info{4} = ['bonf_' num2str(Priv.q)];
Th_info{5} = ['none_' num2str(Priv.q/5)];
Th_info{6} = ['pID_' num2str(Priv.q/5)];
Th_info{7} = ['pN_' num2str(Priv.q/5)];
Th_info{8} = ['bonf_' num2str(Priv.q/5)];
Th_info{9} = ['none_' num2str(Priv.q/50)];
Th_info{10} = ['pID_' num2str(Priv.q/50)];
Th_info{11} = ['pN_' num2str(Priv.q/50)];
Th_info{12} = ['bonf_' num2str(Priv.q/50)];

if mapType ~= 2
% remove nan-values and redundant values from the data of interest:
valsTmp = vals;
M = length(vals);
vals = vals(:)';
crapVals = isnan(vals) | vals < critVal;
vals(crapVals) = [];
% check empty condition for vals (no any voxels showing p < 0.05):
if isempty(vals)
    Th = zeros(1,12);
    pvals_Th = NaN*ones(1,12);
    % save thresholds:
    switch mapType
        case 1
            save([Priv.statsDestination 'Th' Priv.prefixFreqBand num2str(nrBand) Priv.prefixSession...
                num2str(sessionNr) 'win' num2str(winOn)],'Th','pvals_Th')
        case 2
            save([Priv.statsDestination 'ThT' Priv.prefixFreqBand num2str(nrBand) Priv.prefixSession...
                num2str(sessionNr) 'win' num2str(winOn)],'Th','pvals_Th')
        case 3
            save([Priv.statsDestination 'ThPhase' Priv.prefixFreqBand num2str(nrBand) Priv.prefixSession...
                num2str(sessionNr) 'win' num2str(winOn)],'Th','pvals_Th')
    end
    return
end

% Sort values into descending order:
[vals sinds] = sort(vals);
[trash sinds] = sort(sinds); % sinds keeps track of original indexing
vals = fliplr(vals);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get look-up table containing p-values and statistic values 
% from permutation distribution:
XX=[];
FF=[];
for k = 1:length(limits)-1
    switch mapType
        case 1
            if winOn
                load([Priv.statsDestination Priv.prefixLUT num2str(k) 'win' Priv.prefixSession num2str(sessionNr)])
            else
                load([Priv.statsDestination Priv.prefixLUT num2str(k) Priv.prefixSession num2str(sessionNr)])
            end
            XX = [X' XX];
            FF = [FF F];
        case 2
            load([Priv.statsDestination Priv.prefixLUT num2str(k) Priv.prefixSession num2str(sessionNr)])
            XX = [XFis' XX];
            FF = [FF FFis];
        case 3
            load([Priv.statsDestination Priv.prefixLUT num2str(k) Priv.prefixSession num2str(sessionNr)])
            XX = [XPh' XX];
            FF = [FF FPh];
    end
end
X = XX;
F = fliplr(FF);
ss = find(diff(X)<=0);
X(ss) = [];
F(ss) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the end of the range higher than the highest observed sample value 
% to avoid extrapolation (and thus NaNs):
X(end) = max(max(vals),max(X))+1e-6;
X = [critVal X];
F = [0.05 F];

% Interpolation to obtain p-values for observed values:  
p = interp1(X(1:end),F(1:end),vals);

% perform multiple comparison corrections:
[pID,pN] = FDR(p,Priv.q);
[pID2,pN2] = FDR(p,Priv.q/5);
[pID3,pN3] = FDR(p,Priv.q/50);

% Set critical p-values (alpha-levels):
NN = length(valsTmp(~isnan(valsTmp)));
pvals_Th = [pID pN Priv.q/NN Priv.q/5 pID2 pN2 (Priv.q/5)/NN Priv.q/50 pID3 pN3 (Priv.q/50)/NN];

% Interpolation to obtain the corresponding critical thresholds from 
% critical p-values:
Th = interp1(F(1:end),X(1:end),pvals_Th);

% Save also the "loosest" uncorrected threshold whos critical value was
% computed earlier:
Th = [critVal Th];
pvals_Th = [Priv.q pvals_Th];

% save critical thresholds and p-values:
switch mapType
    case 1
        save([Priv.statsDestination 'Th' Priv.prefixFreqBand num2str(nrBand) Priv.prefixSession...
            num2str(sessionNr) 'win' num2str(winOn)],'Th','Th_info','pvals_Th')
    case 2
        save([Priv.statsDestination 'ThFis' Priv.prefixFreqBand num2str(nrBand) Priv.prefixSession...
            num2str(sessionNr) 'win' num2str(winOn)],'Th','Th_info','pvals_Th')
    case 3
        save([Priv.statsDestination 'ThPhase' Priv.prefixFreqBand num2str(nrBand) Priv.prefixSession...
            num2str(sessionNr)],'Th','Th_info','pvals_Th')
end

% compute t-map thresholds:
else
    V = Priv.nrSubjects*(Priv.nrSubjects-1)/2 -1; % DOF
    if V == 0
        vals = sqrt(Priv.dataSize(sessionNr,end)-3)*vals;        
        valsTmp = vals;
        vals = vals(:)';
        critVal = norminv(0.95,0,1);
        crapVals = isnan(vals) | vals < critVal;
        vals(crapVals) = [];
        M = length(vals);
        pvals = 1-normcdf(vals,0,1);
    else
        critVal = tinv(0.95,V);
        valsTmp = vals;
        vals = vals(:)';
        crapVals = isnan(vals) | vals < critVal;
        vals(crapVals) = [];
        M = length(vals);
        pvals = 1-tcdf(vals,V);
    end
    if plotFigs
        figure,plot(pvals(:),vals(:),'.');grid on;xlabel('p-value');
        ylabel('observation');hold on;
    end
    % multiple comparisons correction:
    [pID,pN] = FDR(pvals,Priv.q);
    [pID2,pN2] = FDR(pvals,Priv.q/5);
    [pID3,pN3] = FDR(pvals,Priv.q/50);
    % compute and save critical p-vals and thresholds:
    pvals_Th = [Priv.q pID pN Priv.q/length(pvals) Priv.q/5 pID2 pN2 (Priv.q/5)/length(pvals) Priv.q/50 pID3 pN3 (Priv.q/50)/length(pvals)];
    if V == 0
        Th = norminv(1-pvals_Th,0,1);
        Th = Th/sqrt(Priv.dataSize(sessionNr,end)-3);
    else
        Th = tinv(1-pvals_Th,V);
    end
    save([Priv.statsDestination 'ThT' Priv.prefixFreqBand num2str(nrBand) Priv.prefixSession...
        num2str(sessionNr) 'win' num2str(winOn)],'Th','Th_info','pvals_Th')
end

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
