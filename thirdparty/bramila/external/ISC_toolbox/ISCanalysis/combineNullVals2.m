function combineNullVals2(Params,nrBand,nrSession,winOn)

% This function combines results of the permutation test that has been run
% in several batches. All samples of the right-tail (p < 0.05) are 
% saved (permTestVals.mat,permTestValswin.mat). Also selected samples across
% the entire null-distribution are saved for visualization purposes 
% (nullVis.mat, nullViswin.mat).
% 
% inputs:
% Params - struct containing all necessary parameters
% nrSession - session index
% winOn - perform calculations for windowed (winOn=1) or across-session data (winOn=0)
%
% See also:
%
% PERMUTATIONTEST 
% CALCULATETHRESHOLDS
% INITPARAMS
% RUNLSA
%
% last modified 2.11.2009
% Jukka-Pekka Kauppi
% Tampere University of Technology

Pub = Params.PublicParams;
Priv = Params.PrivateParams;

% load null-distribution values:
nullVals = [];
nullValsFis = [];
for permSetIdx = 1:Priv.nrPermutationSets
  if winOn
    load([Priv.statsDestination 'vals' num2str(permSetIdx) 'win'])
  else
    load([Priv.statsDestination 'vals' num2str(permSetIdx) Priv.prefixSession num2str(nrSession) Priv.prefixFreqBand num2str(nrBand)])
  end
  nullVals = [nullVals val];
  nullValsFis = [nullValsFis valFish];
  clear val valFish
end

% remove nan-values from null distribution if any:
nullVals = nullVals(:)';
nullVals(isnan(nullVals)) = [];
nullValsFis = nullValsFis(:)';
nullValsFis(isnan(nullValsFis)) = [];

% save some data for visualization purposes:
if length(nullVals) >= 100000
  rP = randperm(length(nullVals));
  nullVis = nullVals(rP(1:100000));
  nullFisVis = nullValsFis(rP(1:100000));
else
  nullVis = nullVals;
  nullFisVis = nullValsFis;
end
if winOn
  save([Priv.statsDestination 'nullViswin'],'nullVis','nullFisVis')
else
  save([Priv.statsDestination 'nullVis' Priv.prefixSession num2str(nrSession)],'nullVis','nullFisVis')
end

% sort values into descending order:
nullVals = sort(nullVals);
nullVals = fliplr(nullVals);
nullValsFis = sort(nullValsFis);
nullValsFis = fliplr(nullValsFis);

N = length(nullVals);
q = Params.PrivateParams.q;
nullVals = nullVals(1:(1+floor(N*q)));

NFis = length(nullValsFis);
q = Params.PrivateParams.q;
nullValsFis = nullValsFis(1:(1+floor(N*q)));

% save the most important part (>=q) of the null distribution:
if winOn
  save([Priv.statsDestination 'permTestValswin'],'nullVals','nullValsFis')
else
  save([Priv.statsDestination 'permTestVals' Priv.prefixSession num2str(nrSession) Priv.prefixFreqBand num2str(nrBand)],'nullVals','nullValsFis')
end

% delete unimportant part of the data from directory:
%for permSetIdx = 1:Pub.nrPermutationSets
%  if winOn
%    delete([Priv.statsDestination 'vals' num2str(permSetIdx) 'win.mat'])
%  else
%    delete([Priv.statsDestination 'vals' num2str(permSetIdx) Priv.prefixSession num2str(nrSession) '.mat'])
%  end
%end

% save LUT general parameters:
if winOn
  save([Priv.statsDestination Priv.prefixLUT 'generalwin'],'N','NFis')
else
  save([Priv.statsDestination Priv.prefixLUT 'general' Priv.prefixSession num2str(nrSession) Priv.prefixFreqBand num2str(nrBand)],'N','NFis')
end
