function combineNullVals(Params,nrSession,winOn)

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
% RUNANALYSIS


% last modified 6.10.2010
% Jukka-Pekka Kauppi
% Tampere University of Technology

Pub = Params.PublicParams;
Priv = Params.PrivateParams;

% load null-distribution values:
nullVals = [];
nullValsFis = [];
if ~winOn && Pub.calcPhase
    nullValsPh = [];
end

for permSetIdx = 1:Priv.nrPermutationSets
  if winOn %%% KORJAA!!
    load([Priv.statsDestination 'vals' num2str(permSetIdx) Priv.prefixSession num2str(nrSession) 'win'])
  else
    load([Priv.statsDestination 'vals' num2str(permSetIdx) Priv.prefixSession num2str(nrSession)])
  end
  nullVals = [nullVals val];
  nullValsFis = [nullValsFis valFish];
  clear val valFish
  if ~winOn && Pub.calcPhase
      nullValsPh = [nullValsPh ph];    
  else
      nullValsPh = []; 
  end
end

% remove nan-values from null distribution if any:
nullVals = nullVals(:)';
nullVals(isnan(nullVals)) = [];
nullValsFis = nullValsFis(:)';
nullValsFis(isnan(nullValsFis)) = [];
if ~winOn && Pub.calcPhase
    nullValsPh = nullValsPh(:)';
    nullValsPh(isnan(nullValsPh)) = [];
end
% save some data for visualization purposes:
if length(nullVals) >= 100000
  rP = randperm(length(nullVals));
  nullVis = nullVals(rP(1:100000));
  nullFisVis = nullValsFis(rP(1:100000));
if ~winOn && Pub.calcPhase
  nullPhVis = nullValsPh(rP(1:100000));
end
else
  nullVis = nullVals;
  nullFisVis = nullValsFis;
  if ~winOn && Pub.calcPhase
      nullPhVis = nullValsPh;
  end
end
if winOn
  save([Priv.statsDestination 'nullViswin' Priv.prefixSession num2str(nrSession)],'nullVis','nullFisVis')
else
    if Pub.calcPhase
        save([Priv.statsDestination 'nullVis' Priv.prefixSession num2str(nrSession)],'nullVis','nullFisVis','nullPhVis')  
    else
        save([Priv.statsDestination 'nullVis' Priv.prefixSession num2str(nrSession)],'nullVis','nullFisVis')  
    end
end

% sort values into descending order:
nullVals = sort(nullVals);
nullVals = fliplr(nullVals);
nullValsFis = sort(nullValsFis);
nullValsFis = fliplr(nullValsFis);
if ~winOn && Pub.calcPhase
    nullValsPh = sort(nullValsPh);
    nullValsPh = fliplr(nullValsPh);    
end

N = length(nullVals);
q = Params.PrivateParams.q;
nullVals = nullVals(1:(1+floor(N*q)));

NFis = length(nullValsFis);
nullValsFis = nullValsFis(1:(1+floor(NFis*q)));

if ~winOn && Pub.calcPhase
    NPh = length(nullValsPh);
    if ~isempty(nullValsPh)
        nullValsPh = nullValsPh(1:(1+floor(NPh*q)));    
    end
end

% save the most important part (>=q) of the null distribution:
if winOn
  save([Priv.statsDestination 'permTestValswin' Priv.prefixSession num2str(nrSession)],'nullVals','nullValsFis')
else
    if Pub.calcPhase
        save([Priv.statsDestination 'permTestVals' Priv.prefixSession num2str(nrSession)],'nullVals','nullValsFis','nullValsPh')        
    else
        save([Priv.statsDestination 'permTestVals' Priv.prefixSession num2str(nrSession)],'nullVals','nullValsFis')
    end
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
  save([Priv.statsDestination Priv.prefixLUT 'generalwin' Priv.prefixSession num2str(nrSession)],'N','NFis')
else
    if Pub.calcPhase
        save([Priv.statsDestination Priv.prefixLUT 'general' Priv.prefixSession num2str(nrSession)],'N','NFis','NPh') 
    else
        save([Priv.statsDestination Priv.prefixLUT 'general' Priv.prefixSession num2str(nrSession)],'N','NFis')
    end
end
