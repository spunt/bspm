function generateLUT(Params,nrSession,winOn)

% This function generates look-up tables (LUTs)
% based on null distribution.
% 
% inputs:
% Params - struct containing all necessary parameters
% nrSession - session index
% winOn - perform calculations for windowed (winOn=1) or across-session data (winOn=0)
%
% See also:
%
% COMBINENULLVALS
% CALCULATETHRESHOLDS
% EVALUATEPVALS
% INITPARAMS
% RUNANALYSIS
%
% last modified 10.11.2010
% Jukka-Pekka Kauppi
% Tampere University of Technology

Pub = Params.PublicParams;
Priv = Params.PrivateParams;
q = Priv.q;
intVals = Priv.intVals;
acc = Priv.acc;

% load LUT general parameters:
if winOn
  load([Priv.statsDestination Priv.prefixLUT 'generalwin' Priv.prefixSession num2str(nrSession)])
else
  load([Priv.statsDestination Priv.prefixLUT 'general' Priv.prefixSession num2str(nrSession)])
end
% load the most important part (>=q) of the null distribution:
if winOn
  load([Priv.statsDestination 'permTestValswin' Priv.prefixSession num2str(nrSession)])
else
  load([Priv.statsDestination 'permTestVals' Priv.prefixSession num2str(nrSession)])
end
plotFigs = 0; % plot figures (on/off)

% get decimation factors for each LUT. Accuracy for each LUT can be changed
% by modifying "acc" field in PrivateParams -struct.
acc = floor(N*diff(intVals)./acc);
acc(acc<1)=1;

% set maximal accuracy:
acc(acc>=1)=1;

% find critical value according to:
% Nichols and Holmes: "Nonparametric Permutation Tests For 
% Functional Neuroimaging: A Primer with Examples", HBM 15, 2001.
critVal = nullVals(1+floor(N*q));
critValFis = nullValsFis(1+floor(NFis*q));

if ~winOn && Pub.calcPhase
    critValPh = nullValsPh(1+floor(NPh*q));
end
% generate look-up tables:
inds = [1 1+floor(intVals(2:end)*N)];
limits = nullVals(inds)
for k = 1:length(inds)-1
    v = nullVals(inds(k):inds(k+1)-1);
    v = v(1:acc(k):length(v));
    % calculate empirical cumulative cdf in a specified p-value range:
    %[notUsed XXXX] = ecdf(v);
    X = [min(v) unique(v)];
    X = X(:);
   % startV(k) = v(1);
   F = linspace(intVals(k),intVals(k+1),length(X));

   if winOn
     save([Priv.statsDestination Priv.prefixLUT num2str(k) 'win' Priv.prefixSession num2str(nrSession)],'X','F')
   else
     save([Priv.statsDestination Priv.prefixLUT num2str(k) Priv.prefixSession num2str(nrSession)],'X','F')
   end   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VISUALIZE LOOK-UP TABLES:
    if plotFigs
        if k == 1;figure;end
        subplot(ceil(sqrt(length(inds)-1)),ceil(sqrt(length(inds)-1)),k);
        plot(X,F,'.b:');ylim([intVals(k),intVals(k+1)]);
        
        set(gca,'YTickLabelMode','manual','YTickMode','manual',...
            'YTickLabel',intVals(k) + (intVals(k+1) - get(gca,'YTick')))
        xlabel(['ISC (total number of samples: ' num2str(length(X)) ')'])
        ylabel('p-value')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% generate look-up tables for Fisher transformed data:
inds = [1 1+floor(intVals(2:end)*NFis)];
%length(nullVals);
limitsFis = nullValsFis(inds)
for k = 1:length(inds)-1
    v = nullValsFis(inds(k):inds(k+1)-1);
    v = v(1:acc(k):length(v));
    % calculate empirical cumulative cdf in a specified p-value range:
    %[notUsed XFis] = ecdf(v);
    XFis = [min(v) unique(v)];
    XFis = XFis(:);
   % startV(k) = v(1);
   FFis = linspace(intVals(k),intVals(k+1),length(XFis));
   
   if winOn
     load([Priv.statsDestination Priv.prefixLUT num2str(k) 'win' Priv.prefixSession num2str(nrSession)])
     save([Priv.statsDestination Priv.prefixLUT num2str(k) 'win' Priv.prefixSession num2str(nrSession)],'XFis','FFis','X','F')
   else
     load([Priv.statsDestination Priv.prefixLUT num2str(k) Priv.prefixSession num2str(nrSession)])
     save([Priv.statsDestination Priv.prefixLUT num2str(k) Priv.prefixSession num2str(nrSession)],'XFis','FFis','X','F')
   end
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VISUALIZE LOOK-UP TABLES:
    if plotFigs
        if k == 1;figure;end
        subplot(ceil(sqrt(length(inds)-1)),ceil(sqrt(length(inds)-1)),k);
        plot(XFis,FFis,'.b-');ylim([intVals(k),intVals(k+1)]);
        set(gca,'YTickLabelMode','manual','YTickMode','manual',...
            'YTickLabel',intVals(k) + (intVals(k+1) - get(gca,'YTick')))
        xlabel(['Samples: ' num2str(length(XFis))])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

if ~winOn && Pub.calcPhase
% generate look-up tables for Fisher transformed data:
inds = [1 1+floor(intVals(2:end)*NPh)];
%length(nullVals);
limitsPh = nullValsPh(inds)
for k = 1:length(inds)-1
    v = nullValsPh(inds(k):inds(k+1)-1);
    v = v(1:acc(k):length(v));
    % calculate empirical cumulative cdf in a specified p-value range:
    %[notUsed XFis] = ecdf(v);
    XPh = [min(v) unique(v)];
    XPh = XPh(:);
   % startV(k) = v(1);
   FPh = linspace(intVals(k),intVals(k+1),length(XPh));
   
   load([Priv.statsDestination Priv.prefixLUT num2str(k) Priv.prefixSession num2str(nrSession)])
   save([Priv.statsDestination Priv.prefixLUT num2str(k) Priv.prefixSession num2str(nrSession)],'XFis','FFis','X','F','XPh','FPh')
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VISUALIZE LOOK-UP TABLES:
    if plotFigs
        if k == 1;figure;end
        subplot(ceil(sqrt(length(inds)-1)),ceil(sqrt(length(inds)-1)),k);
        plot(XPh,FPh,'.b-');ylim([intVals(k),intVals(k+1)]);
        set(gca,'YTickLabelMode','manual','YTickMode','manual',...
            'YTickLabel',intVals(k) + (intVals(k+1) - get(gca,'YTick')))
        xlabel(['Samples: ' num2str(length(XPh))])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
end

% save LUT general parameters:
if winOn
    save([Priv.statsDestination Priv.prefixLUT 'generalwin' Priv.prefixSession num2str(nrSession)]...
    ,'limits','intVals','q','critVal','N','critValFis','NFis','limitsFis')
else
    if Pub.calcPhase 
        save([Priv.statsDestination Priv.prefixLUT 'general' Priv.prefixSession num2str(nrSession)],'limits',...
        'intVals','q','critVal','N','critValFis','NFis','limitsFis','critValPh','NPh','limitsPh')
    else
        save([Priv.statsDestination Priv.prefixLUT 'general' Priv.prefixSession num2str(nrSession)],'limits',...
        'intVals','q','critVal','N','critValFis','NFis','limitsFis')
    end
end
