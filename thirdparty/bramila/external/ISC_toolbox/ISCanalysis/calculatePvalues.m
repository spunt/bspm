function pvals = calculatePvalues(vals,nullVals)

% This function returns p-values according to null distribution (generated
% by permutation testing). Returned p-values can then be used to perform,
% e.g. FDR correction. p-values are obtained via look-up tables from empirical distribution function. Several
% look-up tables are used to provide specific accuracy for different p-value 
% ranges.
%
% inputs: 
% vals - vector of true distribution values
% nullVals - vector of null distribution values
%
% outputs:
% pvals - p-values of the observations from true distribution. Values are 
% set to nan if they do not exceed the false discovery rate q. 

% Jukka-Pekka Kauppi 1.9.2009
% Tampere University of Technology


q = 0.05; % false discovery rate level
plotFigs = 1; % plot figures (on/off)

N = length(nullVals);
M = length(vals);

% Set look-up table intervals (from 0 to q):
intVals = [0 0.005 0.01 0.02 0.03 0.04 q];

% set approximate number of samples in look-up tables (accuracy):
acc = [floor(N*0.005) 200 750 100 100 100];




% get corresponding decimation factors:
acc = floor(N*diff(intVals)./acc);
acc(acc<1)=1;

% remove nan-values from null distribution and sort values:
nullVals = nullVals(:)';
nullVals(isnan(nullVals)) = [];
nullVals = sort(nullVals);
nullVals = fliplr(nullVals);
% find critical value according to:
% Nichols and Holmes: "Nonparametric Permutation Tests For 
% Functional Neuroimaging: A Primer with Examples", HBM 15, 2001.
critVal = nullVals(1+floor(N*q));
% remove nan-values and other non-interesting values from true distribution:
valsTmp = vals;
vals = vals(:)';
crapVals = isnan(vals) | vals < critVal;
vals(crapVals) = [];

% sort values into descending order and save original ordering:
[vals sinds] = sort(vals);
[trash sinds] = sort(sinds);
vals = fliplr(vals);

% generate look-up tables:
inds = [1 1+floor(intVals(2:end)*N)];
for k = 1:length(inds)-1
    v = nullVals(inds(k):inds(k+1));
    v = v(1:acc(k):length(v));
    % calculate empirical cumulative cdf in a specified p-value range:
%    [notUsed X{k}] = ecdf(v);
    X{k} = [min(v) unique(v)];
    X{k} = X{k}(:);
    startV(k) = v(end);
    F{k} = linspace(intVals(k),intVals(k+1),length(X{k}));
    if k == 1
        % to avoid nans in interpolation, set last point in the highest-end 
        % look-up table at least as large as in true distribution:
        X{k}(end+1) = max(max(vals),max(X{k}))+eps;
        F{k}(end+1) = F{k}(end);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VISUALIZE LOOK-UP TABLES:
    if plotFigs
        if k == 1;figure;end
        subplot(ceil(sqrt(length(inds)-1)),ceil(sqrt(length(inds)-1)),k);
        plot(X{k},F{k},'.b-');ylim([intVals(k),intVals(k+1)]);
        set(gca,'YTickLabelMode','manual','YTickMode','manual',...
            'YTickLabel',intVals(k) + (intVals(k+1) - get(gca,'YTick')))
        xlabel(['Samples: ' num2str(length(X{k}))])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

startV = [1+eps startV];

% calculate p-values from look-up tables using linear interpolation:
p = zeros(size(vals));
for k = 1:length(startV)-1
    idx = find(vals >= startV(k+1) & vals < startV(k));
    p(idx) = interp1(X{k}(2:end),F{k}(2:end),vals(idx));
    p(idx) = intVals(k) + (intVals(k+1) - p(idx));
end

% return p-values using original indexing:
pvals = NaN*ones(M,1);
p = fliplr(p);
p = p(sinds);
pvals(~crapVals) = p;
pvals = pvals(:)';

if plotFigs
    figure,plot(pvals(:),valsTmp(:),'.');grid on;xlabel('p-value');
    ylabel('observation');hold on;
    % n'th highest p-value should correspond to n'th lowest observation 
    % so the following command should produce the same result as above:
    plot(flipud(sort(pvals(:))),sort(valsTmp(:)),'ro');grid on;hold off
end
