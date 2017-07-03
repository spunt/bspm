function [PCALoadings,PCAScores,var_expl,mse_data,yfitPCR] = bramila_pcaregress(X,Y,ncomp,folds,mcreps)

% perform PCA regression with k-fold cross validation. Results are
% comparable with thos given by PLSREGRESS

warning('off','MATLAB:rankDeficientMatrix');

[PCALoadings,PCAScores,latent,tsquared,PCAVar,mu] = pca(X,'Economy',false);

var_expl=PCAVar';

Y0 = bsxfun(@minus, Y, mean(Y,1));

const = ones(size(X,1),1);

for ord=1:size(PCAScores,2)
    yfitPCR=nan(size(Y0));
    for i=1:size(Y0,2)
        y=Y0(:,i);
        my=mean(y);
        betaPCR = regress(y-my, PCAScores(:,1:ord));
        betaPCR = PCALoadings(:,1:ord)*betaPCR;
        betaPCR = [my - mean(X)*betaPCR; betaPCR];
        yfitPCR(:,i) = [const,X]*betaPCR;
    end
    var_expl(2,ord) = sum(sum(abs(yfitPCR).^2,1))/sum(sum(abs(Y0).^2,1));
end

if any(var_expl(2,:)>1+10*eps)
    error('Bad variance calculation (over 100%%)!')
end

fun = @(Xtrain,ytrain,Xtest,ytest,maxNumComp) pcrsse(Xtrain,ytrain,Xtest,ytest,ncomp);
sumsqerr = crossval(fun,X,Y,'KFold',folds,'mcreps',mcreps);

data_x = 2*sumsqerr/mcreps;

mse(1:ncomp+1) = mean(data_x);

se = nan(size(mse));

se(1:ncomp+1) = std(data_x)/sqrt(size(data_x,1));

minMSE1 = min(mse);
minIx1 = find(mse==minMSE1,1);
minplus1 = mse(minIx1) + se(minIx1);
seIx1= find((mse(1:minIx1) <= minplus1),1,'first');

mse_data.MSE = mse;
mse_data.SE           = se;
mse_data.IndexMinMSE  = minIx1;
mse_data.Index1SE     = seIx1;

end

function sse = pcrsse(Xtrain,ytrain,Xtest,ytest,maxNumComp)
%PCRSSE SSE for Principal Components Regression cross-validation.
% This function is used in the demo PLSPCRDEMO.
%
% SSE = PCRSSE(XTRAIN,YTRAIN,XTEST,YTEST,NCOMP) returns a vector of sum of
% squared prediction errors for principal components regression models with
% 0:10 components, with XTRAIN and YTRAIN as training data, and XTEST and
% YTEST as test data.

%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/04/20 19:46:40 $

sse = zeros(1,maxNumComp+1);

% The 0'th model is just the mean of the training response data.
yfit0 = mean(ytrain);
a = bsxfun(@minus, ytest, yfit0);
sse(1) = sum(sum(abs(a).^2, 2));
%sse(1) = sum((ytest - yfit0).^2);

% Compute PCA loadings from the training predictor data, and regress the first
% 10 principal components on the centered traiing response data.
[Loadings,Scores] = pca(Xtrain,'Economy',false);
beta = zeros(maxNumComp,size(ytrain,2));
for i=1:size(ytrain,2)
    beta(:,i) = regress(ytrain(:,i)-yfit0(i), Scores(:,1:maxNumComp));
end

% Compute predictions for the 1st through 10th model.
const = ones(size(Xtest,1),1);
for ncomp = 1:maxNumComp
    yfit = zeros(size(ytest));
    for i=1:size(ytrain,2)
        beta0 = Loadings(:,1:ncomp)*beta(1:ncomp,i);
        beta1 = [yfit0(i) - mean(Xtrain)*beta0; beta0];
        yfit(:,i) = [const,Xtest]*beta1;
    end
    sse(ncomp+1) = sum(sum(abs(ytest - yfit).^2, 2));
end

end

