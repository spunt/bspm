function [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse_data,stats] = bramila_plsregress(X,Y,ncomp,varargin)

% PLSREGRESS of Matlab2012b with standard error computation for CV scheme
% (taken from LASSO function)

if nargin < 2
    error(message('stats:plsregress:TooFewInputs'));
end

[n,dx] = size(X);
ny = size(Y,1);
if ny ~= n
    error(message('stats:plsregress:SizeMismatch'));
end

% Return at most maxncomp PLS components
maxncomp = min(n-1,dx);
if nargin < 3
    ncomp = maxncomp;
elseif ~isscalar(ncomp) || ~isnumeric(ncomp) || (ncomp~=round(ncomp)) || (ncomp<=0)
    error(message('stats:plsregress:BadNcomp'));
elseif ncomp > maxncomp
    error(message('stats:plsregress:MaxComponents', maxncomp));
end

names = {'cv'                  'mcreps'            'options'};
dflts = {'resubstitution'        1                      []   };
[cvp,mcreps,ParOptions] = internal.stats.parseArgs(names, dflts, varargin{:});

if isnumeric(cvp) && isscalar(cvp) && (cvp==round(cvp)) && (0<cvp) && (cvp<=n)
    % ok, cvp is a kfold value. It will be passed as such to crossval.
elseif isequal(cvp,'resubstitution')
    % ok
elseif isa(cvp,'cvpartition')
    if strcmp(cvp.Type,'resubstitution')
        cvp = 'resubstitution';
    else
        % ok
    end
else
    error(message('stats:plsregress:InvalidCV'));
end

if ~(isnumeric(mcreps) && isscalar(mcreps) && (mcreps==round(mcreps)) && (0<mcreps))
    error(message('stats:plsregress:InvalidMCReps'));
elseif mcreps > 1 && isequal(cvp,'resubstitution')
    error(message('stats:plsregress:InvalidResubMCReps'));
end

% Center both predictors and response, and do PLS
meanX = mean(X,1);
meanY = mean(Y,1);
X0 = bsxfun(@minus, X, meanX);
Y0 = bsxfun(@minus, Y, meanY);

if nargout <= 2
    [Xloadings,Yloadings] = simpls(X0,Y0,ncomp);
    
elseif nargout <= 4
    [Xloadings,Yloadings,Xscores,Yscores] = simpls(X0,Y0,ncomp);
    
else
    % Compute the regression coefs, including intercept(s)
    [Xloadings,Yloadings,Xscores,Yscores,Weights] = simpls(X0,Y0,ncomp);
    beta = Weights*Yloadings';
    beta = [meanY - meanX*beta; beta];
    
    % Compute the percent of variance explained for X and Y
    if nargout > 5
        pctVar = [sum(abs(Xloadings).^2,1) ./ sum(sum(abs(X0).^2,1));
                  sum(abs(Yloadings).^2,1) ./ sum(sum(abs(Y0).^2,1))];
    end
    
    if nargout > 6
        if isequal(cvp,'resubstitution')
            % Compute MSE for models with 0:ncomp PLS components, by
            % resubstitution.  CROSSVAL can handle this, but don't waste time
            % fitting the whole model again.
            mse = zeros(2,ncomp+1,class(pctVar));
            mse(1,1) = sum(sum(abs(X0).^2, 2));
            mse(2,1) = sum(sum(abs(Y0).^2, 2));
            for i = 1:ncomp
                X0reconstructed = Xscores(:,1:i) * Xloadings(:,1:i)';
                Y0reconstructed = Xscores(:,1:i) * Yloadings(:,1:i)';
                mse(1,i+1) = sum(sum(abs(X0 - X0reconstructed).^2, 2));
                mse(2,i+1) = sum(sum(abs(Y0 - Y0reconstructed).^2, 2));
            end
            mse_data = mse / n;
            % We now have the reconstructed values for the full model to use in
            % the residual calculation below
        else
            % Compute MSE for models with 0:ncomp PLS components, by cross-validation
            mse_data = plscv(X,Y,ncomp,cvp,mcreps,ParOptions);
            if nargout > 7
                % Need these for the residual calculation below
                X0reconstructed = Xscores*Xloadings';
                Y0reconstructed = Xscores*Yloadings';
            end
        end
    end
    
    if nargout > 7
        % Save the PLS weights and compute the T^2 values.
        stats.W = Weights;
        stats.T2 = sum( bsxfun(@rdivide, abs(Xscores).^2, var(Xscores,[],1)) , 2);
        
        % Compute X and Y residuals
        stats.Xresiduals = X0 - X0reconstructed;
        stats.Yresiduals = Y0 - Y0reconstructed;
    end
end


%------------------------------------------------------------------------------
%SIMPLS Basic SIMPLS.  Performs no error checking.
function [Xloadings,Yloadings,Xscores,Yscores,Weights] = simpls(X0,Y0,ncomp)

[n,dx] = size(X0);
dy = size(Y0,2);

% Preallocate outputs
outClass = superiorfloat(X0,Y0);
Xloadings = zeros(dx,ncomp,outClass);
Yloadings = zeros(dy,ncomp,outClass);
if nargout > 2
    Xscores = zeros(n,ncomp,outClass);
    Yscores = zeros(n,ncomp,outClass);
    if nargout > 4
        Weights = zeros(dx,ncomp,outClass);
    end
end

% An orthonormal basis for the span of the X loadings, to make the successive
% deflation X0'*Y0 simple - each new basis vector can be removed from Cov
% separately.
V = zeros(dx,ncomp);

Cov = X0'*Y0;
for i = 1:ncomp
    % Find unit length ti=X0*ri and ui=Y0*ci whose covariance, ri'*X0'*Y0*ci, is
    % jointly maximized, subject to ti'*tj=0 for j=1:(i-1).
    [ri,si,ci] = svd(Cov,'econ'); ri = ri(:,1); ci = ci(:,1); si = si(1);
    ti = X0*ri;
    normti = norm(ti); ti = ti ./ normti; % ti'*ti == 1
    Xloadings(:,i) = X0'*ti;
    
    qi = si*ci/normti; % = Y0'*ti
    Yloadings(:,i) = qi;
    
    if nargout > 2
        Xscores(:,i) = ti;
        Yscores(:,i) = Y0*qi; % = Y0*(Y0'*ti), and proportional to Y0*ci
        if nargout > 4
            Weights(:,i) = ri ./ normti; % rescaled to make ri'*X0'*X0*ri == ti'*ti == 1
        end
    end

    % Update the orthonormal basis with modified Gram Schmidt (more stable),
    % repeated twice (ditto).
    vi = Xloadings(:,i);
    for repeat = 1:2
        for j = 1:i-1
            vj = V(:,j);
            vi = vi - (vj'*vi)*vj;
        end
    end
    vi = vi ./ norm(vi);
    V(:,i) = vi;

    % Deflate Cov, i.e. project onto the ortho-complement of the X loadings.
    % First remove projections along the current basis vector, then remove any
    % component along previous basis vectors that's crept in as noise from
    % previous deflations.
    Cov = Cov - vi*(vi'*Cov);
    Vi = V(:,1:i);
    Cov = Cov - Vi*(Vi'*Cov);
end

if nargout > 2
    % By convention, orthogonalize the Y scores w.r.t. the preceding Xscores,
    % i.e. XSCORES'*YSCORES will be lower triangular.  This gives, in effect, only
    % the "new" contribution to the Y scores for each PLS component.  It is also
    % consistent with the PLS-1/PLS-2 algorithms, where the Y scores are computed
    % as linear combinations of a successively-deflated Y0.  Use modified
    % Gram-Schmidt, repeated twice.
    for i = 1:ncomp
        ui = Yscores(:,i);
        for repeat = 1:2
            for j = 1:i-1
                tj = Xscores(:,j);
                ui = ui - (tj'*ui)*tj;
            end
        end
        Yscores(:,i) = ui;
    end
end


%------------------------------------------------------------------------------
%PLSCV Efficient cross-validation for X and Y mean squared error in PLS.
function mse_data = plscv(X,Y,ncomp,cvp,mcreps,ParOptions)

[n,dx] = size(X);

% Return error for as many components as asked for; some columns may be NaN
% if ncomp is too large for CV.
mse = NaN(2,ncomp+1);

% The CV training sets are smaller than the full data; may not be able to fit as
% many PLS components.  Do the best we can.
if isa(cvp,'cvpartition')
    cvpType = 'partition';
    maxncomp = min(min(cvp.TrainSize)-1,dx);
else
    cvpType = 'Kfold';
%    maxncomp = min(min( floor((n*(cvp-1)/cvp)-1), dx));
    maxncomp = min( floor((n*(cvp-1)/cvp)-1), dx);
end
if ncomp > maxncomp
    warning(message('stats:plsregress:MaxComponentsCV', maxncomp));
    ncomp = maxncomp;
end

% Cross-validate sum of squared errors for models with 1:ncomp components,
% simultaneously.  Sum the SSEs over CV sets, and compute the mean squared
% error
CVfun = @(Xtr,Ytr,Xtst,Ytst) sseCV(Xtr,Ytr,Xtst,Ytst,ncomp);
sumsqerr = crossval(CVfun,X,Y,cvpType,cvp,'mcreps',mcreps,'options',ParOptions);
%mse(:,1:ncomp+1) = reshape(sum(sumsqerr,1)/(n*mcreps), [2,ncomp+1]);

data_x = 2*sumsqerr(:,1:2:end)/mcreps;
data_y = 2*sumsqerr(:,2:2:end)/mcreps;

mse(1,1:ncomp+1) = mean(data_x);
mse(2,1:ncomp+1) = mean(data_y);

se = nan(size(mse));

se(1,1:ncomp+1) = std(data_x)/sqrt(size(data_x,1));
se(2,1:ncomp+1) = std(data_y)/sqrt(size(data_y,1));

minMSE1 = min(mse(1,:));
minMSE2 = min(mse(2,:));
minIx1 = find(mse(1,:)==minMSE1,1);
minIx2 = find(mse(2,:)==minMSE2,1);
minplus1 = mse(1,minIx1) + se(1,minIx1);
minplus2 = mse(2,minIx2) + se(2,minIx2);
seIx1= find((mse(1,1:minIx1) <= minplus1),1,'first');
seIx2= find((mse(2,1:minIx2) <= minplus2),1,'first');

mse_data.MSE = mse;
mse_data.SE           = se;
mse_data.IndexMinMSE  = [minIx1,minIx2];
mse_data.Index1SE     = [seIx1,seIx2];

%------------------------------------------------------------------------------
%SSECV Sum of squared errors for cross-validation
function sumsqerr = sseCV(Xtrain,Ytrain,Xtest,Ytest,ncomp)

XmeanTrain = mean(Xtrain);
YmeanTrain = mean(Ytrain);
X0train = bsxfun(@minus, Xtrain, XmeanTrain);
Y0train = bsxfun(@minus, Ytrain, YmeanTrain);

% Get and center the test data
X0test = bsxfun(@minus, Xtest, XmeanTrain);
Y0test = bsxfun(@minus, Ytest, YmeanTrain);

% Fit the full model, models with 1:(ncomp-1) components are nested within
[Xloadings,Yloadings,~,~,Weights] = simpls(X0train,Y0train,ncomp);
XscoresTest = X0test * Weights;

% Return error for as many components as the asked for.
outClass = superiorfloat(Xtrain,Ytrain);
sumsqerr = zeros(2,ncomp+1,outClass); % this will get reshaped to a row by CROSSVAL

% Sum of squared errors for the null model
sumsqerr(1,1) = sum(sum(abs(X0test).^2, 2));
sumsqerr(2,1) = sum(sum(abs(Y0test).^2, 2));

% Compute sum of squared errors for models with 1:ncomp components
for i = 1:ncomp
    X0reconstructed = XscoresTest(:,1:i) * Xloadings(:,1:i)';
    sumsqerr(1,i+1) = sum(sum(abs(X0test - X0reconstructed).^2, 2));

    Y0reconstructed = XscoresTest(:,1:i) * Yloadings(:,1:i)';
    sumsqerr(2,i+1) = sum(sum(abs(Y0test - Y0reconstructed).^2, 2));
end

