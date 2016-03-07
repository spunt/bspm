function [C,h,Ph,F,Fa,Fc,k] = spm_rwls_reml(YY,X,Q,N,hE,hC,A,K)
% ReML estimation of covariance components from y*y' 
% Exploits special structure of first and second derivative for the
% weighted regression model: The first m components are assumed to be 
% a variance scaling parameters for each image
% 
% FORMAT [C,h,Ph,F,Fa,Fc,k] = spm_rwls_reml(YY,X,Q,N,[hE,hC,A,K]);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
% N   - number of samples
%
% hE  - hyperprior expectation in log-space [default = 0]
% hC  - hyperprior covariance  in log-space [default = 256]
% A   - proportional hyperpriors [default = 0, yes]
% K   - number of iterations [default = 32]
%
% C   - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) ReML hyperparameters h
% Ph  - (q x q) conditional precision of log(h)
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
%
% Fa  - accuracy
% Fc  - complexity (F = Fa - Fc)
%
% k   - number of iterations required
%
% Performs a Fisher-Scoring ascent on F to find MAP variance parameter
% estimates.  
% It assumes that the first m components are indivudual variances of 
% the observations (as required by the rWLS toolbox)
% Other components can be estimated
% It estimates parameter in logspace (ensuring positive definite outcomes) 
% and uses weakly informative log-normal hyperpriors.
% 
% Most of the code is taken from the original spm_reml_sc
%__________________________________________________________________________
%
% SPM ReML routines:
%
%      spm_reml:    no positivity constraints on covariance parameters
%      spm_reml_sc: positivity constraints on covariance parameters
%      spm_sp_reml: for sparse patterns (c.f., ARD)
%
%__________________________________________________________________________
% 
% 
% 
% 
% 
 
% Joern Diedrichsen, Karl Friston 
% $Id: spm_rwls_reml.m joern  $

% assume proportional hyperpriors not specified
%--------------------------------------------------------------------------
try, A; catch, A  = 0;  end
 
% assume a single sample if not specified
%--------------------------------------------------------------------------
try, N; catch, N  = 1;  end
 
% default number of iterations
%--------------------------------------------------------------------------
try, K; catch, K  = 64; end
 
% initialise h
%--------------------------------------------------------------------------
n     = length(Q{1});
m     = length(Q);
h     = zeros(m,1);
dh    = zeros(m,1);
dFdh  = zeros(m,1);
dFdhh = zeros(m,m);
 
% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(n,0);
    R = speye(n,n);
else
    X = orth(full(X));
    R = speye(n,n) - X*X';
end
 
 
% initialise and specify hyperpriors
%==========================================================================

% scale Q and YY
%--------------------------------------------------------------------------
if A
    sY = trace(R*YY)/N/n;
    YY = YY/sY;
    for i = 1:m
        sh(i,1) = trace(R*Q{i})/n;
        Q{i}    = Q{i}/sh(i);
    end
    keyboard; 
else
    sY = 1;
    sh = 1;
end

% hyperpriors
%--------------------------------------------------------------------------
try, hE = hE(:);                               catch, hE = 0;   end
try, hP = inv(hC + speye(length(hC))/exp(16)); catch, hP = 1/256; end
 
% check sise
%--------------------------------------------------------------------------
if length(hE) < m, hE = hE(1)*ones(m,1);   end
if length(hP) < m, hP = hP(1)*speye(m,m);  end

% intialise h: so that sum(exp(h)) = 1
%--------------------------------------------------------------------------
if any(diag(hP) >  exp(16))
    h = hE;
end
 
% ReML (EM/VB)
%--------------------------------------------------------------------------
dF    = Inf;
ds    = 1:m;
for k = 1:K
 
    % compute current estimate of covariance
    %----------------------------------------------------------------------
    C     = sparse(n,n);
    for i = 1:m
        C = C + Q{i}*exp(h(i));
    end
    iC    = inv(C + speye(n,n)/exp(32));
 
    % E-step: conditional covariance cov(B|y) {Cq}
    %======================================================================
    iCX    = iC*X;
    if ~isempty(X)
        Cq = inv(X'*iCX);
    else
        Cq = sparse(0);
    end
 
    % M-step: ReML estimate of hyperparameters
    %======================================================================
 
    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    P     = iC - iCX*Cq*iCX';
    dFdh=-diag(P)/2+1/2*sum((P*YY).*P',2);
 
    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    dFdhh=-0.5*(P.*P');
 
    % Check if other components are present 
    %----------------------------------------------------------------------
    dFdhh=-0.5*(P.*P');
    
    
    % modulate
    %----------------------------------------------------------------------
    dFdh  = dFdh.*exp(h);
    dFdhh = dFdhh.*(exp(h)*exp(h)');
 
    % add hyperpriors
    %----------------------------------------------------------------------
    e     = h     - hE;
    dFdh  = dFdh  - hP*e;
    dFdhh = dFdhh - hP;
 
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    dh    = spm_dx(dFdhh,dFdh)*exp(-k/(K/2));
    h = h + dh;
    
    % Convergence (1% change in log-evidence)
    %======================================================================
    % bar(h);drawnow
    
    % convergence
    %----------------------------------------------------------------------
    dF    = dFdh'*dh;
    fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(dF));
    if dF < 1e-2
        break
    end;
end

% log evidence = ln p(y|X,Q) = ReML objective = F = trace(R'*iC*R*YY)/2 ...
%--------------------------------------------------------------------------
Ph    = -dFdhh;
if nargout > 3
 
    % tr(hP*inv(Ph)) - nh + tr(pP*inv(Pp)) - np (pP = 0)
    %----------------------------------------------------------------------
    Ft = trace(hP*inv(Ph)) - length(Ph) - length(Cq);
 
    % complexity - KL(Ph,hP)
    %----------------------------------------------------------------------
    Fc = Ft/2 + e'*hP*e/2 + spm_logdet(Ph*inv(hP))/2 - N*spm_logdet(Cq)/2;
 
    % Accuracy - ln p(Y|h)
    %----------------------------------------------------------------------
    Fa = Ft/2 - trace(C*P*YY*P)/2 - N*n*log(2*pi)/2  - N*spm_logdet(C)/2;
 
    % Free-energy
    %----------------------------------------------------------------------
    F  = Fa - Fc - N*n*log(sY)/2;
 
end


% return exp(h) hyperpriors and rescale
%--------------------------------------------------------------------------
h  = sY*exp(h)./sh;
C  = sY*C;
