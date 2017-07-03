function [FF,MM,WW] = periodic_track_imm(Y,dt,freqlist,nharm,BF,BQ,BL,BH,R,qr,ptrans,poverall)
% PERIODIC_TRACK_IMM - Track periodic signal with IMM
%
% Syntax:
%   [FF,MM,WW] = periodic_track_imm(Y,dt,freqlist,nharm,BF,BQ,BL,BH,R,qr,ptrans)
%
% In:
%   Y  - Signal as Dx1 or 1xD vector
%   dt - Sampling period in seconds (e.g. 0.1)
%   freqlist - List of candidate frequencies in beats-per-min (e.g. 65:75).
%   nharm - Number of harmonics to be estimated (e.g. 3)
%   BF - Feedback matrix for the bias model (e.g. [0 1; 0 0])
%   BQ - Process spectral density for the bias model (e.g. 0.01)
%   BL - Noise multiplier matrix for the bias model (e.g. [0;1])
%   BH - Measurement matrix for the bias model (e.g. [1 0])
%   R  - Measurement variance (e.g. 0.1^2)
%   qr - Resonator's process noise spectral density (e.g. 0.01)
%   ptrans - Transition probability for model change (e.g. 0.001)
%   poverall - Overall transition jump probability (e.g. 0.00)
%
% Out:
%   FF - MMSE estimates of fundamental frequencies as 1xD vector in BPM
%   MM - Filtered means of state components as (2*nharm*dim(bias))xD matrix
%   WW - Smoothed posterior probabilities for models as NxD matrix
%
% Description:
%   Tracks cardiac signal, which is buried in signal Y using the
%   following kind of model:
%
%     Y(k) = xb(k) + xr1(k) + ... + xrn(k) + e(k)
%
%   where e(k) is the measurement noise with variance R and
%
%   - xb(k) is the bias signal, which is modeled the process
%
%        X(k) = BA X(k-1) + w_k,  w_k ~ N(0,BQ)
%       xb(k) = BH X(k),
%
%   - xrj(k) are resonators (nharm of them), each having frequency j*f0,
%     where f0 is the fundamental frequency. The model can be written as
%
%        dRj/dt = [0 2*pi*j*f0; -2*pi*j*f0 0] Rj + wj(t),
%        xrj(k) = [1 0] Rj
%
%     where each wj has spectral density qr.
%
%   - The possible frequencies f0(j) given in freqlist form a HMM with
%     transition probabilities P(j+1|j) = P(j-1|j) = ptrans.
%
%   - The estimation is based on the Interacting Multiple Models
%     (IMM) algorithm, which used Kalman filters for tracking
%     each of the modes and mixing between them.
%

% Copyright (C) 2010 Simo Särkkä and Arno Solin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

%% 

    % Report what we are doing
    fprintf('%-36s:%+35s\n','IMM Tracking of frequencies','Preparing..')

    % Change frequencies to rad/s
    ww = 2*pi*freqlist/60; % To rads/s

    %
    % Compute transition probability matrix
    %
    NM = length(ww);
    pp = 1-ptrans;
    PIJ = poverall*ones(NM) + (pp-NM*poverall)*eye(NM) + ...
                       (1-pp)/2*diag(ones(NM-1,1),1) + ...
                       (1-pp)/2*diag(ones(NM-1,1),-1);
    PIJ(1,1)   = PIJ(1,1) + (1-pp)/2;
    PIJ(NM,NM) = PIJ(NM,NM) + (1-pp)/2;
    
    % Check if the PIJ transition matrix is consistent
    if (sum(PIJ)-NM) > 1e-10 | max(abs(PIJ)) > 1, 
        warning(['The transition matrix Pij might be inconsistent. ' ...
                 'Check your transition probabilities.'])
    end
    
    %
    % Compute the IMM parameters (ML,PL,AL,QL,HL,RL)
    %
    N = nharm;

    % Prior means and covariances
    m0 = zeros(2*N+size(BF,1),1);
    P0 = 100*eye(size(m0,1));

    RL = zeros(1,1,length(ww));
    HL = zeros(1,size(m0,1),length(ww));
    AL = zeros(size(m0,1),size(m0,1),length(ww));
    QL = zeros(size(m0,1),size(m0,1),length(ww));

    ML = repmat(m0,[1 length(ww)]);
    PL = repmat(P0,[1 1 length(ww)]);

    for i=1:length(ww)
        H  = zeros(1,size(m0,1));
        F  = zeros(size(m0,1));
        L  = zeros(size(m0,1),N+size(BQ,1));
        Qc = zeros(N+size(BQ,1));
        
        %
        % Form the resonators
        %
        for j=1:N
            i1 = 1+2*(j-1);
            i2 = 2+2*(j-1);
            F(i1:i2,i1:i2) = [0 j*ww(i); -j*ww(i) 0];
            L(i2,j) = 1;
            Qc(j,j) = qr;
            H(1,i1) = 1;
        end
        
        %
        % Form the bias model
        %
        i1 = 1+2*N;
        i2 = 2+2*N;
        j1 = N+2-size(BQ,1);
        j2 = N+1;
        F(i1:i2,i1:i2)  = BF;
        L(i1:i2,j1:j2)  = BL;
        H(1,i1:i2)      = BH;    
        Qc(j1:j2,j1:j2) = BQ;
 
        %
        % Discretize
        %
        [AL(:,:,i),QL(:,:,i)] = lti_disc(F,L,Qc,dt);
        HL(:,:,i) = H;    
        RL(:,:,i) = R;
        
    end

    % Show waitbar
    %handle = waitbar(0,'Please wait...');  
    loopstart = tic;
    
    %
    % Estimate with IMM
    %
    WW = zeros(NM,length(Y));
    W  = ones(1,NM)/NM;

    FF = zeros(1,length(Y));
    MM = zeros(size(m0,1),length(Y));
    for k=1:length(Y)
        %
        % IMM prediction
        %
        [ML,PL,W] = imm_predict0(ML,PL,W,PIJ,AL,QL);
        
        %
        % IMM update
        %
        [ML,PL,W,K,IM,IS,LH] = imm_update0(ML,PL,W,Y(k),HL,RL);
        WW(:,k) = W;
        FF(k) = sum(W .* ww) / sum(W);
        MM(:,k) = imm_estimate0(ML,PL,W);
        
        % Update waitbar and show time remaining
        if rem(k,40)==0 
          secondsleft = min(toc(loopstart),Inf)/k*(length(Y)-k);
          %waitbar(k/length(Y),handle,sprintf('Estimating with IMM\nTime left: %.0f min %.0f s.', ...
          %  floor(secondsleft/60),rem(secondsleft,60)))
          fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
              sprintf('%.0f min %.0f s', ...
                floor(secondsleft/60),rem(secondsleft,60)))
        end    
    end
    
    % Dump
    % save FF_1.mat FF

    %
    % Smooth the mode probilities
    %
    thr = eps;
    SWW = WW;
    SFF = FF;
    for k=length(Y)-1:-1:1
        PW = PIJ*WW(:,k);   % p(xk+1|Y1:k)
        PW(PW < thr) = thr;
        SW = SWW(:,k+1)./PW; % p(xk+1|Y1:T) / p(xk+1|Y1:k)
        SW = PIJ'*SW;      % int p(xk+1|xk) p(xk+1|Y1:T) / p(xk+1|Y1:k) dxk+1
        SW = SW .* WW(:,k); % * p(xk|y1:k)
        SWW(:,k) = SW / sum(SW);
        SFF(k) = sum(SW' .* ww) / sum(SW);
    end

    % Dump
    % save FF_2.mat FF SWW SFF PW thr PIJ SW ww WW

    % Get rid of the waitbar
    %close(handle)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
            'Done.')
                 
    % Check if we have NaN values
    if sum(isnan(SFF)) == 0
      % The smoother estimate had no problems
      FF = 60*SFF/2/pi;
    elseif sum(isnan(FF)) == 0
      % Use the filter estimate instead (it is ok)
      FF = 60*FF/2/pi;  
      warning('NaN in IMM smoother, using filter result instead.')  
    else
        
      % There were problems in the reference, identify them
      ind = find(isnan(FF),1,'first');

      % Show figure of problem
      figure;
        T = dt*(0:length(FF)-1);
        plot(T,Y,'-k',T(ind:end),Y(ind:end),'-r','LineWidth',2)
        xlabel('Time [s]'); legend('Ok signal','Problematic part starts')
        title('Trying to point out the problems in your reference signal')
        
      % Throw error
      error(['Problem in IMM estimation around time point: %.1f s. ' ...
          'There are probably problems with your reference signal'],dt*ind)

    end