function [S,CS,QCS,SMM,SPP,FMM,FPP] = periodic_separation_kfs(Y,dt,FF,nharm,BF,BQ,BL,BH,R,qr)
%CARDIAC_SEPARATION_KFS - Separate periodic signal from other signals
%
% Syntax:
%   [S,CS,QCS,SMM,SPP,FMM,FPP] = periodic_separation_kfs(Y,dt,FF,nharm,BF,BQ,BL,BH,R,qr)
%
% In:
%   Y  - Signal as 1xD vector
%   dt - Sampling period in seconds (e.g. 0.1)
%   FF - Fundamental frequencies as 1xD vector in BPM (from the IMM)
%   nharm - Number of harmonics to be estimated (e.g. 3)
%   BF - Feedback matrix for the bias model (e.g. [0 1; 0 0])
%   BQ - Process spectral density for the bias model (e.g. 0.01)
%   BL - Noise multiplier matrix for the bias model (e.g. [0;1])
%   BH - Measurement matrix for the bias model (e.g. [1 0])
%   R  - Measurement variance (e.g. 0.1^2)
%   qr - Resonator's process noise spectral density (e.g. 0.1)
%
% Out:
%   S   - Cleaned signal as 1xD vector
%   CS  - Fundamental signal and its harmonics as (nharm)xD vector
%   QCS - Quadrature periodic signals as (nharm)xD vector
%   SMM - The smoother means NxD
%   SPP - The smoother covariances NxNxD
%   FMM - The filter means NxD
%   FPP - The filter covariances NxNxD
%
% Description:
%   Separate given signal into true signal and periodic
%   signal by using the known frequency of the periodic signal.
%   See PERIODIC_TRACK_IMM for the model used.

% Copyright (C) 2011 Simo Särkkä and Arno Solin
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

    % Number of harmonics (including fundamental)
    N = nharm;
    
    % Prior means and covariances
    m0 = zeros(2*N+size(BF,1),1);  % XXX: Could be function parameter
    P0 = 10*eye(size(m0,1));       % XXX: Could be function parameter

    % Mean to measurement at step 1
    m0(1) = Y(1);
    
    %
    % Form the constant matrices
    %
    F  = zeros(size(m0,1));
    H  = zeros(1,size(m0,1));
    L  = zeros(size(m0,1),N+size(BQ,1));
    Qc = zeros(N+size(BQ,1));
    
    for j=1:N
        i1 = 1+2*(j-1);
        i2 = 2+2*(j-1);
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
  
    
%% Run filter
    
    m = m0;
    P = P0;

    FMM = zeros(size(m,1),length(Y));
    FPP = zeros(size(P,1),size(P,2),length(Y));
    
    for k=1:length(Y)
        
        %
        % Update the resonators
        %
        w  = 2*pi*FF(k)/60;
        for j=1:N
            i1 = 1+2*(j-1);
            i2 = 2+2*(j-1);
            F(i1:i2,i1:i2) = [0 j*w; -j*w 0];
        end
        
        %
        % Discretize
        %
        [A,Q] = lti_disc(F,L,Qc,dt);

        %
        % Estimate with KF
        %
        [m,P] = kf_predict(m,P,A,Q);
        [m,P] = kf_update(m,P,Y(k),H,R);
    
        FMM(:,k) = m;
        FPP(:,:,k) = P;
    end
    
    
%% Run smoother

    % Allocate and set initial
    SMM = zeros(size(FMM));
    SPP = zeros(size(FPP));
    SMM(:,end)   = m;
    SPP(:,:,end) = P;
            
    for k=length(Y)-1:-1:1
        
        %
        % Update each resonator
        %
        w  = 2*pi*FF(k+1)/60;
        for j=1:N
            i1 = 1+2*(j-1);
            i2 = 2+2*(j-1);
            F(i1:i2,i1:i2) = [0 j*w; -j*w 0];
        end
        
        %
        % Discretize
        %
        [A,Q] = lti_disc(F,L,Qc,dt);
        
        % Smoother
        m_pred = A * FMM(:,k);
        P_pred = A * FPP(:,:,k) * A' + Q;
        D  = FPP(:,:,k) * A' / P_pred;
        m  = FMM(:,k) + D * (m - m_pred);
        P  = FPP(:,:,k) + D * (P - P_pred) * D';
        
        % Store results
        SMM(:,k)   = m;
        SPP(:,:,k) = P;
        
    end
    
    %
    % Extract the signals
    %
    S   = SMM(end-size(BF,1)+1,:);
    CS  = SMM(1:2:(2*nharm),:);
    QCS = SMM(2:2:(2*nharm+1),:);
    