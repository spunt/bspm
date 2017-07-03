function [SMM,SPP,FMM,FPP] = periodic_separation_kfs3D(Y,dt,FF,N,BF,BQ,BL,BH,R,qr)
% PERIODIC_SEPARATION_KFS3D - Separate periodic signals from other signals in 3D
%
% Syntax:
%   [SMM,SPP,FMM,FPP] = periodic_separation_kfs3D(Y,dt,FF,nharm,BF,BQ,BL,BH,R,qr)
%
% In:
%   Y  - Signal as D1xD2xD3xT matrix
%   dt - Sampling period in seconds (e.g. 0.1)
%   FF - Fundamental frequencies as 1xD vector in bpm (from the IMM)
%   N  - Number of harmonics to be estimated (e.g. [3 3])
%   BF - Feedback matrix for the bias model (e.g. [0 1; 0 0])
%   BQ - Process spectral density for the bias model (e.g. 0.01)
%   BL - Noise multiplier matrix for the bias model (e.g. [0;1])
%   BH - Measurement matrix for the bias model (e.g. [1 0])
%   R  - Measurement variance (e.g. 0.1^2)
%   qr - Resonator's process noise spectral density (e.g. 0.1)
%
% Out:
%   SMM - The smoother means NxT
%   SPP - The smoother covariances NxNxT
%   FMM - The filter means NxD
%   FPP - The filter covariances NxNxT
%
% Description:
%   Separate given signal into bias and periodic signals
%   by using the known frequency of the periodics.
%   See CARDIAC_TRACK_IMM for the model used.
%     This version uses data of dimension D1xD2xD3 with T time samples.
%   The harmonic components are separately separated for each voxel.
%   The linear Kalman filter uses the same covariance and gain term
%   for each voxel.

% Copyright (C) 2011-2012 Simo Sarkka, Arno Solin
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

%% Set up

    % Report what we are doing
    fprintf('%-36s:%+35s\n','Kalman filtering','Preparing..')

    % Get data dimensions
    ydim  = size(Y,1);
    xdim  = size(Y,2);
    zdim  = size(Y,3);
    steps = size(Y,4);

    % Number of harmonics and components for a single point
    nsize = 2*sum(N)+size(BF,1);
    
    % Initialize means and covariances
    m0 = zeros(nsize,xdim*ydim*zdim);
    y = Y(:,:,:,1);
    if (size(BF,1) > 0), m0(2*sum(N)+1:nsize:end) = y(:); end
    P0 = 100*eye(nsize);

    % Form the constant matrices
    F  = zeros(nsize);
    H  = zeros(1,nsize);
    L  = zeros(nsize,sum(N)+size(BQ,1));
    Qc = zeros(sum(N)+size(BQ,1));

    % All periodics
    for k=1:length(N)
      pos = sum(N(1:k-1));
      for j=1:N(k)
        i1 = 1+2*(j-1)+2*pos;
        i2 = 2+2*(j-1)+2*pos;
        L(i2,j+pos) = 1;
        Qc(j+pos,j+pos) = qr(k)/j;
        H(1,i1) = 1;
      end
    end

    % Form the bias model
    if (size(BF,1)>0)
      i1 = 1+2*sum(N);
      i2 = size(BF,1)+2*sum(N);
      j1 = sum(N)+1;
      j2 = sum(N)+size(BQ,1);
      F(i1:i2,i1:i2)  = BF;
      L(i1:i2,j1:j2)  = BL;
      H(1,i1:i2)      = BH;
      Qc(j1:j2,j1:j2) = BQ; 
    end
    
%% Track with KF

    % Initial values
    m = m0;
    P = P0;

    % Allocate space and set initial values
    FMM = zeros(size(m,1),size(m,2),steps);
    FPP = zeros(size(P,1),size(P,2),steps);
    
    % Show waitbar
    %handle = waitbar(0,'Please wait...');  
    loopstart = tic;
    
    for k=1:steps
        
        %
        % Update each resonator
        %
        for l=1:length(N)
          pos = sum(N(1:l-1));
          w  = 2*pi*FF(k,l)/60;
          for j=1:N(l)
            i1 = 1+2*(j-1)+2*pos;
            i2 = 2+2*(j-1)+2*pos;
            F(i1:i2,i1:i2) = [0 j*w; -j*w 0];
          end
            
        end
        
        %
        % Discretize
        %
        [A,Q] = lti_disc(F,L,Qc,dt);
        
        %
        % Estimate with KF
        %
        %[m,P] = kf_predict(m,P,A,Q);
        %[m,P] = kf_update(m,P,Y(k),H,R);
    
        % Prediction
        m = A * m;
        P = A * P * A.' + Q;
        
        % Update
        y = Y(:,:,:,k);
        S = H*P*H.' + R;
        K = P*H.' / S;
        m = m + K * (y(:).'-H*m);
        P = P - K*S*K.';
        
        % Store results
        FMM(:,:,k) = m;
        FPP(:,:,k) = P;
        
        % Update waitbar and show time remaining
        if rem(k,20)==0 
          secondsleft = min(toc(loopstart),Inf)/k*(steps-k);
          %waitbar(k/steps,handle,sprintf('Filtering\nTime left: %.0f min %.0f s.', ...
          %  floor(secondsleft/60),rem(secondsleft,60)))
        
          fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
              sprintf('%.0f min %.0f s', ...
                floor(secondsleft/60),rem(secondsleft,60)))
        
        end
    end

    % Get rid of the waitbar
    %close(handle)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
            'Done.')
        
    
%% Run smoother

    % Report what we are doing
    fprintf('%-36s:%+35s\n','RTS smoother','Preparing..')

    % Allocate and set initial
    SMM = zeros(size(m,1),size(m,2),steps);
    SPP = zeros(size(P,1),size(P,2),steps);
    SMM(:,:,end) = m;
    SPP(:,:,end) = P;
        
    % Show waitbar
    %handle = waitbar(0,'Please wait...');  
    loopstart = tic;
    
    for k=steps-1:-1:1
        
        %
        % Update each resonator
        %
        for l=1:length(N)
          pos = sum(N(1:l-1));
          w  = 2*pi*FF(k+1,l)/60;
          for j=1:N(l)
            i1 = 1+2*(j-1)+2*pos;
            i2 = 2+2*(j-1)+2*pos;
            F(i1:i2,i1:i2) = [0 j*w; -j*w 0];
          end
            
        end
        
        %
        % Discretize
        %
        [A,Q] = lti_disc(F,L,Qc,dt);
        
        % Smoother
        m_pred = A * FMM(:,:,k);
        P_pred = A * FPP(:,:,k) * A.' + Q;
        D  = FPP(:,:,k) * A.' / P_pred;
        m  = FMM(:,:,k) + D * (m - m_pred);
        P  = FPP(:,:,k) + D * (P - P_pred) * D.';
        
        % Store results
        SMM(:,:,k) = m;
        SPP(:,:,k) = P;
        
        % Update waitbar and show time remaining
        if rem(k,20)==0 
          secondsleft = min(toc(loopstart),Inf)/(steps-k)*k;
          %waitbar(k/steps,handle,sprintf('Smoothing\nTime left: %.0f min %.0f s.', ...
          %floor(secondsleft/60),rem(secondsleft,60)))
          
          fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
              sprintf('%.0f min %.0f s', ...
                floor(secondsleft/60),rem(secondsleft,60)))
          
        end
    end
    
    % Get rid of the waitbar
    %close(handle)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
            'Done.')
    
        
%% Extract the signals
    
    FMM = reshape(FMM,nsize*ydim*xdim*zdim,[]);
    SMM = reshape(SMM,nsize*ydim*xdim*zdim,[]);
