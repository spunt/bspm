function [clean,card,resp] = retroicor(data_epi,data_ecg,data_resp,dt,mc,mr)
% RETROICOR - Clean EPI signal with RETROICOR
%
% Syntax:
%   [clean,card,resp] = RETROICOR(data_epi,data_ecg,data_resp,dt,mc,mr)
%
% In:
%   data_epi  - XxYxZxN matrix with EPI data
%   data_ecg  - 1xN vector of external cardiac signal
%   data_resp - 1xN vector of external respiratory signal
%   dt - 3x1 vector of sampling freq. e.g. (dt=1/10,dtC=1/1000,dtR=1/1000)
%   mc - Order of cardiac coefficients (e.g. mc=2)
%   mr - Order of respiratory coefficients (e.g. mr=2)
%
% Out:
%   clean - XxYxZxN matrix of cleaned EPI data (y-y_delta)
%   card  - XxYxZxN matrix of cardiac contribution
%   resp  - XxYxZxN matrix of respiratory contribution
%
% Description:
%     This method is based on the article:
%   Glover et al. (2000) Image-Based Method for Retrospective Correction 
%   of Physiological Motion Effects in fMRI: RETROICOR. Magnetic Resonance
%   in Medicine, 44, 162-167.
%
% See also:
%   retroRR, retroRESP

% Copyright:
%   Arno Solin, 2011
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

%% Setup
  
  % Use SPM style progress reporting (no graphical output)
  fprintf('%-36s:%+35s\n','RETROICOR','Preparing..')

  % Sampling intervals of the signals
  dtC = dt(2);
  dtR = dt(3);
  dt  = dt(1);
  
  % Downsample respiratory signal 
  % (to make the method match Glover et al.)
  ecg = data_ecg;
  R   = interp1(dtR*(0:length(data_resp)-1),data_resp, ...
                0:1/40:dtR*(length(data_resp)-1));
  
  % There is no need to remove the mean here          
  data_mean = mean(data_epi,4);
  data_epi = bsxfun(@minus,data_epi,data_mean); % Remove mean
  
  % Initialize time vector and dimensions
  xdim = size(data_epi,1); 
  ydim = size(data_epi,2);
  zdim = size(data_epi,3);
  T = 0:dt:dt*(size(data_epi,4)-1);
  

%% Find the cardiac and respiratory phases

  % Cardiac phase
  [phiC] = retroRR(ecg,dtC);
  
  % Respiratory phase
  [phiR] = retroRESP(R,1/40);
  
  % Interpolate to match fMRI sampling TR
  phiC = interp1(0:dtC:dtC*(length(phiC)-1),phiC,T,'linear','extrap');
  phiR = interp1(0:1/40:1/40*(length(phiR)-1),phiR,T,'linear','extrap');
  
  
  
  
%% Do the retroicor estimation per voxel

  % Allocate matrices for the Fourier coefficients
  ac = zeros(mc,xdim,ydim,zdim);
  bc = ac;
  ar = zeros(mr,xdim,ydim,zdim);
  br = ar;

  % Show waitbar
  %handle = waitbar(0,'Please wait...');  
  loopstart = tic;
  
  % Loop through all pixels
  for i=1:xdim
    for j=1:ydim
      for k=1:zdim
      
        % Mean
        ybar = mean(data_epi(i,j,k,:));
      
        % Cardiac
        for m=1:mc
          ac(m,i,j,k) = sum(squeeze(data_epi(i,j,k,:)-ybar)'.*cos(m*phiC))/ ...
            sum(cos(m*phiC).^2);
          bc(m,i,j,k) = sum(squeeze(data_epi(i,j,k,:)-ybar)'.*sin(m*phiC))/ ...
            sum(sin(m*phiC).^2);
        end
      
        % Respiratory
        for m=1:mr
          ar(m,i,j,k) = sum(squeeze(data_epi(i,j,k,:)-ybar)'.*cos(m*phiR))/ ...
            sum(cos(m*phiR).^2);
          br(m,i,j,k) = sum(squeeze(data_epi(i,j,k,:)-ybar)'.*sin(m*phiR))/ ...
            sum(sin(m*phiR).^2);
        end
      
        % Show time remaining
          l = (i-1)*ydim*zdim + (j-1)*zdim + k*zdim;
          if rem(l,100)==0 
            secondsleft = min(toc(loopstart),Inf)/l*(xdim*ydim*zdim-l);
            %waitbar(l/xdim/ydim/zdim,handle,sprintf('RETROICOR\nTime left: %.0f min %.0f s.', ...
            %  floor(secondsleft/60),rem(secondsleft,60)))
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
              sprintf('%.0f min %.0f s', ...
                floor(secondsleft/60),rem(secondsleft,60)))
          end
      end  
    end
  end

  % Get rid of the waitbar
  %close(handle)
  fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
          'Done.') 

%% Visualize Fourier coefficients (for debugging only)

  if (nargout == 0)

    % Cardiac
    figure(1); clf
    for m=1:mc
      subplot(mc,2,2*m-1)
        imagesc(squeeze(ac(m,:,:))); axis square; colorbar
        title(sprintf('\\bf Cardiac Fourier Sine Coefficients m=%i',m)) 
      subplot(mc,2,2*m)
        imagesc(squeeze(bc(m,:,:))); axis square; colorbar
        title(sprintf('\\bf Cardiac Fourier Cosine Coefficients m=%i',m))      
    end

    % Respiratory
    figure(2); clf
    for m=1:mr
      subplot(mc,2,2*m-1)
        imagesc(squeeze(ar(m,:,:))); axis square; colorbar
        title(sprintf('\\bf Respiratory Fourier Sine Coefficients m=%i',m))      
      subplot(mc,2,2*m)
        imagesc(squeeze(br(m,:,:))); axis square; colorbar
        title(sprintf('\\bf Respiratory Fourier Cosine Coefficients m=%i',m))      
    end

  end
  
%% Visualize separate components (for debugging only)

  if (nargout == 0)

    Zpos = 1;
      
    % The standard deviation brain image 
    amp = zeros(xdim,ydim);
    for i=1:xdim, 
      for j=1:ydim
        amp(i,j) = std(data_epi(i,j,Zpos,200:end));
      end 
    end
  
    % Initial point
    Xpos = 32; Ypos = 32;
  
    % Wait for user interaction
    but = 1;
    while (but == 1)

      % Calculate separate signals (Cardiac)
      card = zeros(size(T));
      for m=1:mc
        card = card + ac(m,Xpos,Ypos,Zpos)*sin(m*phiC) + ...
                      bc(m,Xpos,Ypos,Zpos)*cos(m*phiC); 
      end
      
      % Calculate separate signals (Cardiac)
      resp = zeros(size(T));
      for m=1:mr
        resp = resp + ar(m,Xpos,Ypos,Zpos)*sin(m*phiR) + ...
                      br(m,Xpos,Ypos,Zpos)*cos(m*phiR); 
      end
      
      
      figure(3);clf
    
      % Brain image
      subplot(3,2,1)
        hold on; 
        imagesc(amp)
        plot(Xpos,Ypos,'ow','LineWidth',2)
        hold off
        axis square; axis tight; box on
      
      % Original signal and cleaned
      subplot(3,2,2)
        hold on
        plot(T,squeeze(data_epi(Xpos,Ypos,Zpos,:)),'-k')
        plot(T,squeeze(data_epi(Xpos,Ypos,Zpos,:))'-card-resp,'-b','LineWidth',1)
        hold off
        legend('Original','Cleaned','Location','best')

      % Cardiac
      subplot(3,2,4)
        plot(T,card,'-r')
        title('\bf Cardiac')

      % Cardiac
      subplot(3,2,6)
        plot(T,resp,'-b')
        title('\bf Respiratory')

      % Mouse input
      [xi,yi,but] = ginput(1);
      Xpos = round(xi);
      Ypos = round(yi);
      
    end
  end

  
%% Output signals
  
  clean = zeros(size(data_epi));
  card  = clean;
  resp  = clean;
  

  % Loop through all pixels
  for i=1:xdim
    for j=1:ydim
      for k=1:zdim
     
        % Cardiac
        for m=1:mc
          card(i,j,k,:) = card(i,j,k,:) + shiftdim(ac(m,i,j,k)*cos(m*phiC) + ...
                                             bc(m,i,j,k)*sin(m*phiC),-2);
        end
      
        % Respiratory
        for m=1:mr
          resp(i,j,k,:) = resp(i,j,k,:) + shiftdim(ar(m,i,j,k)*cos(m*phiR) + ...
                                             br(m,i,j,k)*sin(m*phiR),-2); 
        end

      end
    end
  end

  % Output clean
  clean = data_epi - card - resp;
  
  % Add mean back
  clean = bsxfun(@plus,clean,data_mean);


  
  
  
