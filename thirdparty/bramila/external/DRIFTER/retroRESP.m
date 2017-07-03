function [phi] = retroRESP(R,dt)
% RETRORR - Finds the respiratory phase
%
% Syntax:
%   [phi] = retroRESP(data,dt)
%
% In:
%   R  - external respiratory signal as a 1xN vector
%   dt - optional (default 0.001)
%
% Out:
%   phi - the respiratory phase as in the RETROICOR article
%
% Description:
%   See the articles for details: 
%     Glover et al. (2000) Image-Based Method for Retrospective Correction 
%   of Physiological Motion Effects in fMRI: RETROICOR. Magnetic Resonance
%   in Medicine, 44, 162-167.
%     Savitzky A, Golay M. (1964) Smoothing and differentiation of data by 
%   simplified least squares procedures. Anal Chem 36:1627–1639.
%

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

%% Simulate input

  if (nargin==0)
    dt = 1/40;
    T = 0:dt:124;
    R = sin(T)+0.1*rand(size(T));
    plot(T,R)
  end

%% Set up
  
  % Time for plots
  T = 0:dt:dt*(length(R)-1);
    
  % Normalize to (0,Rmax)
  R    = R(:)-min(R);
  Rmax = max(R);

  % Create bins (N = number of bins)
  N = 100;
  b = linspace(0,Rmax,N+1);

  % Calculate derivative by first smoothing the data
  % Use Savitzky-Golay Filtering as in the Glover article
  sR = sgolayfilt(R,2,39);
  dR = diff(sR);

  
%% Calculate histogram H(b)

  H = zeros(N,1);

  for k=1:N
      
     % Count occurances
     H(k) = sum(R >= b(k) & R < b(k+1));
     
  end
  
  
%% Calculate the phase  
  
  phi = zeros(size(R));

  for k=1:length(R)-1
      
      phi(k) = pi*sum(H(1:floor(N*R(k)/Rmax)))/sum(H) * sign(dR(k));
      
  end

%% Calculate frequency estimate
%{
  ff=diff(phi);
  ind = ff > 1;
  Ti = T(ind);
  ff=Ti(2:end)-Ti(1:end-1);
  ff=interp1(Ti(2:end),1./ff,T,'linear','extrap');
%}
  
%% Plot original, histogram, phase

  if (nargout == 0)

    subplot(311)
      plot(T,R,'-b',T,sR,'-g')
      title('\bf Data'); xlabel('T [s]'); ylabel('Amplitude')
      legend('Original','Savitzky-Golay Filtered','location','best')
    subplot(312)
      stairs(b(1:end-1),H)
      title('\bf Histogram'); xlabel('b'); ylabel('H(b)')
    subplot(313)
      plot(T,phi)
      title('\bf Phase'); xlabel('T [s]'); ylabel('\phi(t)')
  end
  
% Plot histograms
  if (nargout == 0)
    figure(3)
    subplot(121); hold on
      specgram(sin(phi), 512, 1/(1/40));
      plot(T,ff)
      ylim([0 2])
      title('\bf From Estimated Phase')
    subplot(122); hold on
      specgram(R, 512, 1/(1/40));
      plot(T,ff)
      ylim([0 2])
      title('\bf Original Data')
  end
  
  
  
  
  
  
  
  