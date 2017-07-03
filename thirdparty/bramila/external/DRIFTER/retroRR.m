function [phi] = retroRR(ecg,dt)
% RETRORR - Finds the cardiac phase with the help of the R-R time
%
% Syntax:
%   [phi] = retroRR(ecg,dt)
%
% In:
%   ecg - external cardiac signal as a 1xN vector
%   dt  - optional (default 0.001)
%
% Out:
%   phi - the cardiac phase as in the retroicor article
%
% Description:
%   See the article for details: 
%     Glover et al. (2000) Image-Based Method for Retrospective Correction
%   of Physiological Motion Effects in fMRI: RETROICOR. Magnetic Resonance
%   in Medicine, 44, 162-167.
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

%% Calculate heart rate

  % Set dt
  if nargin < 2, dt = 0.001; end

  % Data in ecg
  T  = 0:dt:dt*(length(ecg)-1);
  
  % Smooth data with moving average
  %secg = tsmovavg(ecg, 's', 100,1);
  secg = sgolayfilt(ecg,3,99); 
  
  % Find peak indices by derivatives
  %ind1 = find(diff(diff(secg)>0)<0)-50;
  ind1 = find(diff(diff(secg)>0)<0);

  % We only want the R-R peaks (remove smaller)
  ind2 = find(ecg > 1.2*std(ecg));

  % The peaks are the intersection of these sets
  ind = intersect(ind1,ind2);
  
  % Show results
  if (nargout == 0), figure(1); plot(T,ecg,'-',T(ind),ecg(ind),'*-g'); end
  
  
%% Visualize R-R times

  if (nargout == 0)
    figure(2)
    T1 = T(ind);
    dT = T1(2:end) - T1(1:end-1);
    plot(dT,'.-')
  end
  
  % Calculate frequaency estimate
  T1 = T(ind);
  dT = T1(2:end) - T1(1:end-1);
  ff = interp1(T1(1:end-1),1./dT,T,'nearest','extrap');
  
  
%% Calculate Phase

  % phi(t) = 2*pi*(t-t1)/(t2-t1)
  
  phi = zeros(1,numel(T));
  %ff  = zeros(size(T));
  
  for k=1001:length(T)-1001
     
    % Find preceding R peak
    dT = T(ind)-T(k);
    t2 = min(dT(dT>0)) + T(k);
    
    % Find next R peak
    dT = T(k)-T(ind);
    t1 = T(k) - min(dT(dT>0));
    
    % The phase is
    phi(k) = 2*pi*(T(k)-t1)./(t2-t1);

    % The frequency is
    %ff(k) = 1/(t2-t1);
  end
      
  % Show result
  %plot(phi)

%% Check result

  if (nargout == 0)
    figure(3)
    subplot(121); hold on
      specgram(sin(phi), 16384, 1/0.001);
      plot(T,ff)
      ylim([0 2])
    subplot(122); hold on
      specgram(ecg, 16384, 1/0.001);
      plot(T,ff)
      ylim([0 2])
  end
  