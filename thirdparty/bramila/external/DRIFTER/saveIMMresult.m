function saveIMMresult(epidata,refdata)
% saveIMMresult - Visualize the IMM estimation result
%
% Syntax:
%   saveIMMresult(epidata,refdata)
%
% In:
%   epidata - Structure
%              (...)
%   refdata - D-dim array of structures
%              .name      = String with reference signal name
%              .dt        = Sampling interval
%              .data      = 1xN vector
%              .downdt    = dt to downsample to
%              .FF        = Estimated frequency time series
%              .freqlist  = List of possible frequencies
%              
% Out:
%   ...     - Saves plot to disk.
%
% Description:
%   Save figures of reference signals to disk.
%

% Copyright:
%   Arno Solin (2011)
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

%% Visualize refernces and IMM estimation results
% 
%     % New figure handle
%     handle=figure; clf
%     
%     % Show each reference
%     for i=1:length(refdata)
%         
%       % Is there any data
%       if ~isfield(refdata{i},'data'), continue; end
%       
%       % Plot actual signal and estimate
%       subplot(2,length(refdata),i); hold on
%         
%         T = 0:refdata{i}.dt:refdata{i}.dt*(length(refdata{i}.data)-1);
%         plot(T,refdata{i}.data,'Color',[.7 .7 .7])
%         T = 0:refdata{i}.downdt:refdata{i}.downdt*(length(refdata{i}.downsampled)-1);
%         plot(T,refdata{i}.S+sum(refdata{i}.CS,1),'-r')
%         plot(T,refdata{i}.S,'-b')
%         
%         xlabel('Time [s]')
%         title(sprintf('\\bf Data %s',refdata{i}.name),'FontSize',9)
%         lgnd=legend(refdata{i}.name,'Estimate', ...
%             'Bias without periodic signal', 'Location','best');
%         set(lgnd,'FontSize',6)
%         axis tight; box on
%       
%       % Plot spectrogram with estimates
%       subplot(2,length(refdata),length(refdata)+i)
%         hold on
%         specgram(refdata{i}.downsampled, 256, 1/refdata{i}.downdt)
%         T = 0:epidata.dt:epidata.dt*(length(refdata{i}.FF)-1);
%         plot(T,refdata{i}.FF/60,'b','LineWidth',2);
%         hold off
%         lgnd=legend('IMM Frequency Estimate');
%         set(lgnd,'FontSize',6)
%         ylabel('Frequency [1/s]'); xlabel('Time [s]')
%         title(sprintf('\\bf Spectrogram of %s',refdata{i}.name),'FontSize',9)
%         axis tight; box on
%         ylim([refdata{i}.freqlist(1)-10 refdata{i}.freqlist(end)+10]/60);
%     end
% 
%   % New filename
%   filename = sprintf('%s.%s', datestr(now,30),'png');
% 
%   % Save figure as png
%   saveas(handle,filename,'png')
% 
%   % Close figure handle
%   close(handle)
%   
%   
%% Further figures

  % New figure handle
  handle=figure; clf  
  
  % Show each reference
  for i=1:length(refdata)
        
    % Is there any data
    if ~isfield(refdata{i},'data'), continue; end
            
    % Plot spectrogram with estimates
    subplot(1,length(refdata),i); hold on
        
      % Plot spectrogram
      specgram(refdata{i}.data,1024*16,1/refdata{i}.dt)
        
      % Plot IMM frequency estimate
      T = 0:epidata.dt:epidata.dt*(length(refdata{i}.FF)-1);
      plot(T,refdata{i}.FF/60,'-k','LineWidth',2);
        
      % Modify plot
      title(refdata{i}.name, 'Interpreter','none','FontWeight','bold')
      ylabel('Frequency [1/s]'); xlabel('Time [s]')
      ylim([0 4]); xlim([min(T) max(T)])
      box on, set(gca,'layer','top')
      
  end
  
  % Set paper size and options
  set(handle,'PaperUnits','inches');
  set(handle,'PaperSize',[10 3])
  set(handle,'PaperPosition',[0.0 0.0 10 3])
  set(handle,'Color','w')
  
  % New filename
  filename = sprintf('%s.%s', datestr(now,30),'png');
  
  % Save figure as png
  saveas(handle,filename,'png') %'epsc2')
  
  % Close figure handle
  close(handle)


