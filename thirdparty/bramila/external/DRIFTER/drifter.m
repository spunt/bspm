function [data,refdata,SMM,SPP,FMM,FPP] = drifter(data,refdata)
% DRIFTER - Estimate model with the DRIFTER method
%
% Syntax:
%   [data,refdata,SMM,SPP,FMM,FPP] = drifter(data,refdata)
%
% In:
%   data    - Structure
%              .dt = sampling interval
%              .data = xdim x ydim x numslices x N matrix
%              .sd = (optional) Noise standard deviation
%              .BF = (optional) Feedback matrix
%              .BQ = (optional) Process spectral density for the bias model
%              .BL = (optional) Noise multiplier matrix for the bias model
%              .BH = (optional) Measurement matrix for the bias model
%              .mode = (optional) Estimation mode for SPM output
%   refdata - D-dim array of structures
%              .name      = (optional) String with reference signal name
%              .dt        = sampling interval
%              .data      = (optional) 1xN vector
%              .N         = (optional) Number of harmonics to be estimated
%              .Nimm      = (optional) The same, but can be different in IMM
%              .downdt    = (optional) dt to downsample to
%              .freqlist  = (optional) List of possible frequencies for IMM
%                            in beats/min.
%              .frequency = (optional) 1xN vector of frequencies instead of
%                            using IMM
%              .ptrans    = (optional) IMM transition probability
%              .poverall  = (optional) IMM jump anywhere transition probability
%              .sd = (optional) Noise standard deviation
%              .BF = (optional) Feedback matrix
%              .BQ = (optional) Process spectral density for the bias model
%              .BL = (optional) Noise multiplier matrix for the bias model
%              .BH = (optional) Measurement matrix for the bias model
%              .qr = (optional) Resonator's process noise spectral density.
%              .filter = 1/0 if 0 - only bandpass filter, 1 is for amplitude normalisation 
% Out:
%   data    - Structure
%              [the ones specified above with defaults added and...]
%              .estimate = xdim x ydim x numslices x N cleaned signal 
%              .noise    = noise estimate
%   refdata - D-dim array of structures
%              [the ones specified above with defaults added and...]
%              .estimate = xdim x ydim x numslices x N signal in EPI 
%              .FF  = frequency with sampling matching the EPI data
%              .S   = Cleaned signal
%              .CS  = Fundamental signal and its harmonics
%              .QCS = Quadrature periodic signals 
%   SMM - Smoother output mean matrix (returned if size of data small enough)
%   SPP - Smoother output covariance matrix (returned if size of data small enough)
%   FMM - Filter output mean matrix (returned if size of data small enough)
%   FPP - Filter output covariance matrix (returned if size of data small enough)
%
% Description:
%   Clean signal using the DRIFTER method. The estimation is split into
%   multiple passes if the size of the data is larger than 64x64x32x100.
%   In this case the matrices SMM, SPP, FMM and FPP are not returned.
%
% References:
%   Sarkka, S., Solin, A., Nummenmaa, A., Vehtari, A., Auranen, T., 
%   Vanni, S., and Lin, F.-H. (2012). Dynamical retrospective filtering of 
%   physiological noise in BOLD fMRI: DRIFTER. NeuroImage, 60:1517?1527.
%   
% Update history: (original version 20110526)
%   20110809 - Support for large matrices by splitting into passes
%   20110817 - Modified warning message not to be shown if N=1
%   20110829 - Fixed bug with wrong values being returned
%   20110901 - Reference data not required anymore
%   20110930 - Fixed indexing bug in smoother
%   20120313 - Fixed bug when running without SPM
%   20120417 - The method now accepts complex-valued data.
%              Changed the way data is loaded in SPM (support for larger
%              matrix sizes without running out of memory).
%              Empty matrices now interpreted as unset.
%              Performs an update check on each run.
%              Accepts 1D data for filtering.
%   20120425 - Fixed broken RETROICOR support. 
%
% Copyright:
%   Arno Solin, 2011-2012
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
  
%% Report version and check online for updates

  % This is version:
  version = '20120425';
  
  % Check for updates
  try
    current_version = urlread(...
      'http://becs.aalto.fi/en/research/bayes/drifter/version.txt');
    current_version = str2num(current_version);  
  catch
    current_version = nan;
  end
  
  % Report version and update status
  if current_version > str2num(version)
    fprintf('%-36s:%+35s\n',['DRIFTER version: ' version], ...
        'New version available online!')
  elseif isnan(current_version)
    fprintf('%-36s:%+35s\n',['DRIFTER version: ' version], ...
        'Update check failed.')
  else
    fprintf('%-36s:%+35s\n',['DRIFTER version: ' version], ...
        'Up to date.')
  end

  
%% Convert inputs to cell arrays

  % Check if any inputs have been assigned in the first place
  if nargin ~= 2 , error('Inputs arguments not specified.'); end
  
  % We accept both structure arrays and cell arrays to be submitted to
  % DRIFTER. This is because SPM does not support cell arrays. All
  % the refdata elements are converted to cell arrays to ease the
  % rest of the implementation.

  % Convert to cell arrays
  if ~iscell(refdata)
    foo=num2cell(refdata);
  else
    foo=refdata;  
  end

  % The user might have left values unassigned or empty to revert to
  % the defaults. We only check for unassigned field, so we remove 
  % all fields with empty values.
  for j=1:length(foo)
    fnames = fieldnames(foo{j});
    for i=1:length(fnames)
      if isempty(getfield(foo{j},fnames{i}))   
        foo{j} = rmfield(foo{j},fnames{i});
      end
    end
  end
  refdata = foo;

  % We do the same for the 'data' structure
  foo = data;
  fnames = fieldnames(data);
  for i=1:length(fnames)
    if isempty(getfield(foo,fnames{i}))   
        foo{j} = rmfield(foo{j},fnames{i});
    end
  end
  data = foo;

  
%% Check input values and set defaults

  % We check each input value and if they are left unassigned
  % (or empty) we assign default values to the fields.
    
  if ~isfield(data,'data'), 
     error('Data field in ''data'' required.') 
  end
   
  % If the data is a one-dimensional vectorm, we pad it with
  % singleton dimensions and make it four-dimensional
  if length(size(data.data)) == 2 && ...
     numel(data.data) == max(size(data.data))
       data.data = reshape(data.data,1,1,1,numel(data.data));
  end
      
  if length(size(data.data)) ~= 4,
      error(['Dimension of the ''data'' field in ''data'' has to be 4 ' ...
               '(xdim x ydim x numslices x N).']);
  end
    
  if ~isfield(data,'dt'), 
      error('Sampling interval dt required in data.')
  end
    
  if ~isfield(data,'sd'), data.sd = 0.1; end
  if ~isfield(data,'BF'), data.BF = [0 1; 0 0]; end
  if ~isfield(data,'BQ'), data.BQ = 0.01;       end
  if ~isfield(data,'BL'), data.BL = [0;1];      end
  if ~isfield(data,'BH'), data.BH = [1 0];      end
  if ~isfield(data,'qr'), data.qr = 0.01;       end
  
  % Set default matrix size limit (split up if larger)
  if ~isfield(data,'limit'), data.limit = 64*64*32*100; end
  
  % Set default mode to return all components
  if ~isfield(data,'mode'), data.mode = -1;     end
  
  % Check each reference signal for all parameters. The only required
  % fields are dt and data/frequency. To all others default values are
  % applied if they are unspecified.

  % Apply defaults or throw errors for each reference/field
  for i=1:length(refdata)

    % Use the average fMRI signal as data if no data supplied
    if ~isfield(refdata{i},'frequency') && ~isfield(refdata{i},'data'), 
        
        % Remove data mean and reshape signal
        meansignal = bsxfun(@minus,data.data,mean(data.data,4));
        dims       = size(meansignal);
        meansignal = reshape(meansignal,prod(dims(1:3)),[]);
        
        % Discard parts with stdev < median and with highest 10% stdev 
        % meansignal = meansignal(std(meansignal,[],2) > median(std(meansignal,[],2)),:);
        % meansignal = meansignal(std(meansignal,[],2) < 0.9*max(std(meansignal,[],2)),:);
        
        % Average and set
        refdata{i}.data   = mean(meansignal,1);
        refdata{i}.dt     = data.dt;
        refdata{i}.downdt = data.dt;

        % Report what we have done
        fprintf('%-36s:%+35s\n',sprintf('Composing data for reference %i',i),'Done.')
        clear meansignal
    end
      
    if ~isfield(refdata{i},'name'), 
        refdata{i}.name = sprintf('Refrence_signal_%i',i); 
    end
    
    if isfield(refdata{i},'data') && ~isfield(refdata{i},'dt'), 
        error('Sampling interval dt required in refrence %i.',i)
    end
 
    if ~isfield(refdata{i},'N'), 
        refdata{i}.N = 1; % Apply default
    end
    
    if ~isfield(refdata{i},'Nimm'), 
        refdata{i}.Nimm = refdata{i}.N; % Apply default
    end
    
    if ~isfield(refdata{i},'downdt'), 
        refdata{i}.downdt = refdata{i}.dt; % No downsampling
    end
    
    if ~isfield(refdata{i},'freqlist') && isfield(refdata{i},'data'), 
        error('Freqlist required in reference %i.',i);
    end
    
    if ~isfield(refdata{i},'ptrans'), 
        refdata{i}.ptrans = 0.01; % Apply default
    end
    
    if ~isfield(refdata{i},'poverall'), 
        refdata{i}.poverall = 0; % Apply default
    end
    
    if ~isfield(refdata{i},'sd'), refdata{i}.sd = 0.1;        end
    if ~isfield(refdata{i},'BF'), refdata{i}.BF = [0 1; 0 0]; end
    if ~isfield(refdata{i},'BQ'), refdata{i}.BQ = 0.01;       end
    if ~isfield(refdata{i},'BL'), refdata{i}.BL = [0;1];      end
    if ~isfield(refdata{i},'BH'), refdata{i}.BH = [1 0];      end
    if ~isfield(refdata{i},'qr'), refdata{i}.qr = 0.01;       end
    
    % Check if the reference data is complex-valued. If so, discard 
    % the imaginary part and warn the user about this.
    if isfield(refdata{i},'data') && sum(abs(imag(refdata{i}.data))) > 1e-9
       warning(['Reference %i is complex-valued. Discarding the ' ...
           'imaginary part in frequency estimation.'],i)
       refdata{i}.data = real(refdata{i}.data);
    end
    
    % Check if Nyquist frequency is met during IMM
    if (isfield(refdata{i},'freqlist') && ...
       refdata{i}.Nimm*max(refdata{i}.freqlist)/60 > 1/refdata{i}.downdt/2)
       warning(['In %s: maximum freqlist frequency %.2f Hz with N=%i ' ...
           'harmonics might hit Nyqvist frequency at %.2f Hz. ' ...
           'Consider changing N to %i.'], ...
           refdata{i}.name,max(refdata{i}.freqlist)/60, ...
           refdata{i}.Nimm, 1/refdata{i}.downdt/2, ...
           floor(60/max(refdata{i}.freqlist)/refdata{i}.downdt/2))
    end
  end

  
%% Run reference frequency estimation

  % This is the first actual part of the DRIFTER method. We use IMM to
  % find the optimal frequency trajectory of the oscillators with 
  % respect to a given frequency interval (i.e. freqlist). If the
  % 'frequency' field is supplied by the user, the IMM step is not
  % needed and the frequency is only interpolated to the same time
  % instants as the fMRI data.

  % For each reference signal we do the IMM step
  for i=1:length(refdata)
      
    % Use IMM if we have data else we should already have the frequencies
    if isfield(refdata{i},'data') && ~isfield(refdata{i},'frequency')
      
      % Downsample data
      T  = 0:refdata{i}.dt:refdata{i}.dt*(length(refdata{i}.data)-1);
      Ti = 0:refdata{i}.downdt:T(end);
      refdata{i}.downsampled = interp1(T,refdata{i}.data,Ti,'linear');
      % filtering 0=only BP filter, 1=normalisation
      if isfield(refdata{i},'filter')
      refdata{i}.downsampled = drifterfilter(refdata{i});
      end
      
      % Find mean and scale factors if not supplied by user
      if ~isfield(refdata{i},'mean')
        refdata{i}.mean = mean(refdata{i}.downsampled);
      end
      if ~isfield(refdata{i},'scalefactor')
        refdata{i}.scalefactor = ...
            std(refdata{i}.downsampled(floor(end/2):end));
      end
      
      % Remove mean and normalize scale
      refdata{i}.downsampled = (refdata{i}.downsampled-refdata{i}.mean)./ ...
          refdata{i}.scalefactor;
      
      % Run the IMM tracking
      [refdata{i}.frequency] = periodic_track_imm( ...
                          refdata{i}.downsampled, ...
                          refdata{i}.downdt,      ...
                          refdata{i}.freqlist,    ...
                          refdata{i}.Nimm,        ...
                          refdata{i}.BF,          ...
                          refdata{i}.BQ,          ...
                          refdata{i}.BL,          ...
                          refdata{i}.BH,          ...
                          refdata{i}.sd^2,        ...
                          refdata{i}.qr,          ...
                          refdata{i}.ptrans,      ...
                          refdata{i}.poverall);
                    
      % Run single filter and smoother with the MMSE frequency
      [S,CS,QCS] = periodic_separation_kfs( ...
                          refdata{i}.downsampled, ...
                          refdata{i}.downdt,      ...
                          refdata{i}.frequency,   ...
                          refdata{i}.Nimm,        ...
                          refdata{i}.BF,          ...
                          refdata{i}.BQ,          ...
                          refdata{i}.BL,          ...
                          refdata{i}.BH,          ...
                          refdata{i}.sd^2,        ...
                          refdata{i}.qr);      
      
      % Scale back and store for output
      refdata{i}.S = (S*refdata{i}.scalefactor) + refdata{i}.mean;
      refdata{i}.CS = CS*refdata{i}.scalefactor;
      refdata{i}.QCS = QCS*refdata{i}.scalefactor;
                      
    end
    
    % Sample frequency to match the size of the actual data
    T  = 0:refdata{i}.downdt:refdata{i}.downdt*(length(refdata{i}.frequency)-1);
    Ti = 0:data.dt:data.dt*(size(data.data,4)-1);
    refdata{i}.FF = interp1(T,refdata{i}.frequency,Ti,'linear','extrap');

    % Check if Nyqvist frequency will be met in the actual estimation,
    % and throw a warning to the user is this is the case. This is only
    % to make the use more convenient for the less experienced users
    if (refdata{i}.N > 1 && ...
      refdata{i}.N*max(refdata{i}.FF)/60 > 1/data.dt/2)
       warning(['When doing the estimation in the EPI signal, the ' ...
           'maximum estimated frequency %.2f Hz with N=%i ' ...
           'harmonics hits Nyqvist frequency at %.2f Hz. ' ...
           'Consider changing N to %i for %s.'], ...
           max(refdata{i}.FF)/60, refdata{i}.N, 1/data.dt/2, ...
           max([floor(60/max(refdata{i}.FF)/data.dt/2) 1]), ...
           refdata{i}.name)
    end
    
  end

  
%% Plot estimated frequencies with spectrogram

  % Plot only if no output arguments given
  %if (nargout == 0)
  %always plot
  
   h=figure('Visible', 'off');
    for i=1:length(refdata)
      % Is there any data
      if ~isfield(refdata{i},'data'), continue; end
      
      % Plot actual signal and estimate

      subplot(2,length(refdata),i); hold on
        
        T = 0:refdata{i}.dt:refdata{i}.dt*(length(refdata{i}.data)-1);
        plot(T,refdata{i}.data,'Color',[.7 .7 .7])
        T = 0:refdata{i}.downdt:refdata{i}.downdt*(length(refdata{i}.downsampled)-1);
        plot(T,refdata{i}.S+sum(refdata{i}.CS,1),'-r')
        plot(T,refdata{i}.S,'-b')
        
        xlabel('Time [s]')
        title(sprintf('\\bf Data %s',refdata{i}.name),'FontSize',14)
        legend(refdata{i}.name,'Estimate','Bias without periodic signal')
        axis tight; box on
      
      % Plot spectrogram with estimates
      subplot(2,length(refdata),length(refdata)+i)
        hold on
        specgram(refdata{i}.downsampled, 256, 1/refdata{i}.downdt)
        T = 0:data.dt:data.dt*(length(refdata{i}.FF)-1);
        plot(T,refdata{i}.FF/60,'b','LineWidth',2);
        hold off
        legend('IMM Frequency Estimate')
        ylabel('Frequency [1/s]'); xlabel('Time [s]')
        title(sprintf('\\bf Spectrogram %s',refdata{i}.name),'FontSize',14)
        axis tight; box on
        ylim([refdata{i}.freqlist(1)-10 refdata{i}.freqlist(end)+10]/60);
    end
    [pathstr,name,~] = fileparts(refdata{i}.name);
    saveas(h,[pathstr,'/refdata' num2str(i) '_signal_estimate - spectrogram_' datestr(clock)],'jpg');
  %end

  
%% Normalize the fMRI original data

  % Remove mean and save it for later use if not supplied by the user
  if ~isfield(data,'mean')
    data.mean = mean(data.data,4);
  end
  data.data = bsxfun(@minus, data.data, data.mean); 
  
  % We have to separate scale factors for the real and imaginary parts
  if ~isfield(data,'scalefactor')
    data.scalefactor = std(real(data.data(floor(end/2):end))) + ...
        1i*std(imag(data.data(floor(end/2):end)));
  end
  
  % Normalize scale. This makes choosing the paramteres simpler. 
  if abs(imag(data.scalefactor)) > 1e-9
    data.data = real(data.data)./real(data.scalefactor) + ...
        1i*imag(data.data)./imag(data.scalefactor);
  else
    data.data = data.data./data.scalefactor;  
  end
  
  % The number of time steps is
  steps = size(data.data,4);
  T = 0:data.dt:data.dt*(steps-1);
  
  
%% Do the 3D linear Kalman filter estimation with known frequency

  % This code runs the linear Kalman filter with the help of the estimated
  % frequency and separates the original signal into a (i) brain signal,
  % (ii) periodic components, and (iii) noise

  % Clear old values / allocate space
  data.FF = zeros(steps,length(refdata));
  data.N  = zeros(1,length(refdata));
  data.qr = zeros(1,length(refdata));
  
  % We have frequency estimates from the IMM method
  for i=1:length(refdata)
    data.FF(:,i) = refdata{i}.FF;
    data.N(i)    = refdata{i}.N;
    data.qr(i)   = refdata{i}.qr;
  end
   
  % Split into multiple passes if the dataset is too large
  npass = min([size(data.data,1), ceil(numel(data.data)/data.limit)]);
  
  % Allocate space for results
  data.estimate = zeros(size(data.data));
  
  % If mode = -1, all components are returned
  if (data.mode < 0)
    for i=1:length(refdata)
      refdata{i}.estimate = zeros(size(data.data));
    end
  end
  
  % Run the passes
  for i=1:npass
    
    % The indices in this pass
    inds = ceil(linspace(1,size(data.data,1),npass+1)); 
    inds(end)=inds(end)+1;
    inds = inds(i):inds(i+1)-1;
    
    % Run single filter and smoother with the MMSE frequency
    [SMM,SPP,FMM,FPP] = periodic_separation_kfs3D( ...
        data.data(inds,:,:,:), ...
        data.dt,   ...
        data.FF,   ...
        data.N,    ...
        data.BF,   ...
        data.BQ,   ...
        data.BL,   ...
        data.BH,   ...
        data.sd^2, ...
        data.qr);
    
    % The estimate for the brain (bold) signal
    dims = size(data.data); dims(1) = numel(inds);
    H=zeros(1,2*sum(data.N)+size(data.BH,2)); 
    H(end-size(data.BH,2)+1:end) = data.BH;
    MULTIPLIER = kron(speye(prod(dims(1:3))),H);
    
    % This gets a bit messy, because we try to avoid large marices being
    % allocated if they are not needed
    
    % (i) Everything is returned
    if (data.mode < 0) 
      data.estimate(inds,:,:,:) = reshape(MULTIPLIER*SMM,dims);
    
      % The periodic noise estimates
      for i=1:length(refdata)
      
        % Form the multiplier matrix
        H=zeros(1,2*sum(data.N)+size(data.BH,2));
        H(2*sum(data.N(1:i-1))+(1:2:2*data.N(i))) = 1;
        MULTIPLIER = kron(speye(prod(dims(1:3))),H);
      
        refdata{i}.estimate(inds,:,:,:) = ...
            real(reshape(MULTIPLIER*SMM,dims))*real(data.scalefactor) + ...
            1i*imag(reshape(MULTIPLIER*SMM,dims))*imag(data.scalefactor); 
      end
      
    % The code is called from SPM and only the relevant data is returned
    % (ii) The noise-free bias is returned 
    elseif (data.mode == 0)  
      data.estimate(inds,:,:,:) = reshape(MULTIPLIER*SMM,dims);

    % The code is called from SPM and only the relevant data is returned
    % (iii) The method runs in RETROICOR-compatible mode
    else
        
      % Start with the original data
      data.estimate(inds,:,:,:) = data.data(inds,:,:,:);    
        
      % The periodic noise estimates
      for i=1:length(refdata)
      
        % Form the multiplier matrix
        H=zeros(1,2*sum(data.N)+size(data.BH,2));
        H(2*sum(data.N(1:i-1))+(1:2:2*data.N(i))) = 1;
        MULTIPLIER = kron(speye(prod(dims(1:3))),H);
      
        data.estimate(inds,:,:,:) = data.estimate(inds,:,:,:) - ...
            reshape(MULTIPLIER*SMM,dims);
      end
    end
    
    % If there are multiple passes, do not return the matrices
    if (npass > 1 && data.mode > -1)
      SMM = []; SPP = []; FMM = []; FPP = [];  
    end
  end
  
  % The output result is now in the matrices SMM (smoother mean) and FMM 
  % (filter mean). The data is stored as vectors [(2*N+2)*X*Y] x numslices
  % x [number of timesteps], where 2*N first terms for each voxel are the
  % periodic components (component and derivative) and the last two values
  % the bold signal and its derivative. X and Y are the resolution of the
  % layer (eg. 64 x 64) and here we use a number of [numslices] layers.
  
  
%% Resize, rescale and return
  
  % Rescale and move
  data.estimate = real(data.estimate)*real(data.scalefactor) + ...
      1i*imag(data.estimate)*imag(data.scalefactor);
  data.estimate = bsxfun(@plus, data.estimate, data.mean);
  
  % Rescale and move original data
  data.data = real(data.data)*real(data.scalefactor) + ...
      1i*imag(data.data)*imag(data.scalefactor);
  data.data = bsxfun(@plus, data.data, data.mean);
  
  % Also return the noise estimate, if needed
  if (data.mode < 0)
  
    % Allocate noise matrix
    data.noise = data.data - data.estimate;
  
    % The reference signals
    for i=1:length(refdata)
      
      % Change in noise
      data.noise = data.noise - refdata{i}.estimate;
    
    end
  end
    
