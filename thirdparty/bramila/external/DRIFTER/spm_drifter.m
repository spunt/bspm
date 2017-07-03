function spm_drifter(job)
% SPM_DRIFTER - Configuration file for SPM toolbox
%
% Syntax:
%   spm_drifter(job)
%
% In:
%   job - A valid SPM job structure
%
% Description:
%     This file defines a SPM (Statistical Parametric Mapping) toolbox that
%   can be used for removing periodic noise components from fMRI data, most
%   typically respiration- and cardiac-induced noises.
%     The basic idea of the toolbox is to use externally recorded reference
%   signals of respiration and cardiac activity to estimate the frequency
%   time series of the noise components using an IMM (interacting multiple 
%   model) approach. The frequencies are used in separating the noise-
%   induced components from the fMRI data. The method is originally focused
%   on fast acquisition data with TRs spanning some 100 ... 250 ms.
%     If the TR is small enough, reference signals are not necessarily
%   needed at all, and even if the TR ranges from 500 ... >2000 ms the
%   method can be used as longa as there are reference signals available.
%
% Refrences:
%   Sarkka, S., Solin, A., Nummenmaa, A., Vehtari, A., Auranen, T., 
%   Vanni, S., and Lin, F.-H. (2012). Dynamical retrospective filtering of 
%   physiological noise in BOLD fMRI: DRIFTER. NeuroImage, 60:1517-1527.
%
%   Glover et al. (2000) Image-Based Method for Retrospective Correction 
%   of Physiological Motion Effects in fMRI: RETROICOR. Magnetic Resonance
%   in Medicine, 44, 162-167.
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

%% Report what we are doing

  % This is version:
  version = '20120425';

  % Use SPM style progress reporting (no graphical output)
  fprintf('\n%-40s%+32s\n', ...
    ['DRIFTER v' version ' (arno.solin@aalto.fi)'], ...
    datestr(now,'HH:MM:SS - dd/mm/yyyy'))
  fprintf(['========================================' ...
    '================================\n'])


%% Read volumes

  % Report what we are doing
  fprintf('%-36s:%+35s\n','Reading volumes','Preparing..')
  
  % Read epifiles into a matrix array
  V=spm_vol(job.epidata.files);

  if (~isempty(V) && isfield(V{1},'dim'))
      dims = V{1}.dim;
  else
      error('No files actually loaded. File list empty?') 
  end
  
  % Data writing done
  fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
          'Done.') 
  
%% Do the filtering per slice to save memory
  
  % To avoid allocating all the data into memory at once, we load just
  % one slice at a time and run the DRIFTER/RETROICOR method on that
  % data. Thereafter the results are saved and the data corresponding
  % to the next slice is loaded and filtered.

  for slicenum=1:dims(3)
      
    % Allocate space for the time series of one slice
    job.epidata.data = zeros(dims(1),dims(2),1,length(V));
    
    % Report what is being done
    fprintf('%-36s:%+35s\n', ...
        sprintf('Reading slice %3i/%-3i',slicenum,dims(3)),'Preparing..')
    
    % Loop through all datafiles
    for k=1:length(V)
      
      % Read volumes using SPM functions
      foo=spm_read_vols(V{k});
      
      % Just take the data corresponding the current slice
      job.epidata.data(:,:,1,k) = foo(:,:,slicenum); 
      
      % Report volume number to user
      fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
        sprintf('%i/%i',k,length(V)))
    
    end
  
    % Data reading done
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
            'Done.') 
    
    % Set job details
    job.epidata.dt = job.epidata.tr / 1000;
    job.epidata.mode = job.mode;
    
    % Run the chosen method, mode '2' is for RETROICOR and
    % mode '0' for DRIFTER returning the cleaned BOLD signal only, 
    % and '1' for DRIFTER in RETROICOR-compatible mode, where the
    % measurement noise is added back to the estimate
    if job.mode == 2

      % Check that the dimension of refdata is two:
      if length(job.refdata) ~= 2, 
        error('RETROICOR requires cardiac and respiartory signals.'); 
      end
        
      % Run RETROICOR if chosen
      [out.epidata.estimate]=retroicor(job.epidata.data, ...
          job.refdata(1).data, ...
          job.refdata(2).data, ...
          [job.epidata.dt job.refdata(1).dt job.refdata(2).dt], ...
          job.refdata(1).N, ...
          job.refdata(2).N);
  
    else
      
      % Run estimation using the DRIFTER method
      [out.epidata,out.refdata]=drifter(job.epidata,job.refdata);

      % Use the same scaling for all slices, if not zero
      if (out.epidata.scalefactor ~= 0),
        job.epidata.scalefactor = out.epidata.scalefactor;
      end
      
      % Check results for NaNs and replace by zeros
      out.epidata.estimate(isnan(out.epidata.estimate)) = 0;
      
      % Use the same frequency estimates for all slices (i.e. do not
      % run the IMM step more than once)
      for j=1:length(job.refdata)
        job.refdata(j).frequency = out.refdata{j}.frequency;  
      end      
      
    end
    
    % Get defined prefix for output
    prefix = job.prefix;
  
    % Report what is being done
    fprintf('%-36s:%+35s\n', ...
        sprintf('Writing slice %3i/%-3i',slicenum,dims(3)),'Preparing..')
    
    % Loop through all datafiles
    for k=1:length(V)

      % If this is the first slice, we load the original data and
      % replace the data corresponding to the first slice. Thereafter
      % we change the data filenames to point to thsese new files and
      % the rest of the data will be loaded from the new files.
      if (slicenum==1)
        
        % Read volumes using SPM functions
        foo=spm_read_vols(V{k});
        
        % Modify the filename by adding the prefix
        [pathstr, name, ext] = fileparts(V{k}.fname);
        V{k}.fname = [pathstr '/' prefix name ext];

      else
        % Read volumes using SPM functions
        foo=spm_read_vols(V{k});  
      end
        
      % Only modify the current slice
      foo(:,:,slicenum) = out.epidata.estimate(:,:,:,k);
      
      % Write modified volume to disk 
      spm_write_vol(V{k},foo);
              
      % Report volume number to user
      fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
        sprintf('%i/%i',k,length(V)))
     
    end

    % Data writing done
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%+20s\n', ...
            'Done.') 
    
  end
  
  
%% Save visualization of IMM result to disk

  if (job.visual == 1 && job.mode < 2)
       saveIMMresult(out.epidata,out.refdata);
  end

  
