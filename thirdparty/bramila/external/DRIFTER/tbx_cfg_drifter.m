function drifter = tbx_cfg_drifter
% tbx_cfg_drifter - Configuration file for the SPM toolbox
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
%     Sarkka, S., Solin, A., Nummenmaa, A., Vehtari, A., Auranen, T., 
%   Vanni, S., and Lin, F.-H. (2012). Dynamical retrospective filtering of 
%   physiological noise in BOLD fMRI: DRIFTER. NeuroImage, 60:1517-1527.
%     Glover et al. (2000) Image-Based Method for Retrospective Correction 
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

%% Specify the reference data parameters

% Reference data: name
  name              = cfg_entry;
  name.tag          = 'name';
  name.name         = 'Reference name';
  name.help         = {'Reference signal name (e.g. ''External cardiac'' ' ...
                       'signal). This iformation can be provided for ' ...
                       'debugging puposes.'};
  name.val          = {'Refrence signal'};
  name.strtype      = 's';
  name.num          = [1 Inf];
% Reference data: dt
  dt                = cfg_entry;
  dt.tag            = 'dt';
  dt.name           = 'Sampling interval';
  dt.help           = {'(Required) The sampling interval, dt, of the '...
                       'reference signal (in seconds).'};
  dt.val            = {[]};
  dt.strtype        = 'r';
  dt.num            = [0 Inf];
% Reference data: data
  data              = cfg_entry;
  data.tag          = 'data';
  data.name         = 'Reference signal data';
  data.help         = {'Reference signal data as a 1xN vector ' ...
                       'sampled the specified resolution and spanning ' ...
                       'the same time interval as the (fMRI) scans.' ...
                       'If no data is specified, the average fMRI signal is used.'};
  data.val          = {[]};
  data.strtype      = 'r';
  data.num          = [0 Inf];
% Reference data: N
  N                 = cfg_entry;
  N.tag             = 'N';
  N.name            = 'Number of periodics';
  N.help            = {'(Required) The total number of periodics (fundamental + ' ...
                       'harmonics) that are estimated from the fMRI data.'};
  N.val             = {1};
  N.strtype         = 'r';
  N.num             = [1 1];
% Reference data: Nimm
  Nimm              = cfg_entry;
  Nimm.tag          = 'Nimm';
  Nimm.name         = 'Number of periodics in IMM';
  Nimm.help         = {'The total number of periodics (fundamental + ' ...
                       'harmonics) that are estimated during the IMM stage ' ...
                       'where the frequency trajectory is estimated from ' ...
                       'the reference signal.'};
  Nimm.val          = {1};
  Nimm.strtype      = 'r';
  Nimm.num          = [1 1];
% Reference data: downdt
  downdt            = cfg_entry;
  downdt.tag        = 'downdt';
  downdt.name       = 'Sampling interval to downsample to';
  downdt.help       = {'The sampling interval (in seconds) to downsample ' ...
                       'the reference signal to. Considerable speed-up can ' ...
                       'be achieved in this way as often the fMRI TR >> ' ...
                       'the reference dt.'};
  downdt.val        = {1/10};
  downdt.strtype    = 'r';
  downdt.num        = [1 1];
% Reference data: freqlist
  freqlist          = cfg_entry;
  freqlist.tag      = 'freqlist';
  freqlist.name     = 'Array of possible frequencies in bpm';
  freqlist.help     = {'An array of possible frequencies in beats per ' ...
                       'minute (e.g. 40:120). The interacting Markov ' ...
                       'chain model uses these for finding the ' ...
                       'frequency trajectory in the external signal.'};
  freqlist.val      = {[]};
  freqlist.strtype  = 'r';
  freqlist.num      = [0 Inf];
% Reference data: frequency
  frequency         = cfg_entry;
  frequency.tag     = 'frequency';
  frequency.name    = 'Frequency time series';
  frequency.help    = {'Alternatively a frequency time series spanning ' ...
                       'over the same time interval as the scans can be ' ...
                       'defined.'};
  frequency.val     = {[]};
  frequency.strtype = 'r';
  frequency.num     = [0 Inf];
% Reference data: ptrans
  ptrans            = cfg_entry;
  ptrans.tag        = 'ptrans';
  ptrans.name       = 'Markov chain transition probability';
  ptrans.help       = {'Transition probability between ' ...
                       'consecutive steps of frequencies (i.e. the ' ...
                       'probability of a jump from e.g. 70 bpm to 71 bpm).'};
  ptrans.val        = {0.01};
  ptrans.strtype    = 'r';
  ptrans.num        = [1 1];
% Reference data: poverall
  poverall          = cfg_entry;
  poverall.tag      = 'poverall';
  poverall.name     = 'Overall transition probability';
  poverall.help     = {'Transition probability between ' ...
                       'all steps of frequencies (i.e. the ' ...
                       'probability of a jump from e.g. 70 bpm to 80 bpm).'};
  poverall.val      = {0};
  poverall.strtype  = 'r';
  poverall.num      = [1 1];
% Reference data: sd
  sd                = cfg_entry;
  sd.tag            = 'sd';
  sd.name           = 'Noise standard deviation';
  sd.help           = {'Estimate of measurement noise standard deviation ' ...
                       'in the IMM model.'};
  sd.val            = {0.1};
  sd.strtype        = 'r';
  sd.num            = [Inf Inf]; % Should be [1 1]
% Reference data: BF
  BF                = cfg_entry;
  BF.tag            = 'BF';
  BF.name           = 'Feedback matrix';
  BF.help           = {'Feedback matrix for the IMM bias model of the external ' ...
                       'signal. The default value corresponds to a Wiener ' ...
                       'velocity model.'};
  BF.val            = {[0 1; 0 0]};
  BF.strtype        = 'r';
  BF.num            = [Inf Inf];
% Reference data: BQ
  BQ                = cfg_entry;
  BQ.tag            = 'BQ';
  BQ.name           = 'Process spectral density for the bias model';
  BQ.help           = {'The process spectral density for the bias model ' ...
                       'represents the continuous time noise in the bold ' ...
                       'signal'};
  BQ.val            = {0.01};
  BQ.strtype        = 'r';
  BQ.num            = [Inf Inf];
% Reference data: BL
  BL                = cfg_entry;
  BL.tag            = 'BL';
  BL.name           = 'Noise multiplier matrix for the bias model';
  BL.help           = {'The noise multiplier matrix for the bias model defines ' ...
                       'how the noise affects the bias evolution.'};
  BL.val            = {[0;1]};
  BL.strtype        = 'r';
  BL.num            = [Inf Inf];
% Reference data: BH
  BH                = cfg_entry;
  BH.tag            = 'BH';
  BH.name           = 'Measurement matrix for the bias model';
  BH.help           = {'The measurement matrix for the bias model defines the ' ...
                       'observation of the bold signal bias process.'};
  BH.val            = {[1 0]};
  BH.strtype        = 'r';
  BH.num            = [Inf Inf];
% Reference data: qr
  qr                = cfg_entry;
  qr.tag            = 'qr';
  qr.name           = 'Resonator process noise spectral density';
  qr.help           = {'The resonator process noise spectral density defines ', ...
                       'the continuous-time variation of the resonator signals. ', ...
                       'Adjust primarily this parameter to control the behavior ', ...
                       'of the periodic signals.'};
  qr.val            = {0.01};
  qr.strtype        = 'r';
  qr.num            = [Inf Inf]; % Should be [1 1]

% Reference branch setup
  data1             = cfg_branch;
  data1.tag         = 'refdata';
  data1.name        = 'Reference Signal';
  data1.val         = {name dt data freqlist frequency N Nimm downdt ...
                       ptrans poverall sd BF BQ BL BH qr};
  data1.help        = {'Specify the reference signal setup for each ' ...
                       'periodic noise component (e.g. pulse and respiration).'};
             
% Array of reference branches
  refdata           = cfg_repeat;
  refdata.tag       = 'refdata';
  refdata.name      = 'Reference Signals';
  refdata.help      = {'Specify reference signal -- e.g. pulse or respiration.'};
  refdata.values    = {data1};
  refdata.num       = [0 Inf];


%% Specify epidata section

% EPI data: dt/TR
  epiTR             = cfg_entry;
  epiTR.tag         = 'tr';
  epiTR.name        = 'TR';
  epiTR.help        = {'(Required) Interscan interval, TR, (specified in ' ...
                       'milliseconds). This is the time between the plane ' ...
                       'scans (typically 100-4000 ms).'};
  epiTR.val         = {};
  epiTR.strtype     = 'n';
  epiTR.num         = [1 1];
% EPI data: files
  files             = cfg_files;
  files.tag         = 'files';
  files.name        = 'EPI files';
  files.help        = {'(Required) Select the (fMRI) scans for this ' ...
                       'session. They all must have identical dimensions/' ...
                       'orientation/voxel size etc.'};
  files.filter      = 'image';
  files.ufilter     = '.*';
  files.num         = [1 Inf];
% EPI data: sd
  episd             = cfg_entry;
  episd.tag         = 'sd';
  episd.name        = 'Noise standard deviation';
  episd.help        = {'Estimate of measurement noise standard deviation.'};
  episd.val         = {0.1};
  episd.strtype     = 'r';
  episd.num         = [1 1];
% EPI data: BF
  epiBF             = cfg_entry;
  epiBF.tag         = 'BF';
  epiBF.name        = 'Feedback matrix';
  epiBF.help        = {'Feedback matrix for the bias model of the bold ' ...
                       'signal. The default value corresponds to a Wiener ' ...
                       'velocity model which assumes that the signal try to' ...
                       'continue in the same direction as on the previous ' ...
                       'step.'};
  epiBF.val         = {[0 1; 0 0]};
  epiBF.strtype     = 'r';
  epiBF.num         = [Inf Inf];
% EPI data: BQ
  epiBQ             = cfg_entry;
  epiBQ.tag         = 'BQ';
  epiBQ.name        = 'Process spectral density';
  epiBQ.help        = {'The process spectral density for the bias model ' ...
                       'represents the continuous time noise in the bold ' ...
                       'signal'};
  epiBQ.val         = {0.01};
  epiBQ.strtype     = 'r';
  epiBQ.num         = [Inf Inf];
% EPI data: BL
  epiBL             = cfg_entry;
  epiBL.tag         = 'BL';
  epiBL.name        = 'Noise multiplier matrix for the bias model';
  epiBL.help        = {'The noise multiplier matrix for the bias model defines ' ...
                       'how the noise affects the bias evolution.'};
  epiBL.val         = {[0;1]};
  epiBL.strtype     = 'r';
  epiBL.num         = [Inf Inf];
% EPI data: BH
  epiBH             = cfg_entry;
  epiBH.tag         = 'BH';
  epiBH.name        = 'Measurement matrix';
  epiBH.help        = {'The measurement matrix for the bias model defines the ' ...
                       'observation of the bold signal bias process.'};
  epiBH.val         = {[1 0]};
  epiBH.strtype     = 'r';
  epiBH.num         = [Inf Inf];
  
  
% EPI data branch setup
  epidata           = cfg_branch;
  epidata.tag       = 'epidata';
  epidata.name      = 'Data';
  epidata.val       = {files epiTR episd epiBF epiBQ epiBL epiBH};
  epidata.help      = {'Specify the fMRI data.'};
  
  
%% Common settings

% Prefix for output files
  prefix            = cfg_entry;
  prefix.tag        = 'prefix';
  prefix.name       = 'Prefix';
  prefix.val        = {'n'};
  prefix.help       = {'File name prefix for output files.'};
  prefix.strtype    = 's';
  
% Mode either with noise of without or with RETROICOR
  mode              = cfg_menu;
  mode.tag          = 'mode';
  mode.name         = 'Output mode';
  mode.val          = {0};
  mode.labels       = {'0: DRIFTER (BOLD)'
                       '1: DRIFTER (BOLD+Noise)'
                       '2: RETROICOR'};
  mode.values       = {0 1 2};
  mode.help         = {'Choose data output mode. Either only the cleaned bias ' ...
                       'BOLD signal is returned or alternatively also the ' ...
                       'noise is added back to the output (i.e. only the ' ...
                       'periodic components are removed).' ...
                       'The RETROICOR method is also included for comparison.' ...
                       '(see Glover et al. (2000) for details). Note that ', ...
                       'most parameters won''t have any effect, when using RETROICOR.'};
  
% Wether or not to save visualization of the IMM result to disk
  visual            = cfg_menu;
  visual.tag        = 'visual';
  visual.name       = 'IMM visualization output';
  visual.val        = {0};
  visual.labels     = {'0: No figure'
                       '1: Save figure'};
  visual.values     = {0 1};
  visual.help       = {'Choose ''Save figure'' to save a visualization of the IMM ' ...
                       'estimation results. The figure is saved under the current' ...
                       'path as a PNG file.'};
                   
% Method main branch setup
  drifter           = cfg_exbranch;
  drifter.tag       = 'drifter';
  drifter.name      = 'DRIFTER';
  drifter.val       = {epidata prefix mode visual refdata};
  drifter.help      = {'A SPM toolbox for removing periodic noise in fMRI data', ...
                       'by Arno Solin and Simo Sarkka', ...
                       '', ...
                       ['This SPM toolbox is an implementation of the DRIFTER ', ...
                       'algorithm [1], which is a Bayesian method for physiological ', ...
                       'noise modeling and removal allowing accurate dynamical ', ...
                       'tracking of the variations in the cardiac and respiratory ', ...
                       'frequencies by using Interacting Multiple Models (IMM), ', ...
                       'Kalman Filter (KF) and Rauch-Tung-Striebel (RTS) smoother ', ...
                       'algorithms.'], ...
                       '', ...
                       ['The frequency trajectories can be either estimated ', ...
                       'from reference signals (e.g., by using pulse meters ', ...
                       'and respiratory belts), or if the time resolution allows, ', ...
                       'directly from the fMRI signal. The estimated frequency ', ...
                       'trajectory is used for accurate model based separation ', ...
                       'of the fMRI signal into activation, physiological ', ...
                       'noise and white noise components using Kalman filter ', ...
                       'and RTS smoother algorithms. This separation is done ', ...
                       'for each voxel in the image separately. The basic idea ', ...
                       'of the method is to build a stochastic model for each ', ...
                       'component of the signal: the cleaned BOLD signal ', ...
                       '(including haemodynamical effects) is a  relatively ', ...
                       'slowly varying signal, cardiac and respiration induced ', ...
                       'noise components are stochastic resonators with multiple ', ...
                       'harmonics, and the rest of the signal is white noise.'], ...
                       '', ...
                       ['To get started using the DRIFTER toolbox, the user is ', ...
                       'required to specify at least the following parameters ', ...
                       'for the data section:'], ...
                       '  - EPI files : the fMRI data files', ...
                       '  - TR : Repetition time (in ms)', ...
                       'And for each reference signal:', ...
                       '  - Reference signal data', ...
                       '  - Sampling interval', ...
                       '  - List of possible frequencies', ...
                       '', ...
                       ['The physiological noise components are usually cardiac- ', ...
                       'and respiration-induced signals. However, the toolbox ', ...
                       'is not limited to these two in any way. Furthermore, ', ...
                       'the toolbox interface lets the user modify virtually ', ...
                       'every parameter in the estimation process. The default ', ...
                       'setup is not optimized for any particular purpose, ', ...
                       'and the user is encouraged to test the effects and ', ...
                       'abilities of the toolbox for the particular data at ', ...
                       'hand.'], ...
                       '', ...
                       'For further details on the toolbox, refer to the article [1].', ...
                       '', ...
                       'References:', ...
                       ['  [1] Sarkka, S., Solin, A., Nummenmaa, A., Vehtari, , ', ...
                       'A., Auranen, T., Vanni, S., and Lin, F.-H. (2012). ', ...
                       'Dynamical retrospective filtering of physiological ', ...
                       'noise in BOLD fMRI: DRIFTER. NeuroImage, 60:1517-1527.'], ...
                       '', ...
                       '', ...
                       'Copyright (C) 2011-2012  Arno Solin and Simo Sarkka', ...
                       '', ...                       
                       ['This program is free software: you can redistribute it and/or modify ', ...
                       'it under the terms of the GNU General Public License as published by ', ...
                       'the Free Software Foundation, either version 3 of the License, or ', ...
                       '(at your option) any later version.'], ...
                       '', ...
                       ['This program is distributed in the hope that it will be useful, ', ...
                       'but WITHOUT ANY WARRANTY; without even the implied warranty of ', ...
                       'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ', ...
                       'GNU General Public License for more details.'], ...
                       '', ...
                       ['You should have received a copy of the GNU General Public License ', ...
                       'along with this program.  If not, see <http://www.gnu.org/licenses/>.']};
  drifter.prog      = @spm_local_drifter;


%% Define batch run method
  
function spm_local_drifter(job)

  if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','DRIFTER')); end
  spm_drifter(job);




%% For documentation description writing only (not affecting SPM)

  %{

  % Set values
  cfg_entry = struct;
  cfg_branch = struct;
  cfg_files = struct;
  cfg_entry = struct;
  cfg_repeat = struct;
  cfg_menu = struct;
  cfg_exbranch = struct;

  % Get all variables
  vars = who;
  
  str = '';
  
  for i=1:length(vars)
      
      % Get variable value
      foo = evalin('base',vars{i});
      
      % Check if valid and output
      if isstruct(foo) && isfield(foo,'tag') && ~strcmpi(vars{i},'epidata')
      
      
        % Fix help value
        if iscell(foo.help)
          helpstr = cell2mat(foo.help);
        else
          helpstr = foo.help; 
        end
      
        % Fix default value
        if isfield(foo,'val')
         if iscell(foo.val)
          val = mat2str(cell2mat(foo.val),3);
         else
          val = mat2str(foo.val,3); 
         end
        else
          val = '';  
        end
      
        if isfield(foo,'tag')
            str = sprintf('%s\n {\\bf %s} \\vfill\n \\texttt{%s} \\vfill\n \\emph{%s} \n & %s  \\\\ \n \\hline', ...
                str,foo.name,foo.tag,val,helpstr);
        end
      end
      
  end

  str
  %}
