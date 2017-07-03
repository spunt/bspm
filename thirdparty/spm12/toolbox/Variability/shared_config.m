%_______________________________________________________________________
%
% Configuration and helper functions
%_______________________________________________________________________
%
% Output
%
% cfg | toolbox configuration (struct)
%
%_______________________________________________________________________
%
% This file is part of the Variability Toolbox for SPM
% Published by the Lifespan Neural Dynamics Group
% Provided under the GNU General Public License
% See LICENSE for terms and conditions

function cfg = shared_config

	%%
	%% release informaton
	%%
  cfg.name = 'Variability Toolbox for SPM';
  cfg.version = '0.1';

	%%
	%% formatting
	%%
	cfg.pad = 20;
  cfg.format = '\fontsize{16}';

  %%
  %% modeling
  %%
  cfg.modeltype = {'boxcar'};
  cfg.metric = {'var', 'sd', 'mssd', 'sqrt_mssd'};

  %%
  %% is graphical interface available?
  %%
  if usejava('jvm') && feature('ShowFigureWindows');
    cfg.gui = true;
  else
    cfg.gui = false;
  end

  %%
  %% is the matlab desktop running?
  %%
  if cfg.gui && desktop('-inuse')
    cfg.desktop = true;
  else
    cfg.desktop = false;
  end
  
  %%
  %% is spm running?
  %%
  if cfg.gui && not(isempty(spm_figure('FindWin','Menu')))
    cfg.spm = true;
  else
    cfg.spm = false;
  end

end
