%_______________________________________________________________________
%
% Wrapper for console progress bar
%_______________________________________________________________________
%
% This file is part of the Variability Toolbox for SPM
% Published by the Lifespan Neural Dynamics Group
% Provided under the GNU General Public License
% See LICENSE for terms and conditions

classdef shared_progbar

  properties
    prog_bar;
    maximum;
  end

  methods

    function obj = shared_progbar(maximum)
      prog_bar = shared_proglib();
      prog_bar.setPercentVisible(0)
      prog_bar.setMaximum(maximum);
      prog_bar.setLength(30);
      prog_bar.start();
      obj.prog_bar = prog_bar;
      obj.maximum = maximum;
    end

    function update(obj, value)
      text = sprintf('[%i/%i]', value, obj.maximum);
      progbar = obj.prog_bar;
      progbar.setValue(value);
      progbar.setText(text);
      fprintf('\033[K')
    end

    function stop(obj)
      progbar = obj.prog_bar;
      progbar.stop();
      fprintf('\n')
    end

  end
end
