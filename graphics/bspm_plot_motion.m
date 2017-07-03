function bspm_plot_motion(rpfile)
% BOB_BADSCAN Identifies scans and saves a nuisance regressor file
% 
% USAGE: bspm_plot_motion(rpfile)
%
% ARGUMENTS
%   rpfile = realignment parameter file (cell or char)
%

% ---------------------- Copyright (C) 2014 ----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1
    rpfile = files('rp*txt');
    if isempty(rpfile), mfile_showhelp; return; end
end
if iscell(rpfile), rpfile = char(rpfile); end

% load in the rp file 
rp = load(rpfile);
rp(:,4:6) = rp(:,4:6)*57.3;
rpdiff = diff(rp);

% raw motion
figure('Color','white');
subplot(2,1,1);
plot(rp, 'LineWidth', 1);
set(get(gca,'YLabel'), 'String', 'Movement (mm/degrees)', 'FontSize', 13);
set(get(gca,'XLabel'), 'String', 'Scan', 'FontSize', 13);
legendh = legend({'x' 'y' 'z' 'pitch' 'roll' 'yaw'});
set(legendh, 'Orientation', 'vertical', 'Location', 'bestoutside', 'FontSize', 12);
subplot(2,1,2);
plot(diff(rp), 'LineWidth', 1);
set(get(gca,'YLabel'), 'String', 'Displacement (mm/degrees)', 'FontSize', 13);
set(get(gca,'XLabel'), 'String', 'Scan', 'FontSize', 13);
legendh = legend({'x' 'y' 'z' 'pitch' 'roll' 'yaw'});
set(legendh, 'Orientation', 'vertical', 'Location', 'bestoutside', 'FontSize', 12);
 
 
 
 
