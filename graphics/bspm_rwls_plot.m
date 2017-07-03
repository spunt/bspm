function bspm_rwls_plot(spmmat)
% BSPM_RWLS
%   USAGE: bspm_rwls_plot(spmmat)
%

% ---------------------- Copyright (C) 2014 ----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<1, spmmat = 'SPM.mat'; end
if ~isstruct(spmmat),
    if iscell(spmmat), spmmat = char(spmmat); end
    load(spmmat);
end

% xmatrix = SPM.xX.xKXs.X; % filtered and whitened design matrix
xmatrix = SPM.xX.X; % raw design matrix
xname = SPM.xX.name; % regressor names
nuisidx = find(cellregexp(xname, '.*\sR\d+$'));
MOV = xmatrix(:,nuisidx(1:6));
% MOV(:,4:6) = MOV(:,4:6)*57.2958;
subset = 1:size(xmatrix, 1);
startval=1;


% raw motion
figure('Color','white');
% subplot(2,1,1);
% plot(rp, 'LineWidth', 1);
% set(get(gca,'YLabel'), 'String', 'Movement (mm/degrees)', 'FontSize', 13);
% set(get(gca,'XLabel'), 'String', 'Scan', 'FontSize', 13);
% legendh = legend({'x' 'y' 'z' 'pitch' 'roll' 'yaw'});
% set(legendh, 'Orientation', 'vertical', 'Location', 'bestoutside', 'FontSize', 12);
% subplot(2,1,2);
% plot(diff(rp), 'LineWidth', 1);
% set(get(gca,'YLabel'), 'String', 'Displacement (mm/degrees)', 'FontSize', 13);
% set(get(gca,'XLabel'), 'String', 'Scan', 'FontSize', 13);
% legendh = legend({'x' 'y' 'z' 'pitch' 'roll' 'yaw'});
% set(legendh, 'Orientation', 'vertical', 'Location', 'bestoutside', 'FontSize', 12);

% Plot movementparams
start=[cumsum(SPM.nscan)+1];
start=start(start>subset(1) & start<subset(end));
if (~isempty(MOV))
    subplot(4,1,1);
    plot(subset,MOV(subset,1:3));
    legend({'x','y','z'});
    legend(gca,'boxoff');
    set(gca,'Box','off');
    ylabel('Translation [mm]');
    drawline(start);

    subplot(4,1,2);
    plot(subset,MOV(subset,4:6));
    legend({'Pitch','Roll','Yaw'});
    legend(gca,'boxoff');
    set(gca,'Box','off');
    ylabel('Rotation [deg]');
    drawline(start);
end;

max_sd=[];
if (strcmp(SPM.xVi.form,'wls'))
    subplot(4,1,3);
    sd_pre=sqrt(SPM.xVi.h);
    max_sd=max(sd_pre);
    plot(subset,sd_pre(subset),'k.');
    set(gca,'Box','off','YLim',[0 max_sd*1.1]);
    ylabel('Pre-WLS residual SD');
    drawline(start);
end;


subplot(4,1,4);
sd_post=sqrt(SPM.ResStats.ss./SPM.ResStats.n);
if (isempty(max_sd))
    max_sd=max(sd_post);
end;
plot(subset,sd_post(subset),'k.');
set(gca,'Box','off','YLim',[0 max_sd*1.1]);
ylabel('Final Residual SD');
drawline(start);

function drawline(x)
lim=get(gca,'YLim');
x1=[x;x];
y=[ones(1,length(x));ones(1,length(x))];
y(1,:)=y(1,:).*lim(1);
y(2,:)=y(2,:).*lim(2);
l=line(x1,y,'Color','k');
