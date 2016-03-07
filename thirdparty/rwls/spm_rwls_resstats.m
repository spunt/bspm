function varargout=spm_rwls_resstats(SPM,subset,movparam_name,startval)
% function varargout=spm_rwls_resstats(SPM,subset,movparam_name)
%
% It attempts to find the movemen parameters for each scan
% INPUT:
%      SPM - SPM structure estimated with the RWLS toolbox
%      subset - Subset of scans to be plotted (defaults is all)
%      moveparam_name - Cell array of bames of movement parameter files. If not given the
%      routine will attempt to find the file from the image information in
%      the SPM structure.
%
% Joern Diedrichsen (j.diedrichsen@ucl.ac.uk) v.3.0
%

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Plot residual time series',0);

% Starting value of the movement parameter file: Use if movement parameters
% include more scans than the GLM 
if (nargin<4 || isempty(startval))
    startval=1; 
end; 


% determine the subset: force subset to be in range of 1:T
%-----------------------------------------------------------------------
if (~isfield(SPM,'ResStats'))
    error('SPM needs to be estmated with RWLS toolbox');
end;

T=size(SPM.ResStats.s,1);

if (nargin<2 || isempty(subset))
    subset=[1:T];
else
    subset=subset(subset>=1 & subset<=T);
end;

%-----------------------------------------------------------------------
% Try to find the movement parameter files for each of the blocks
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);
figure(Fgraph);
MOV=[];

if (nargin<3 || isempty(movparam_name) || isempty(movparam_name{1}))
    start=[0 cumsum(SPM.nscan)]+1;
    for i=1:length(SPM.nscan)
        [dir,filename]=spm_fileparts(SPM.xY.VY(start(i)).fname);
        movparam_name{i}=[dir filesep 'rp_' filename(2:end) '.txt'];
    end;
end;

%-----------------------------------------------------------------------
% Load movement parameter files
for i=1:length(movparam_name)
    if (~isempty(movparam_name{i}))
        try
            mov=dlmread(movparam_name{i});
            MOV=[MOV;mov(startval:end,:)];
        catch
            warning(['Movementparameter file ' movparam_name{i} ' not found.']);
        end;
    end;
end;
if (size(MOV,1)~=sum(SPM.nscan))
    warning('Number of scans in movement parameter file do not match information in SPM.mat');
    MOV=[];
end;


%-----------------------------------------------------------------------
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
