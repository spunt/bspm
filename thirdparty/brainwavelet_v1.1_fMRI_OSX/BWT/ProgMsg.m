function [] = ProgMsg(pos,Nout)
%                  
% FUNCTION:     ProMsg -- Shows progress messages.
%                
% USAGE:        ProgMsg(message)
%
% Inputs:       pos    -- ID of Progress message to be displayed.
%               Nout   -- Number of output args.
%
% AUTHOR:       Ameera X Patel
% CREATED:      10-01-2014
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     3
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: ProgMsg.m 3 30-01-2014 BWTv1.1 axpatel


%% check number of inputs

fname=mfilename;
if nargin<1
    help(fname); return
end

if chkninput(nargin,[2,2],nargout,[0,0],fname) >=1 || ...
        chkintype(Nout,'numeric',fname,{'2','3','4'}) >=1 ;
    return;
end

%% main msg struct

msg=struct();

msg.inp='1/%s  Checking inputs                  ...';
msg.modwt='2/%s  Computing scales and MODWT       ...';
msg.talc='3/%s  Temporally aligning coefficients ...';
msg.lmm='4/%s  Searching for local max/min      ...';
msg.nghbr='5/%s  Computing neighbourhood matrices ...';
msg.ch='6/%s  Searching for chains             ...';
msg.pdsh='7/%s  Applying phase-delay shift       ...';
msg.sn='8/%s  Finding signal and noise coefs   ...';
msg.imodwt='9/%s  Recomposing time series          ...';
msg.sp='10/%s Computing Spike Percentage       ...';
msg.dof='11/%s Computing spatial DOF map        ...';
msg.done='%c\n';

steps=num2str(Nout+7);

%% validate message fields

if  chkintype(pos,'char',fname,fieldnames(msg))>=1
    return;
end

%% print message

if strcmpi(pos,'done')
%    fprintf(' done\n');
    fprintf(' %c\n',char(216));
else
    fprintf(msg.(pos),steps);
end

end