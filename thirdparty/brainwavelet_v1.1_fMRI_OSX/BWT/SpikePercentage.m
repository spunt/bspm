function SP = SpikePercentage(mmcs,m)
%                  
% FUNCTION:     SpikePercentage -- Compute SP from max/min chain matrix.
%                
% USAGE:        SP = SpikePercentage(mmcs,m)
%
% Inputs:       mmcs    -- Binary matrix of max/min chains.
%                          Dimensions NX x NJ x NTS.
%               m       -- Original MODWT transform of time series matrix.
%                          Dimensions NX x NJ x Nts.
%
% Output:       sp      -- vector containing values for the Spike 
%                          Percentage at each time point.
%
% EXAMPLE:      sp = SpikePercentage(maxchains,3Dmat)
%
% AUTHOR:       Ameera X Patel
% CREATED:      10-01-2014
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     4
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: SpikePercentage.m 4 30-01-2014 BWTv1.1 axpatel


%% check inputs

fname=mfilename;
if nargin<1
    help(fname); return
end

err=struct();
err(1).inp=chkninput(nargin,[2,2],nargout,[0,1],fname);
err(2).inp=chkintype(mmcs,'double',fname);
err(3).inp=chkintype(m,'double',fname);

if sum(cat(1,err.inp))>=1
    SP=[]; return
end

%% SP calculation

mmcs(:,:,sum(sum(m,2),1)==0)=[];
nnzts=size(mmcs,3);
s1=squeeze(mmcs(:,1,:));
SP=(sum(s1,2)/nnzts)*100; 

end