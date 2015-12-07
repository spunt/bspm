function [NBatch,RAMout] = MemorySolver(NTS,NX,varargin)
%
% FUNCTION:     MemorySolver -- Computes number of batches in which to
%                               process data with WaveletDespike, in order
%                               to conserve memory usage.
%                            
% USAGE:        NBatch = MemorySolver(NVox,NTS,varargin)
%
% Inputs:       NTS          -- Number of time series to process.
%               NX           -- Number of time points in the time series.
%
%               Additional Input Options:
%               (These must be specified as MATLAB string-value pairs).
%
%               LimitRam     -- Maximum number of Gb of RAM to use.
%
% Output:       NBatch       -- Number of batches recommended for 
%                               processing the data.
%
% EXAMPLE:      Batches = MemorySolver(40000,300,'LimitRAM',3)
%
% AUTHOR:       Ameera X Patel
% CREATED:      22-01-2014
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     7
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: MemorySolver.m 7 10-02-2014 BWTv1.1 axpatel


%% check inputs

fname=mfilename;
if nargin<1
    help(fname); return
end

if chkninput(nargin,[2,4],nargout,[0,2],fname)>=1
   NBatch=[]; 
   if nargout==2; RAMout=[]; end
   return;
end    

DefaultOpts=struct('LimitRAM',0);
Opts=parseInOpts(DefaultOpts,varargin);

err=struct();
err(1).inp=chkintype(NTS,'numeric',fname);
err(2).inp=chkintype(NX,'numeric',fname);
err(3).inp=chkintype(Opts.LimitRAM,'numeric',fname);

if sum(cat(1,err.inp))>=1;
    NBatch=[];
    if nargout==2; RAMout=[]; end
    return;
end

FreeRAM=FindFreeRAM;
if isempty(FreeRAM);
	cprintf('_err','*** BrainWavelet Error ***\n');
    cprintf('err','Could not identify amount of free memory. Please enter\n');
    cprintf('err','this manually. Type ''help WaveletDespike'' for more\n');
    cprintf('err','information\n\n');
    NBatch=[]; 
    if nargout==2; RAMout=[]; end
    return;
end

if Opts.LimitRAM > FreeRAM
    cprintf('_[1,0.5,0]','*** BrainWavelet Warning ***\n');
    cprintf([1,0.5,0],'You do not have %s Gb of free RAM.\n',...
        num2str(Opts.LimitRAM));
    cprintf([1,0.5,0],'Using maximum available.\n\n');
    Opts.LimitRam=0;
    RAMout=FreeRAM;
else
    RAMout=Opts.LimitRAM;
end

%% compute number of batches required

if Opts.LimitRAM==0;
    DimMax=(FreeRAM+0.5+0.1564)*(2^30)/(8*94.6116);
else
    DimMax=(Opts.LimitRAM+0.5+0.1564)*(2^30)/(8*94.6116);
end

VoxMax=floor(DimMax/NX);
NBatch=NTS/VoxMax;

if NBatch<=1
    NBatch=0;
else
    NBatch=ceil(NBatch);
end

end