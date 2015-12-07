function Mem = FindFreeRAM()
%
% FUNCTION:     FindFreeRAM -- Computes maximum available memory.
%                            
% USAGE:        Mem = FindFreeRAM()
%
% Output:       Mem         -- Outputs free memory in Gigabytes
%
% EXAMPLE:      Memory = FindFreeRAM
%
% AUTHOR:       Ameera X Patel
% CREATED:      22-01-2014
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     6
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: FindFreeRAM.m 6 10-02-2014 BWTv1.1 axpatel

%% check inputs

fname=mfilename;

if chkninput(nargin,[0,0],nargout,[0,1],fname)>=1
    Mem=[]; return
end

%% detect system and find free memory

if isunix==1 && ismac==1
    [E,Num]=system('top -l 1 | grep -e PhysMem | awk ''{print $10}'' ');
    unit=Num(isstrprop(Num,'alpha'));
    Mem=str2double(Num(isstrprop(Num,'digit')));

    if strcmpi(unit,'k');
        Mem=Mem/(2^20);
    elseif strcmpi(unit,'m');
        Mem=Mem/(2^10);
    elseif strcmpi(unit,'b');
        Mem=Mem/(2^30);
    end
    
else
	cprintf('_err','*** BrainWavelet Error ***\n')
	cprintf('err','You do not appear to be running Mac OS. This is an \n')
	cprintf('err','OSX Mavericks release. Please download the correct \n')
	cprintf('err','version from www.brainwavelet.org.\n\n')
	Mem=[];
	end

end