function [x] = chkninput(nin,limin,nout,limout,fname)
%                  
% FUNCTION:     chkninput -- Error function that checks number of arguments 
%                            given to a function.
%                
% USAGE:        x = chkninput(nin,limin,nout,limout,fname)
%
% Inputs:       nin       -- Number of actual inputs.
%               limin     -- Limits of accepted number of inputs.
%               nout      -- Number of actual outputs.
%               limout    -- Limits of accepted number of outputs.
%               fname     -- Name of function. Input must be a string.
%
% Output:       x         -- Error counter.
%
% EXAMPLE:      errors = chkinput(nargin,[3 4],nargout,[0 1],mfilename)
%
% AUTHOR:       Ameera X Patel
% CREATED:      23-12-2013
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     6
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: chkninput.m 6 30-01-2014 BWTv1.1 axpatel


%% check input/output to target fn

if nargin<1
    help(mfilename); return;
end

x=0;

if nin<limin(1)
    cprintf('_err','*** BrainWavelet Error ***\n')
    cprintf('err','Error using %s: Not enough input arguments\n',fname);
    x=x+1;
elseif nin>limin(2)
    cprintf('_err','*** BrainWavelet Error ***\n')
    cprintf('err','Error using %s: Too many input arguments\n',fname);
    x=x+1;
end

if nout<limout(1)
    if x>=1
        cprintf('err','Error using %s: Not enough output arguments\n',fname);
    else
        cprintf('_err','*** BrainWavelet Error ***\n')
        cprintf('err','Error using %s: Not enough output arguments\n',fname);
    end
    x=x+1;
elseif nout>limout(2)
    if x>=1;
        cprintf('err','Error using %s: Too many output arguments\n',fname);
    else
        cprintf('_err','*** BrainWavelet Error ***\n')
        cprintf('err','Error using %s: Too many output arguments\n',fname);
    end
    x=x+1;
end

%  TODO: create recursive call of chkninput to remomve narg*chk dependency
end