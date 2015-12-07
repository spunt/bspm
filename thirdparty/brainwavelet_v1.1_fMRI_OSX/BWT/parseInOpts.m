function optsOut = parseInOpts(opts,inputs)
%                  
% FUNCTION:     parseInOpts -- Parses varargin options to function and
%                              edits the defaults (saved as a struct).
%                
% USAGE:        optsOut = parseInOpts(opts,inputs)
%
% Inputs:       opts        -- Structure of valid options and defaults.
%               inputs      -- String of varargin inputs.
%
% Output:       optsOut     -- Opts structure with new values.
%
% EXAMPLE:      optsOut = parseInOpts(OptsDefault,varargin)
%
% AUTHOR:       Ameera X Patel
% CREATED:      10-01-2014
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     5
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: parseInOpts.m 5 30-01-2014 BWTv1.1 axpatel


%% check inputs

fname=mfilename;
if nargin<1
    help(fname); return
end

err=struct();
err(1).inp=chkninput(nargin,[2,2],nargout,[0,1],fname)>=1;
err(3).inp=chkintype(opts,'struct',fname);
err(3).inp=chkintype(inputs,'cell',fname);

if sum(cat(1,err.inp))>=1
    optsOut=[]; return
end

%% check fieldnames and edit defaults

optNames=fieldnames(opts);

if mod(length(inputs),2)>0
     cprintf('_err','*** BrainWavelet Error ***\n')
     cprintf('err','Please ensure you have input a value for each\n');
     cprintf('err','field specified.\n\n');
     optsOut=[]; return;
end
      
for inpPair=reshape(inputs,2,[])
   inpName=inpPair{1};

   if sum(strcmpi(inpName,optNames))~=0;
      opts.(inpName)=inpPair{2};
   else
      cprintf('_err','*** BrainWavelet Error ***\n')
      cprintf('err','%s is an invalid parameter name.\n',inpName)
      optsOut=[]; return;
   end
end

optsOut=opts;

end