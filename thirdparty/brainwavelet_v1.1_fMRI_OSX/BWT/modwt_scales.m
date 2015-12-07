function nscales = modwt_scales(N,method,wavelet)
%
% FUNCTION:     modwt_scales -- Calculates maximum number of wavelet
%                               scales computable for any given time
%                               series length.
%                            
% USAGE:        nscales = modwt_scales(N)
%
% Inputs:       N            -- Length of time series / number of time
%                               points.
%               method       -- Method for estimating number of scales.
%                               Options are 'conservative' or 'liberal'.
%                               Input must be a string.
%                               [Default='liberal'].
%                            -- Wavelet used for analysis. This is used
%                               to determine filter length for conservative
%                               estimation of number of scales. Input must
%                               be a string.
%               
% Output:       nscales      -- Number of scales
%
% EXAMPLE:      nscales = modwt_scales(250)
%
% AUTHOR:       Ameera X Patel
% CREATED:      26-12-2013
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     5
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: modwt_scales.m 5 05-02-2014 BWTv1.1 axpatel


%% check inputs

fname=mfilename;
if nargin<1
    help(fname); return
end

wavelets={'d4','d6','d8','d10','d12','d14','d16','d18','d20',...
      'la8','la10','la12','la14','la16','la18','la20',...
      'bl14','bl18','bl20','c6','c12','c18','c24','haar'};

method_opts={'conservative','liberal','extreme'};
err=struct();

if exist('method','var')
    if strcmpi(method,'conservative') && ~exist('wavelet','var')
        cprintf('_err','*** BrainWavelet Error ***\n')
        cprintf('err','Error using %s - please specify a wavelet. \n',...
            fname);
        nscales=[];
        return
    end 
    err(1).inp=chkninput(nargin,[2,3],nargout,[0,1],fname);
    err(3).inp=chkintype(method,'char',fname,method_opts);
    if exist('wavelet','var')
        err(4).inp=chkintype(wavelet,'char',fname,wavelets);
    else
        err(4).inp=0;
    end
else
    err(1).inp=chkninput(nargin,[1,1],nargout,[0,1],fname);
    err(3).inp=0; err(4).inp=0;
    method='liberal';
end

err(2).inp=chkintype(N,'numeric',fname);

if sum(cat(1,err.inp))>=1
    nscales=[]; return
end


%% calculate max scales

nscales=[];

if strcmpi(method,'liberal')
    nscales= nextpow2(N)-1;
elseif strcmpi(method,'extreme')
    nscales= nextpow2(N*1.5)-1;
elseif strcmpi(method,'conservative')
    wtf=modwt_filter(wavelet);
    L=wtf.L;
    nscales=floor(log2(N/(L-1))+1);
end

end