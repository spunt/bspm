function Lj = FiltWidth(wavelet,J)
%
% FUNCTION:     FiltWidth -- Computes wavelet filter width at each scale.
%                            
% USAGE:        Lj = FiltWidth(wavelet,J)
%
% Inputs:       wavelet   -- Wavelet for which to calculate width.
%                            Input must be a string.
%               J         -- Number of scales
%
% Output:       Lj        -- Filter width at each of 1:J scales.
%
% EXAMPLE:      Lj = FiltWidth('d4',7)
%
% AUTHOR:       Ameera X Patel
% CREATED:      11-01-2013
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     4
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: FiltWidth.m 4 30-01-2014 BWTv1.1 axpatel


%% check inputs

fname=mfilename;
if nargin<1
    help(fname); return
end

wavelets={'d4','d6','d8','d10','d12','d14','d16','d18','d20',...
      'la8','la10','la12','la14','la16','la18','la20',...
      'bl14','bl18','bl20','c6','c12','c18','c24','haar'};

err=struct();
err(1).inp=chkninput(nargin,[2,2],nargout,[0,1],fname);
err(2).inp=chkintype(J,'numeric',fname);
err(3).inp=chkintype(wavelet,'char',fname,wavelets);

if sum(cat(1,err.inp))>=1
    Lj=[]; return;
end

%% find filter details

wf=modwt_filter(wavelet);
L=wf.L;

%% compute filter widths

j=(1:J);
Lj=((2.^j-1)*(L-1))+1;

end