function wavcs = calccshift(wavelet,J)
%
% FUNCTION:     calcshift -- Computes circular shift of wavelet
%                            coefficients.
%                            
% USAGE:        wavcs = calccshift(wavelet,J)
%
% Inputs:       wavelet   -- Wavelet for which to calculate circular shift.
%                            Input variable must be a string.
%               J         -- Number of scales.
%
% Output:       wavcs     -- Circular shift of coefficients for J scales.
%
% EXAMPLE:      wavcs = calccshift('d4',7)
%
% AUTHOR:       Ameera X Patel
% CREATED:      26-12-2013
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     4
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: calccshift.m 4 30-01-2014 BWTv1.1 axpatel


%% check input arguments

fname=mfilename;
if nargin<1
    help(fname); return
end

wavelets={'d4','d6','d8','d10','d12','d14','d16','d18','d20',...
      'la8','la10','la12','la14','la16','la18','la20',...
      'bl14','bl18','bl20','c6','c12','c18','c24','haar'};

err=struct();
err(1).inp=chkninput(nargin,[2,2],nargout,[0,1],fname);
err(2).inp=chkintype(wavelet,'char',fname,wavelets);
err(3).inp=chkintype(J,'numeric',fname);

if sum(cat(1,err.inp))>=1
	wavcs=[]; return;
end

clear fname err

%% compute circular shift vector

wavcs=zeros(1,J);

for ii=1:J
    wavcs(1,ii) = advance_wavelet_filter(wavelet,ii);
end

end