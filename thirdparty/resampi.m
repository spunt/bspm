function [xrs,trs] = resampi(x,fs,fsrs,method,antialiasing)

% [xrs,trs] = resampi(x,fs,fsrs,method)
%  
% Resampling to arbitrary sampling rate using linear interpolation or Fourier
% domain resizing. Linear interpolation may cause periodic attenuation of high 
% frequencies depending on where the interpolation points fall relative to
% the original data. Fourier domain interpolation preserves spectral content
% but may create high-frequency artifacts in the time domain. Use whichever
% most suits your requirements.
%
%   x - original signal as a column vector or matrix.
%   fs - original sampling frequency
%   fsrs - target sampling frequency
%   method - 'linear' (default) for linear interpolation, 'fft' for
%            fourier, otherwise methods accepted by interp1, including 
%             'spline' and 'pchip'. 

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C. Kovach 2012

if nargin < 4 || isempty(method)
    method = 'linear';
end

if nargin < 5
    antialiasing = true;
end

N = size(x,1);
torig = (1:N)./fs; % original time points

trs = 1/fsrs:1/fsrs:torig(end); %resampled time points

% Anti-aliasing low-pass filter 
if antialiasing && fsrs < fs 
    
    filt = [fs fsrs/2];

%     sband = .2;
    sband = 1/length(trs);

    if size(x,1) == 1
        x = x';
    end


    n = size(x,1);

    fnorm  = filt(2)./(filt(1)./2);
    %mx = mean(x);
    %x = x-mx(ones(size(x,1),1),:);
    %B = fir1(size(x,1)-1, fnorm);
    B = fir2(size(x,1)-1, [0, fnorm(1), fnorm(1)+sband, 1], [1 1 0 0] );

    F = abs(fft(cat(2,B,  zeros(1,size(x,1) ))))';

    x = real( ifft( fft(cat(1,x,zeros(size(x)))).*F(:,ones(size(x,2),1)) ) );

    x = x(1:n,:);
 
end

switch lower(method)
    
    case 'linear'
        xrs = interp1q(torig',x,trs');
    case 'fft'
        xrs = interpft(x,length(trs));
    otherwise
        xrs = interp1(torig',x,trs',method);
end

       





