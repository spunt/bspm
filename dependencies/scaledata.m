function out = scaledata(in, minmax)
% SCALEDATA
%
% USAGE: out = scaledata(in, minmax)
%
% Example:
% a = [1 2 3 4 5];
% a_out = scaledata(a,0,1);
% 
% Output obtained: 
%            0    0.1111    0.2222    0.3333    0.4444
%       0.5556    0.6667    0.7778    0.8889    1.0000
%
% Program written by:
% Aniruddha Kembhavi, July 11, 2007
if nargin<2, minmax = [0 1]; end
if nargin<1, error('USAGE: out = scaledata(in, minmax)'); end
out = in - repmat(min(in), size(in, 1), 1); 
out = ((out./repmat(range(out), size(out,1), 1))*(minmax(2)-minmax(1))) + minmax(1); 
end