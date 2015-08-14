function y = idwt3D(w, J, sf)

% 3-D Stationary Discrete Wavelet Transform (Un-decimated)
% Author: Siddharth Khullar
% Date Modified: 07/22/2010

y = cell2mat(w{J}{8});
for k = J:-1:1
   y = sfb3D(y, w{k}, sf, sf, sf);
%    y = circshift(y,[-ceil(length(sf(:,1))/8),-ceil(length(sf(:,1))/4),0]);
end

% Routine for analysis

function y = sfb3D(lo, hi, sf1, sf2, sf3)

% 3D Synthesis Filter Bank
%
% USAGE:
%   y = sfb3D(lo, hi, sf1, sf2, sf3);
% INPUT:
%   lo, hi - lowpass subbands
%   sfi - synthesis filters for dimension i
% OUPUT:
%   y - output array
% See afb3D
%

if nargin < 4
   sf2 = sf1;
   sf3 = sf1;
end

LLL = lo;LLL(isnan(LLL)) = 0;
LLH = hi{1};LLH(isnan(LLH)) = 0;
LHL = hi{2};LHL(isnan(LHL)) = 0;
LHH = hi{3};LHH(isnan(LHH)) = 0;
HLL = hi{4};HLL(isnan(HLL)) = 0;
HLH = hi{5};HLH(isnan(HLH)) = 0;
HHL = hi{6};HHL(isnan(HHL)) = 0;
HHH = hi{7};HHH(isnan(HHH)) = 0;


% filter along dimension 3
LL = sfb3D_A(LLL, LLH, sf3, 3);
LH = sfb3D_A(LHL, LHH, sf3, 3);
HL = sfb3D_A(HLL, HLH, sf3, 3);
HH = sfb3D_A(HHL, HHH, sf3, 3);


% filter along dimension 2
L = sfb3D_A(LL, LH, sf2, 2);
H = sfb3D_A(HL, HH, sf2, 2);


% filter along dimension 1
y = sfb3D_A(L, H, sf1, 1);


function y = sfb3D_A(lo, hi, sf, d)

% 3D Synthesis Filter Bank
% (along single dimension only)
%
% y = sfb3D_A(lo, hi, sf, d);
% sf - synthesis filters
% d  - dimension of filtering
% see afb2D_A
[xl,yl,zl] = size(lo);
[xh,yh,zh] = size(hi);

xx = min(xl,xh);
yy = min(yl,yh);
zz = min(zl,zh);

lo = lo(1:xx,1:yy,1:zz);
hi = hi(1:xx,1:yy,1:zz);

lpf = sf(:, 1);     % lowpass filter
hpf = sf(:, 2);     % highpass filter

% permute dimensions of lo and hi so that dimension d is first.
if d==1
    p = [2 1 3];
elseif d==2
    p = [1 2 3];
elseif d==3;
    p = [3 2 1];
end;

lo = permute(lo, p);
hi = permute(hi, p);

[N1, N2, N3] = size(lo);
N = 2*N1;
L = length(sf);
y = zeros(N1, N2, N3);

for k = 1:N3
    if d==3
        y(:, :, k) = sconv2SWT(squeeze(lo(:, :, k)), lpf, 'col') + sconv2SWT(squeeze(hi(:, :, k)), hpf, 'col') ;
    else
        y(:, :, k) =  sconv2SWT(squeeze(lo(:, :, k)), lpf, 'row') + sconv2SWT(squeeze(hi(:, :, k)), hpf, 'row') ;
    end;
    
end

% permute dimensions of y (inverse permutation)
y = ipermute(y, p);


function [y]=sconv2SWT(x,h,direction)
% *************************************************************************
% % *************************************************************************
% INPUT--> x: Input signal
%          h: Filter response
%          direction: 'row' or 'col'
% OUTPUT--> y: Output signal
%**************************************************************************
s=size(x);
m=length(h);

if mod(m,2)==0
    n = (m)/2;
else
    n=(m-1)/2;
end;
y=x;

if strcmp(direction, 'row')
    x  = padarray(x , [0 n] , 0 , 'both');
    x=imfilter(x,h','replicate','conv');
    y=x(:,(m/2):s(2)+m/2-1);
    
elseif strcmp(direction, 'col')
    x  = padarray(x , [n 0] , 0 , 'both');
    x=imfilter(x,h,'replicate','conv');
    y=x((m/2):s(1)+m/2-1, :);
    
end