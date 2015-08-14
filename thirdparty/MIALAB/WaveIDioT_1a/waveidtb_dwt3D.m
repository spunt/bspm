function [app,LLH, LHL, LHH, HLL, HLH, HHL, HHH]  = waveidtb_dwt3D(x, J, af)

% 3-D Stationary Discrete Wavelet Transform (Un-decimated)
% Author: Siddharth Khullar
% Date Modified: 07/22/2010

% Initialize 3-D subbands for speed.
LLH = cell(1,J);LHL = LLH; HLL = LLH; HLH = LLH; HHL = LLH; HHH = LLH;

for k = 1:J
    
    [x , LLH{k}, LHL{k}, LHH{k}, HLL{k}, HLH{k}, HHL{k}, HHH{k}] = afb3D(x, af, af, af);
    
        x = circshift(x,[-ceil(length(af(:,1))/8),-ceil(length(af(:,1))/8),0]);
    app{k} = x;
end

end

%% Function defined to perform the 3D wavelet analysis
function [LLL, LLH, LHL, LHH, HLL, HLH, HHL, HHH] = afb3D(x, af1, af2, af3)

% 3D Analysis Filter Bank
% USAGE:
%    [lo, LLH, LHL, LHH, HLL, HLH, HHL, HHH] = afb3D(x, af1, af2, af3);

if nargin < 3
    af2 = af1;
    af3 = af1;
end

% filter along dimension 1
[L, H] = afb3D_A(x, af1, 1);

% filter along dimension 2
[LL LH] = afb3D_A(L, af2, 2);
[HL HH] = afb3D_A(H, af2, 2);

% filter along dimension 3
[LLL LLH] = afb3D_A(LL, af3, 3);
[LHL LHH] = afb3D_A(LH, af3, 3);
[HLL HLH] = afb3D_A(HL, af3, 3);
[HHL HHH] = afb3D_A(HH, af3, 3);

end

%% Function defined to call the convolution and decide directions
function [lo, hi] = afb3D_A(x, af, d)

% 3D Analysis Filter Bank
% (along one dimension only)

lpf = af(:, 1);     % lowpass filter
hpf = af(:, 2);     % highpass filter

% Check the direction for filtering (1 = x, 2 = y, 3 = z)
if d==1
    p = [2 1 3];
elseif d==2
    p = [1 2 3];
elseif d==3;
    p = [3 2 1];
end

% permute the volume according to p (Filtering direction)
x = permute(x, p);
[N1, N2, N3] = size(x);
lo = zeros(N1, N2, N3);
hi = zeros(N1, N2, N3);

% Low pass
for k = 1:N3
    if d==3
        lo(:, :, k) = sconv2SWT(squeeze(x(:, :, k)), lpf, 'col');
    else
        lo(:, :, k) = sconv2SWT(squeeze(x(:, :, k)), lpf, 'row');
    end;
    
end

% High pass
for k = 1:N3
    if d==3
        hi(:, :, k) = sconv2SWT(squeeze(x(:, :, k)), hpf, 'col');
    else
        hi(:, :, k) = sconv2SWT(squeeze(x(:, :, k)), hpf, 'row');
    end;
end

% permute dimensions of xx (inverse permutation)
lo = ipermute(lo, p);
hi = ipermute(hi, p);
end

%% Function created to perform Convolution with Wavelet functions without
 % perforimg any decimation and also compensating for the shift.
function [y]=sconv2SWT(x,h,direction)
% *************************************************************************
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
x1 = x;
if strcmp(direction, 'row')
    
    x  = padarray(x , [0 n] , 0 , 'both');
    x=imfilter(x,h','replicate','conv');
    y=x(:,(m/2)+1:s(2)+m/2);
    
elseif strcmp(direction, 'col')
    
    x  = padarray(x , [n 0] , 0 , 'both');
    x=imfilter(x,h,'replicate','conv');
    y=x((m/2)+1:s(1)+m/2 , :);
end

end