function W = KendallW(data,opt)
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%% Copyright (C) 2014,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.

if nargin>1 %% Perform ranking
    [~,b] = sort(data);
    [~,d] = sort(b);
    data = d;
    clear b d
end

m = size(data,2);
n = size(data,1);

Ri = sum(data,2);
Rbar = .5*m*(n+1);
S = sum((Ri-Rbar).^2);

W = (12*S)/ (m^2 * (n^3 -n));