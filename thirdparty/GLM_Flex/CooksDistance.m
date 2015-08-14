function CD = CooksDistance(y,xx)
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%%
%%% Copyright (C) 2012,  Aaron P. Schultz
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

df1 = ResidualDFs(xx);
df2 = size(y,1)-df1;

%%% Run a basic GLM
pv = pinv(xx);
beta  = pv*y;
pred = xx*beta;
res   = y-pred;
ResSS = sum(res.^2);
MSE = ResSS./df2;
%%% Compute Cook's Distance
e2 = res.^2;
p = df1;
bc = pinv(xx'*xx);
% bc(abs(bc)<eps*(size(xx,1)))=0;
tmp = diag(xx*bc*xx');
CD = (e2./(p*repmat(MSE,size(e2,1),1))) .* (repmat(tmp./((1-tmp).^2),1,size(e2,2)));