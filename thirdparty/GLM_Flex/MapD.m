function MapD
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

% For use with two-sample t-tests
load I;
i1 = find(I.X(:,1)==1);
i2 = find(I.X(:,2)==1);

h = spm_vol(char(I.Scans(i1)));
m1 = FastRead(h);
% [m1 h] = openIMG(char(I.Scans(i1)));
% m1 = reshapeWholeBrain(size(m1),m1);

m2 = FastRead(I.Scans(i2));
% m2 = openIMG(char(I.Scans(i2)));
% m2 = reshapeWholeBrain(size(m2),m2);

D = cohensD(m1,m2);

hh = h(1);
hh.fname = 'CohensD.nii';
V = zeros(hh.dim);
V(:) = D;

spm_write_vol(hh,V);
