function D = getROI(fn, ss, mni)
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

V1 = spm_vol(fn);
ch = length(V1);
M1 = spm_read_vols(V1);
if length(V1)>1
    V1 = spm_vol([fn ',1']);
end
    
if isnumeric(mni)
    rad = ss/2;
    dd = [];
    for ii = 1:.5:rad
        [X,Y,Z] = sphere(50);
        X = ((X*ii)+mni(1,1));
        Y = ((Y*ii)+mni(1,2));
        Z = ((Z*ii)+mni(1,3));
        cc = [X(1:end)' Y(1:end)' Z(1:end)'];
        cc = [floor(cc); ceil(cc)];
        cc = unique(cc,'rows');
        dd = [dd;cc];
    end
    dd = [dd;mni(1,:)];
    dd = unique(dd,'rows');
    
    dd = [dd ones(size(dd,1),1)]*(inv(V1.mat)');
    dd = [ceil(dd(:,1:3)); floor(dd(:,1:3))];
    dd = unique(dd,'rows');
%     ind = square2vect3D(dd, V1.dim);
    ind = sub2ind(V1.dim,dd(:,1),dd(:,2),dd(:,3));
else
    [mm vv] = openIMG(mni);
    if sum(vv.dim==V1.dim)==3
        ind = find(mm>0);
    else
        nn = resizeVol(vv,V1);
        ind = find(nn>0);
    end
end


D.mean = nanmean(M1(ind));
D.allVals = M1(ind);
   
%%% Singular value decomposition of voxels in sphere
y = [M1(ind)'];
nNaN(1:2) = [numel(find(~isnan(y))) numel(find(isnan(y)))];
D.nNaN = nNaN;

y = y(find(~isnan(y)));
[m n]   = size(y);
[u s u] = svd(y*y');
s       = diag(s);
u       = u(:,1);
v       = y'*u/sqrt(s(1));
d       = sign(sum(v));
u       = u*d;
v       = v*d;
Y       = u*sqrt(s(1)/n);
D.svd = Y;

    
