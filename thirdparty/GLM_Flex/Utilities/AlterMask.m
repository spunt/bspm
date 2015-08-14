function AlterMask(fn,opt,criteria)
% thresh = 26;
% %% criteria = 0:26
% cd /autofs/space/schopenhauer_003/users/BigTemplate;
% fn = 'test.nii';
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

v = spm_vol(fn);
[m xyz] = spm_read_vols(v);
m = (m~=0).*(~isnan(m)).*(isfinite(m));

loc = round([xyz' ones(size(xyz,2),1)]*inv(v.mat'));

X = 2:size(m,1)-1; 
Y = 2:size(m,2)-1;
Z = 2:size(m,3)-1;
xyz = zeros(numel(X)*numel(Y)*numel(Z), 3);
count = 0;
for ii = 1:length(Z)
    for jj = 1:length(Y)
        for kk = 1:length(X)
            count = count+1;
            xyz(count,1:3) = [X(kk) Y(jj) Z(ii)];
        end
    end
end

V = sub2ind(size(m),xyz(:,1),xyz(:,2),xyz(:,3));
%%% In Plane 
comp = [m(V-1) m(V+1) m(V-size(m,1)) m(V-size(m,1)+1) m(V-size(m,1)-1) m(V+size(m,1)) m(V+size(m,1)+1) m(V+size(m,1)-1)];
%%% In back of plane.
B = (V-(size(m,1)*size(m,2)));
comp = [comp m(B) m(B-1) m(B+1) m(B-size(m,1)) m(B-size(m,1)+1) m(B-size(m,1)-1) m(B+size(m,1)) m(B+size(m,1)+1) m(B+size(m,1)-1)];
%%% In front of plane.
F = (V+(size(m,1)*size(m,2)));
comp = [comp m(F) m(F-1) m(F+1) m(F-size(m,1)) m(F-size(m,1)+1) m(F-size(m,1)-1) m(F+size(m,1)) m(F+size(m,1)+1) m(F+size(m,1)-1)];

% %%% In Plane 
% indices = [(V-1) (V+1) (V-size(m,1)) (V-size(m,1)+1) (V-size(m,1)-1) (V+size(m,1)) (V+size(m,1)+1) (V+size(m,1)-1)];
% %%% In back of plane.
% B = (V-(size(m,1)*size(m,2)));
% indices = [indices (B) (B-1) (B+1) (B-size(m,1)) (B-size(m,1)+1) (B-size(m,1)-1) (B+size(m,1)) (B+size(m,1)+1) (B+size(m,1)-1)];
% %%% In front of plane.
% F = (V+(size(m,1)*size(m,2)));
% indices = [indices (F) (F-1) (F+1) (F-size(m,1)) (F-size(m,1)+1) (F-size(m,1)-1) (F+size(m,1)) (F+size(m,1)+1) (F+size(m,1)-1)];

if opt == -1
    a = find(m(V)==1);
    counts = sum(comp(a,:),2);
    ind = find(counts<criteria);
    
    m(V(a(ind)))=0;
    
    v.fname = [v.fname(1:end-4) '_deflated.nii'];
    spm_write_vol(v,m);
end

if opt == 1
    a = find(m(V)==0);
    counts = sum(comp(a,:),2);
    ind = find(counts>criteria);
    
    m(V(a(ind)))=1;
    
    v.fname = [v.fname(1:end-4) '_inflated.nii'];
    spm_write_vol(v,m);
end
