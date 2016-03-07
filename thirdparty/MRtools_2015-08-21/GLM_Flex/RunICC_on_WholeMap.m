function RunICC_on_WholeMap
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

%%%  Runs ICC(2,1) from Shrout & Fleiss
load I;

if ~isempty(I.F.IN.Between) || numel(I.F.IN.Within)>1 || ~isempty(I.F.IN.Covar)
   error('This script can only be run on vanilla one-way repeated measures ANOVA.  No covariates, between subject terms, or multiple repeated factors are allowed.')
end

if exist('FinalDataSet.nii','file')
    h = spm_vol('FinalDataSet.nii');
else
    h = spm_vol(char(I.Scans));
end
FullIndex = 1:prod(I.v.dim);

mh = spm_vol('NN.nii');
msk = resizeVol(mh,I.v);
mskInd = find(msk==max(msk(:)));
mskOut = find(msk~=max(msk(:)));

[x y z] = ind2sub(I.v.dim,mskInd);
Y = zeros(numel(h),numel(x));
for ii = 1:numel(h);
    Y(ii,:) = spm_sample_vol(h(ii),x,y,z,0);
end

n = I.F.IN.Within;
df1 = n-1;
df2 = (size(Y,1)/n)-1;
df3 = size(Y,1)-df1-df2-1;

F = I.F;

X = I.X;
X(isnan(X))=0;

beta = pinv(X)*Y;

[r c cc rr] = MakeContrastMatrix('a',F.XX(:,(n+1):end),[df2+1],F.XX,X);
MSS = LoopEstimate(beta,X,r)/df2;
[r c cc rr] = MakeContrastMatrix('a',F.XX(:,1:n),[n],F.XX,X);
MSF = LoopEstimate(beta,X,r)/df1;

% [r c cc rr] = MakeContrastMatrix2(F,cols,levs,Xind)

er = eye(size(X,1))-(X*pinv(X));
MSE = LoopEstimate(Y,1,er)./df3;

ICC21 = (MSS-MSE) ./ (MSS + (n-1)*MSE + (n*(MSF-MSE)/(df2+1)));

% save tmp.mat MSS MSE MSF n df1 df2 df3

ICC21(1) = 1;

[m,h] = openIMG('NN.nii');
m(:) = 0;
h.dt = [64 0];
h.fname = 'ICC21.nii';
m(mskInd) = ICC21;
spm_write_vol(h,m);

[m,h] = openIMG('NN.nii');
m(:) = 0;
h.dt = [64 0];
h.fname = 'BSV.nii';
m(mskInd) = sqrt(MSS);
spm_write_vol(h,m);

[m,h] = openIMG('NN.nii');
m(:) = 0;
h.dt = [64 0];
h.fname = 'WSV.nii';
m(mskInd) = sqrt(MSF);
spm_write_vol(h,m);

[m,h] = openIMG('NN.nii');
m(:) = 0;
h.dt = [64 0];
h.fname = 'Mean.nii';
m(mskInd) = nanmean(Y);
spm_write_vol(h,m);


