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


if numel(I.MOD.RFMs)~=2 && numel(I.MOD.RFMs.Effect)~=1
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

n = size(I.MOD.EffParts{1},2)+1;
df1 = n-1;
df2 = (size(Y,1)/n)-1;
df3 = size(Y,1)-df1-df2-1;


X = I.X;

er = eye(size(X,1))-(X*pinv(X));
MSE = LoopEstimate(Y,1,er)./df3;

[SSm r1 r2] = makeSSmat(I.MOD.RFMs(2).Effect(1).tx1, I.MOD.RFMs(2).Effect(1).tx2);
% df = ResidualDFs(I.MOD.RFMs(2).Effect(1).tx1)-ResidualDFs(I.MOD.RFMs(2).Effect(1).tx2);
MSF = LoopEstimate(Y,1,SSm)./df1;

[SSm r1 r2] = makeSSmat([ones(size(Y,1),1) I.MOD.EffParts{1} I.MOD.ErrParts{1}], [ones(size(Y,1),1) I.MOD.EffParts{1}]);
df = ResidualDFs([ones(size(Y,1),1) I.MOD.EffParts{1} I.MOD.ErrParts{1}])-ResidualDFs([ones(size(Y,1),1) I.MOD.EffParts{1}]);
MSS = LoopEstimate(Y,1,SSm)./df2;

ICC21 = (MSS-MSE) ./ (MSS + (n-1)*MSE + (n*(MSF-MSE)/(df2+1)));

ICC21(1) = 1;
ICC21(2) = -1;

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


