function MakeHubMap(fn,thresh,reslice,greyThresh)
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
if nargin < 2;
    thresh = .2;
end

if nargin < 3
    reslice = [];
end

if nargin < 4
    greyThresh = [];
end

if ~isempty(reslice)
    if ischar(fn); fn = {fn}; end
    m = [];
    for ii = 1:numel(fn)

        hh = spm_vol(fn{ii});

        [mm mat] = SliceAndDice3(hh,[],reslice,[],[1 0],[]);
        mm = squeeze(mm);
        ss = size(mm);
        
        h = hh(1);
        h.mat = mat;
        h.dim = ss(1:3);

        mm = reshapeWholeBrain(size(mm),mm);
        m = [m; zscore(mm)];
    end
    clear mm;
else
    if ischar(fn); fn = {fn}; end
    h = spm_vol([fn{1} ',1']);
    
    m = [];
    for ii = 1:numel(fn);
        m = [m; zscore((FastRead(fn{ii})))];
    end
end

if ~isempty(greyThresh)
    hm = spm_vol(which('grey.nii'));
    msk = resizeVol2(hm,h,0);
    ind = find(msk>greyThresh);
else
    ind = find(std(m)>0);
end
%%
disp(['Computing matrix across ' num2str(numel(ind)) ' voxels']);

rr = zeros(1,numel(ind));
rr2 = zeros(1,numel(ind));

disp('Creating the correlation matrix');
cols = numel(ind);
rows = size(m,1);
blocks=300;

CC = (1:blocks:cols)'; CC(:,2) = [CC(1:end-1)+blocks-1; cols];
persisText('0% Complete');
for ii=1:size(CC,1);
    %disp([ii size(CC,1)]);
    persisText([sprintf('%2.1f', (ii/size(CC,1))*100) ' % Complete']);
    cor = single((m(:,ind(CC(ii,1):CC(ii,2)))' * m) / (rows-1));
    cor(find(cor==0))=NaN;
    rr(CC(ii,1):CC(ii,2)) =  nanmean(cor'.^2);
    rr2(CC(ii,1):CC(ii,2)) = sum(cor'>thresh);
end
persisText();

mm = nan(h.dim);
mm(ind) = rr;


[fp1 fp2 fp3] = fileparts(fn{end});
stem = [fp2 fp3];

h1 = h(1);
h1.dt = [16 0];
h1.fname = ['Hubs_Weighted_' stem];
spm_write_vol(h1,mm);

mm = nan(h.dim);
mm(ind) = (rr-nanmean(rr))/nanstd(rr);

h1 = h(1);
h1.dt = [16 0];
h1.fname = ['Hubs_WeightedZ_' stem];
spm_write_vol(h1,mm);

nan(h.dim);
mm(ind) = rr2;

h1 = h(1);
h1.dt = [16 0];
h1.fname = ['Hubs_Degree_' stem];
spm_write_vol(h1,mm);

nan(h.dim);
mm(ind) = zscore(rr2);

h1 = h(1);
h1.dt = [16 0];
h1.fname = ['Hubs_DegreeZ_' stem];
spm_write_vol(h1,mm);

