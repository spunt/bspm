function  TBR(fn,targets,controls,subset,tag,outdir)
% This function implements the TBR functional mapping procedure.  An SPM 
% installation is required.
%
% FN - A cell array of input files (expects 4D-Nifti files)
%
% TARGETS - a cell array of file paths to target template maps.
%
% CONTROLS - a cell array of file paths to template maps that should be controled for.
%
% SUBSET - an OPTIONAL index of time points to include. Inputs should be vectors nested within a cell array (e.g. {[1:50 60:100] [1:99]}.
%
% TAG - an optional postfix for output directory names.
%
% OUTDIR - an option input to specify the path to the output directory.
% 
% Written by Aaron P. Schultz - aschultz@martinos.org
%
% Copyright (C) 2013,  Aaron P. Schultz
%
% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

hm = pwd;

if nargin<3
    controls = [];
end

if nargin<5
    tag = '';
end

if nargin<6 || isempty(outdir);
    outdir = hm;
end

if ischar(fn);
    fn = cellstr(fn);
end

if nargin<4 || isempty(subset); 
    subset = {};
    for ii = 1:numel(fn)
        subset{ii} = [];
    end
elseif isnumeric(subset);
    subset = mat2cell(subset);
end

h = spm_vol([fn{1} ',1']);

II = []; M = [];
for ii = 1:numel(fn);
    mm = FastRead(fn{ii});
    mm(isnan(mm))=0;
    
    [a b c] = fileparts(fn{ii});
    if isempty(a); a = pwd; end
    
    ind = 1:size(mm,1);
    ind = setdiff(ind,subset{ii});
    
    II{ii,1} = ind;
    II{ii,2} = 1:numel(ind);
    
    if ii>1
        II{ii,2} = II{ii,2}+II{ii-1,2}(end);
    end
    
    M = [M; zscore(mm(ind,:))];
end
% M(isnan(M))=0;
oM = M;
M = demean(M');

cm = (M'*M)./(size(M,1)-1);
[U,sigma,R] = svd(cm);
    
tmp = abs([max(R); min(R)]);
[trash,i1] = max(tmp); i1(i1==2)=-1;
for kk = 1:size(R,2); R(:,kk) = R(:,kk)*i1(kk); end

RM = R;
E = diag(sigma);
PCs = M*R;

vv = cumsum(E)./sum(E);
vi = find(vv>.9);
pind = 1:vi(1);    

fnpth = [outdir filesep 'TBR_Maps' tag];
if isempty(controls)
    %%
    Tgt = [];
    for ii = 1:numel(targets)
        tmp = resizeVol2(spm_vol(targets{ii}) ,h,3);
        Tgt(:,ii) = tmp(:);
    end
    %Tgt = FastRead(targets)';
    Tgt(isnan(Tgt))=0;
    jx1 = 1:numel(targets);
    jx2 = [];
else
     Tgt = [];
    for ii = 1:numel([targets; controls])
        tmp = resizeVol2(spm_vol(targets{ii}) ,h,3);
        Tgt(:,ii) = tmp(:);
    end    
    %Tgt = [FastRead(targets)' FastRead(controls)'];
    Tgt(isnan(Tgt))=0;
    jx1 = 1:numel(targets);
    jx2 = (numel(targets)+1):size(Tgt,2);
end


b = pinv(PCs(:,pind))*(Tgt);
tcs = RM(:,pind)*b;
F = PCs(:,pind)*b;


try
    rmdir(fnpth, 's')
catch
end
mkdir(fnpth)


save([fnpth filesep 'PCAinfo.mat'], 'tcs','pind','U', 'sigma', 'R', 'E','II','jx1','jx2');

hh = h(1);
hh.dt = [16 0];
V = zeros(hh.dim); 

M = zscore(oM); % Comment this line if you want the B_maps to be rescaled versions of the F_maps.
for gg = jx1
    B = pinv(zscore([tcs(:,gg) tcs(:,jx2)]))*(M);
    [nm1 nm2 nm3] = fileparts(targets{gg});
    hh.fname = [fnpth filesep nm2  '.nii'];
    V(:) = B(1,:);
    spm_write_vol(hh,V);
end


pred = demean(tcs)*(pinv(demean(tcs))*(M));
R2 = SumOfSquares(pred)./SumOfSquares(M);
V(:) = R2;
hh.fname = [fnpth filesep 'R2mapAll.nii'];
spm_write_vol(hh,V);

if ~isempty(controls)
    pred = demean(tcs(:,jx1))*(pinv(demean(tcs(:,jx1)))*(M));
    R2 = SumOfSquares(pred)./SumOfSquares(M);
    V(:) = R2;
    hh.fname = [fnpth filesep 'R2mapTarget.nii'];
    spm_write_vol(hh,V);

    pred = demean(tcs(:,jx2))*(pinv(demean(tcs(:,jx2)))*(M));
    R2 = SumOfSquares(pred)./SumOfSquares(M);
    V(:) = R2;
    hh.fname = [fnpth filesep 'R2mapControl.nii'];
    spm_write_vol(hh,V);
end


for qq = 1:size(II,1)
    if size(II,1)==1
        continue;
    end
    mkdir([fnpth filesep 'Run' sprintf('%0.2d',qq)]);
    
    hh = h(1); hh.dt = [16 0];
    V = zeros(hh.dim);
    for gg = jx1;
        [nm1 nm2 nm3] = fileparts(targets{gg});
        hh.fname = [fnpth filesep 'Run' sprintf('%0.2d',qq) filesep nm2 '.nii'];
        B = pinv(zscore(tcs(II{qq,2},[gg jx2])))*(oM(II{qq,2},:));
        V(:) = B(1,:);
        spm_write_vol(hh,V);
    end
end

cd(hm);
