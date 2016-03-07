function  DualReg(fn,templates,subset,tag,outdir)
% This function implements the TBR functional mapping procedure.  An SPM 
% installation is required.
%
% FN - A cell array of input files (expects 4D-Nifti files)
%
% TEMPLATES - a cell array of file paths to template maps.
%
% SUBSET - an OPTIONAL index of time points to include.
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

if nargin<4
    tag = '';
end

if nargin<5 || isempty(outdir);
    outdir = hm;
end

if ischar(fn);
    fn = cellstr(fn);
end

if nargin<3 || isempty(subset); 
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
M(isnan(M))=0;

% keyboard;
fnpth = [outdir filesep 'OrthoNormOn20_675_DR' tag];
Tgt = FastRead(templates)';

%%% First Regression
tcs = (pinv(zscore(Tgt))*(M'))';
%%% Second regression
b = pinv(zscore(tcs))*M;


if ~exist(fnpth)
    mkdir(fnpth)
end

%%% Save the estimate time-courses
dr_tcs = tcs;
save([fnpth filesep 'TimeCourses.mat'], 'dr_tcs');

hh = h(1);
hh.dt = [16 0];
V = zeros(hh.dim);

%%% Create the beta maps and write them to file.
for gg = 1:size(b,1)
    [nm1 nm2 nm3] = fileparts(templates{gg});
    hh.fname = [fnpth filesep nm2 '.nii'];
    V(:) = b(gg,:);
    spm_write_vol(hh,V);
end
%%% Create an R^2 map
pred = zscore(tcs)*b;
R2 = SumOfSquares(pred)./SumOfSquares(M);
V(:) = R2;
hh.fname = [fnpth filesep 'R2map.nii'];
spm_write_vol(hh,V);

%%%%%%%%%%%%%
for qq = 1:size(II,1)
    if size(II,1)==1
        continue;
    end
    mkdir([fnpth filesep 'Run' sprintf('%0.2d',qq)]);
    
    b = pinv(zscore(tcs(II{qq,2},:)))*(M(II{qq,2},:));
    
    hh = h(1); hh.dt = [16 0];
    V = zeros(hh.dim);
    for gg = 1:size(b,1)
        [nm1 nm2 nm3] = fileparts(templates{gg});
        hh.fname = [fnpth filesep 'Run' sprintf('%0.2d',qq) filesep nm2 '.nii'];
        V(:) = b(gg,:);
        spm_write_vol(hh,V);
    end
end

cd(hm);
