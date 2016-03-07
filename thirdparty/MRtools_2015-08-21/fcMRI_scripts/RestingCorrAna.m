function RestingCorrAna(fn,seed_spec,dv,params,writeLocal)
%%% This Script will perform a seed based resting connectivity analysis.
%%%
%%% fn:             The name of the nifti file to be analyzed. (Can be more
%%%                 than one file if there were multiple acquisitions for
%%%                 the session.  If multiple are specified you will get a
%%%                 seed map for each run and one map across all runs.
%%%
%%% seed_spec:      A cell aray of cell arrays.  Each inner cell array
%%%                 consists of {'Name', [mni locations], [seed diameter]}
%%%                 e.g. 
%%%                 {{'PCC'  [  0 -53  26] [10]} ...
%%%                 {'mPFC' [  0  52 -06] [10]} ...
%%%                 {'lLPC' [-48 -62  36] [8]} ...
%%%                 {'LRPC' [ 46 -62  32] [8]}}
%%%                 The above will create seed maps for each of the four
%%%                 specified seeds.
%%%                 Seeds can also be specified from masks:
%%%                     {'ThisMask' 'path/to/special/mask.nii'}
%%%                 Or times series can be passed in directly:
%%%                     {'TimeCourse' X}
%%%
%%% params:         params (optional) is a string that specifies a param file for use
%%%                 with RestPostProc.m to perform custom preprocessing on
%%%                 the fly.
%%%
%%% Example:
%%%     seed_spec = {{'PCC' [  0 -53  26] [10]} ...
%%%                 {'mPFC' [  0  52 -06] [10]} ...
%%%                 {'lLPC' [-48 -62  36] [8]}  ...
%%%                 {'LRPC' [ 46 -62  32] [8]}};
%%%     RestingCorrAna('res_ff_ss_nn_rr_st_dv_restingState_1.nii',seed_spec);
%%%     or
%%%     RestingCorrAna('ss_nn_rr_st_dv_restingState_1.nii',seed_spec, 'PostParams');
%%%
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%%
%%% Last Updated Dec. 11 2012;
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

if ~iscell(fn)
    fn = {fn};
end

if nargin>2 && ~isempty(dv)
    if ~iscell(dv)
        dv = {dv};
    end
end


h = spm_vol([fn{1} ',1']);
ind = 1:prod(h.dim);
MM = [];
if nargin >= 2 || isempty(params)
    for ii = 1:length(fn)
        [M h] = openIMG(fn{ii});
        ss = size(M);
        MM{ii} = reshape(M, prod(ss(1:3)),ss(4))';
        ind = intersect(ind,find(~isnan(sum(MM{ii}))));
    end
elseif nargin == 4
    hm = pwd;
    for ii = 1:length(fn)
        [a1 a2 a3] = fileparts(fn{ii});
        if ~isempty(a1);
            cd(a1);
            [M tmp] = RestPostProc([a2 a3],params);
            fn{ii} = [a1 filesep tmp];
            cd(hm)
        else
        [M fn{ii}] = RestPostProc(fn{ii},params);
        end
        
        
        ss = size(M);
        MM{ii} = M;
        ind = intersect(ind,find(~isnan(sum(MM{ii}))));
    end
end
    
x = [];
for zz = 1:numel(seed_spec)
    for ii = 1:numel(fn)
        if numel(seed_spec{zz})==2 && ischar(seed_spec{zz}{2})
            mask = openIMG(seed_spec{zz}{2});
            i1 = find(mask>0);
            i1 = intersect(i1,ind);
            x{zz}{ii} = mean(MM{ii}(:,i1),2);
            lab{zz} = 'MaskImg.nii';
        elseif numel(seed_spec{zz})==2 && isnumeric(seed_spec{zz}{2})
            x{zz}{ii} = seed_spec{zz}{2};
            lab{zz} = 'PreSpecTS.nii';
        elseif numel(seed_spec{zz}) == 3
                [ml vi] = getMatCoord(h(1),seed_spec{zz}{2},seed_spec{zz}{3});

                i1 = intersect(ind,vi);
                x{zz}{ii} = mean(MM{ii}(:,i1),2);
                
                t1 = num2str(seed_spec{zz}{2}); t1 = regexprep([t1], '  ', ' '); t1 = regexprep([t1], '  ', ' '); t1 = regexprep([t1], '  ', ' ');
                lab{zz} = [num2str(seed_spec{zz}{3}) 'mm_' regexprep([num2str(t1)], ' ', '_') '.nii'];
        end
    end
end
% keyboard;
%%% Drop Bad Volumes
for ii = 1:length(fn)
    [a b c] = fileparts(fn{ii});
    
    if nargin>2 && ~isempty(dv);
        if isempty(a)
            bv = dv{ii};
        else
            bv = [a filesep dv{ii}];
        end
    end
    
    if exist('bv','var')>0
        load(bv);
        vind{ii} = find(sum(R,2)==0)';
        disp(['Dropping ' num2str(numel(find(sum(R,2)==1))) ' Volumes from the analysis']);
    else
        vind{ii} = 1:size(MM{ii},1);
    end
end

%%%
V = h(1);
vol = nan(V.dim);

Y = nan(numel([vind{:}]),numel(ind));
last = 0;
for ii = 1:length(fn)
    rs = (1:numel(vind{ii}))+last;
    Y(rs,:) = MM{ii}(vind{ii},ind);
    last = rs(end);
end

for zz = 1:numel(seed_spec)
    mod = [];
    for ii = 1:length(fn)
        
        mod(end+1:end+numel(vind{ii}),(((ii-1)*2)+1):ii*2) = [ones(numel(vind{ii}),1) x{zz}{ii}(vind{ii})];
        
        if numel(fn)>1
            bb = pinv(zscore(x{zz}{ii}(vind{ii})))*zscore(MM{ii}(vind{ii},ind));
            vol(:)=NaN;
            vol(ind) = atanh(bb);
            [a b c] = fileparts(fn{ii});
            if isempty(a)
                V.fname = ['zmap_' seed_spec{zz}{1} '_' b '_' lab{zz}];
            else
                V.fname = [a filesep 'zmap' num2str(ii) '_' seed_spec{zz}{1} '_' b '_' lab{zz}];
            end
            V.dt = [16 0];
            spm_write_vol(V,vol);
        end
    end
    
    %%%
    if numel(fn)>1
        bb = pinv(mod)*Y;
        
        SSc = zeros(1,size(bb,2));
        c = 0;
        
        for ii = 2:2:numel(fn)*2
            c = c+1;
            tmp = mod(:,ii)*bb(ii,:);
            tmp = tmp-repmat(mean(tmp),size(tmp,1),1);
            
            SSc = SSc + (sum(tmp.^2).*sign(bb(ii,:)));
        end
        
        
        SSt = Y-repmat(mean(Y),size(Y,1),1);
        SSt = sum(SSt.^2);
        
        sgn = sign(SSc);
        R2 = abs(SSc)./SSt;
        R = atanh(sgn.*sqrt(R2));
        
    else
        R = atanh(pinv(zscore(mod(:,2)))*zscore(Y));
    end
    
    %%%
    V = h(1);
    vol = nan(V.dim);
    vol(ind) = R;
    
    
    [a b c] = fileparts(fn{1});
    
    if isempty(a) || nargin == 5
        V.fname = ['zmap_' seed_spec{zz}{1} '_' b(1:end-1) '_' lab{zz}];
    else
        V.fname = [a filesep 'zmap_' seed_spec{zz}{1} '_' b(1:end-1) '_' lab{zz}];
    end
    V.dt = [16 0];
%     keyboard;
    spm_write_vol(V,vol);
end
