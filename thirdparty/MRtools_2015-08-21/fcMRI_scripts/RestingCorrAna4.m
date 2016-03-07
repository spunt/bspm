function RestingCorrAna4(fn,seed_spec,name)
%%% Multiple simultaneous Seeds
%%%         {'PCC'  [  0 -53  26] [10]} ...
%%%         {'mPFC' [  0  52 -06] [10]} ...
%%%         {'lLPC' [-48 -62  36] [8]} ...
%%%         {'LRPC' [ 46 -62  32] [8]} ...
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
if ~iscell(fn)
    fn = {fn};
end
% keyboard;
h = spm_vol([fn{1} ',1']);
ind = 1:prod(h.dim);
for ii = 1:length(fn)
    [M h] = openIMG(fn{ii});
    ss = size(M);
    MM{ii} = reshape(M, prod(ss(1:3)),ss(4))';
    ind = intersect(ind,find(~isnan(sum(MM{ii}))));
end

x = [];
for zz = 1:numel(seed_spec)
    for ii = 1:numel(fn)
        if numel(seed_spec{zz})==2 && ischar(seed_spec{zz}{2})
            mask = openIMG(seed_spec{zz}{2});
            i1 = find(mask>0);
            i1 = intersect(i1,ind);
            x{ii}(:,zz) = mean(MM{ii}(:,i1),2);
            lab{zz} = 'MaskImg.nii';
        else
            if numel(seed_spec{zz}) == 3
                D = extractTS(fn{ii}, {seed_spec{zz}});
                i1 = intersect(ind,D.(seed_spec{zz}{1}).vec_loc);
                x{ii}(:,zz) = mean(MM{ii}(:,i1),2);
                lab{zz} = D.(seed_spec{zz}{1}).fn(6:end);
            else
                x{ii}(:,zz) = seed_spec{zz}{2};
                lab{zz} = 'PreSpecTS.nii';
            end
        end
    end
end

% %%% Drop Bad Volumes
% for ii = 1:length(fn)
%     [a b c] = fileparts(fn{ii});
%     if isempty(a)
%         bv = 'ExtraRegressors.mat';
%     else
%         bv = [a filesep 'ExtraRegressors.mat'];
%     end
%     
%     if exist(bv)>0
%         load(bv);
%         vind{ii} = find(sum(R,2)==0)';
%         disp(['Dropping ' num2str(numel(find(sum(R,2)==1))) ' Volumes from the analysis']);
%     else
%         vind{ii} = 1:size(MM{ii},1);
%     end
% end
vind{ii} = 1:size(MM{ii},1);
%%%
V = h(1);
vol = zeros(V.dim)*NaN;

Y = nan(numel([vind{:}]),numel(ind));
last = 0;
mod = []; 
ns = numel(seed_spec)+1;
for ii = 1:length(fn)
    rs = (1:numel(vind{ii}))+last;
    Y(rs,:) = MM{ii}(vind{ii},ind);
    last = rs(end);
    
        mod(end+1:end+numel(vind{ii}),(((ii-1)*ns)+1):ii*ns) = [ones(numel(vind{ii}),1) x{ii}(vind{ii},:)];

end

bb = pinv(mod)*Y;


SSt = Y-repmat(mean(Y),size(Y,1),1);
SSt = sum(SSt.^2);

for ii = 2:ns
    vol(:)=NaN;
    BSS = zeros(1,numel(ind));
    for jj = ii:ns:size(bb,1)
        pp = mod(:,jj)*bb(jj,:);
        BSS = BSS+(sum((pp-repmat(mean(pp),size(pp,1),1)).^2).*sign(bb(jj,:)));
    end
    BSS = (BSS./numel(ii:ns:size(bb,1)))./SSt;
    sgn = sign(BSS);
    vol(ind) = atanh(sqrt(abs(BSS)).*sgn);
    
    [a b c] = fileparts(fn{1});
    if isempty(a)
        V.fname = ['zr2map_' name '_beta_' num2str(ii-1) '_' b lab{ii-1}];
    else
        V.fname = [a filesep 'zr2map_' name '_beta_' num2str(ii-1) '_' b lab{ii-1}];
    end
    V.dt = [16 0];
    spm_write_vol(V,vol);
end


pred = mod*bb;
SSp = pred-repmat(mean(pred),size(pred,1),1);
SSp = sum(SSp.^2);


R = atanh(sqrt(SSp./SSt));

%%%
V = h(1);
vol = zeros(V.dim)*NaN;
vol(ind) = R;
[a b c] = fileparts(fn{1});

if isempty(a)
    V.fname = ['zr2map_' name '_' b '.nii'];
else
    V.fname = [a filesep 'zr2map_' name '_' b '.nii'];
end
V.dt = [16 0];
spm_write_vol(V,vol);
