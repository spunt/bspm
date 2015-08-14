function [r c cc rr] = MakeContrastMatrix2(F,cols,levs,Xind)
%%% Generate contrast matricies and mixing matrices
%%%
%%% Inputs:
%%%
%%% type = the method of generating the pre contrasts. In general it is
%%% advisable to always use 'a', this seems to be the most general and
%%% stable method.  The other methods will likely be removed in future
%%% versions
%%%
%%% sDM = Sub Design Matrix. here you pass in the columns of the design
%%% matrix corresponding to the effect that you are interested in (each
%%% condition should be uniquely specified by a column).
%%%
%%% levs = the number and distribution of levels.  A four level on way
%%% anova would be specified as [4] wheras a 2x2 interation would be
%%% specified as [2 2].  The order should be the same as in the model. That
%%% is if you put in the 3 level factor before a 2 level factor then the
%%% levs order would be [3 2] and not [2 3].  Order here implies the order
%%% of looping to create the conditions "for ii = 1:3; for jj = 1:2; ..."
%%% will produce a different design matrix patternd than "for ii = 1:2; for
%%% jj = 1:3; ..." so be sure that the order is correct or the contrast
%%% will not be correct.
%%%
%%% xx = the full design matrix with NaN's in place of zeros.
%%%
%%% x = the whitened designed matrix or if no variance correction, it will
%%% be the same as xx only with zeros instead of NaN's
%%%
%%% CovarCols = a vector of ones and zeros coding whether or not each
%%% column of the design matrix is a covariate.  If not specified it is
%%% assumed that there are no covaraites or that covariates should be
%%% included in the group/condition definitions (e.g. don't remove the
%%% effect of covariate from group).
%%%
%%% Outputs:
%%%
%%% r = the mixing matrix to compute the SS for the contrast. This is used
%%% in SS = b'*x'*m*x*b.  m is the mixing matrix which is r-rr.
%%%
%%% rr =  the residual forming matrix for the full model
%%%
%%% c = the contrasts vector/matrix
%%%
%%% cc = the pre-contrast vector/matrix.  
%%%
%%% Note: c is comuted from cc.  Based on the number of levels and their
%%% distribution across factors, the contrast c is created by folding cc
%%% and the differentiating across the proper dimensinos
%%% 
%%%
%%% Written by Aaron Schultz, May 5th, 2010 (aschultz@martinos.org)
%%% Copyright (C) 2011,  Aaron P. Schultz
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
if nargin<4 || isempty(Xind)
    Xind = 1:size(F.XX,1);
end

xx = F.XX(Xind,:);
% dropcols = find(nansum(~isnan(xx))>0);
x = xx; x(isnan(x))=0;

CovarCols = F.CovarCols;


ndx = sub2ind_aps(size(F.FF),find(F.isBet + F.isWith == 1));
partition = F.FF{ndx};
partition = partition(Xind,:);
% if ~isfield(F,'parts') || isempty(F.parts);
    parts = [];
    for ii = 1:size(partition,2);
        i1 = find(partition(:,ii)==1);
        parts(ii,:) = nansum(xx(i1,:))./numel(i1); %/size(tmp,2)
    end
    parts = parts';
% else
%     parts = F.parts;
% end

parts(find(CovarCols==1),:) = 1;


if isvector(cols) && size(cols,1)==1 && iscell(cols)
    sDM = [];
    for ii = 1:numel(cols)
        ti = find(mean(isnan(xx(:,cols{ii})),2)==1);
        sDM(:,ii) = nansum(xx(:,cols{ii}),2);
        sDM(ti,ii) = NaN;
    end
elseif isvector(cols) && size(cols,1)==1 && ~iscell(cols)
    sDM = xx(:,cols);    
else
    sDM = cols(Xind,:);
end

cc = [];
covCon = 0;  match = [];
for ii = 1:size(sDM,2)
    a1 = sDM(:,ii);
    match = [];
    for jj = 1:size(x,2)
        ti = find(~isnan(xx(:,jj)));
        if all(a1(ti)==xx(ti,jj));
            match(end+1) = jj;
        end
    end
    
    if any(CovarCols(match))
        covCon = 1;
        cc(:,ii) = zeros(size(x,2),1);
        cc(match,ii) = 1;
    end
end

if covCon==1
    for ii = 1:size(cc,2)
        cc(:,ii)= cc(:,ii).*F.weights';
    end
else
    cc = [];
    count = 0;
    cc = zeros(size(x,2),size(cols,2));
    for ii = 1:size(sDM,2)
        %i1 = sDM(:,ii)==1;
        i1 = ~isnan(sDM(:,ii));
        
        i2 = nanmean(partition(i1,:))>0;
        cc(:,ii) = mean(parts(:,i2),2);
    end
    
    % added May 7th 2014.  Need to remove weights for covariates
    cc(CovarCols==1,:)=0;
end

if numel(levs)>1 && levs(1)~=0
    tmp = reshape(cc,[size(cc,1) levs(end:-1:1)]);
else
    tmp = cc;
end
tmp = squeeze(tmp);


if size(tmp,2)==1 || levs(1) == 0
    c = tmp;
else
    c = differencer(tmp)*-1;
end
   
c0 = eye(size(x,2))-(c*pinv(c));
x0 = x*c0;
r = eye(size(x0,1))-(x0*pinv(x0));

if nargout == 4
    c0 = eye(size(x,2))-(cc*pinv(cc));
    x0 = x*c0;
    rr = eye(size(x0,1))-(x0*pinv(x0));
end
end

function out = differencer(tmp,count)

if nargin == 1;
    count = numel(size(tmp));
end

out = diff(tmp,1,count);

if count ~= 2;
    out = differencer(out,count-1);
else
    if numel(size(out))>2   
        ss = size(out);
        out = reshape(out,size(out,1),prod(ss(2:end)));
    end
    return
end
end
