function [F RES] = estimateError2(y,F)
%%% This script will compute the error terms specified in the Design
%%% Structure F.  
%%%
%%% Inputs:
%%% y = a column vector or matrix of data mated to the design specified in F
%%%
%%% F = Design structure from CreateDesign.m
%%%
%%%
%%% Outputs:
%%% The computed statistics are put in the data structure F and returned.
%%%
%%%
%%% Written by Aaron Schultz (aschultz@martinos.org) 
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

if ~isfield(F,'writeRes');
    F.writeRes=0;
end

if isfield(F,'ind');
    Xind = F.ind;
else
    Xind = 1:size(F.XX,1);
end

if size(F.XX,1) == numel(Xind)
	xx = F.XX;
else
	xx = F.XX(Xind,:);
end

if isfield(F,'WX');
	if size(F.WX,1) == numel(Xind)
    	x = F.WX;
    else
    	x = F.WX(Xind,:);
	end
else
    x = xx;
    x(isnan(xx))=0;
end

wh = find(sum(x)~=0);
x = x(:,wh);
xx = xx(:,wh);

eb = pinv(x)*y;
er = eye(size(x,1))-(x*pinv(x));


ind = contains('^S',F.name);
if isempty(ind);
    FM = x;
else
    FM = x(:,1:ind(1)-1);
end

df = [];
SS = {};
RES = {};
first = [];
for ii = 1:length(F.err)
        if F.err(ii)==-1
            a = eye(size(x,1));
        else
            a = F.FF{F.err(ii)}(Xind,:);
            a = a(:,find(nansum(a)>0));
            a(isnan(a))=0;
        end
        
        %%%% This only works since the error terms are computed in order.
        FM = [FM first];
        a = [first a];
        
        tx1 =  FM;    er1 = eye(size(tx1,1))-(tx1*pinv(tx1));
        tx2 = [FM a]; er2 = eye(size(tx2,1))-(tx2*pinv(tx2));
        
        mxm = er1-er2;
        %SS{ii,1} = LoopEstimate(y,1,mxm);
              
        [z1 z2 z3] = svd(tx1);
        tol = max(size(tx1))*max(abs(diag(z2)))*eps;
        if size(z2,2)>1
            df1 = sum(diag(z2)>tol);
        else
            df1 = 1;
        end
       

        [z1 z2 z3] = svd(tx2);
        tol = max(size(tx2))*max(abs(diag(z2)))*eps;
        if size(z2,2)>1
            df2 = sum(diag(z2)>tol);
        else
            df2 = 1;
        end
        
        df(ii) = df2-df1;
        
        
%%% Will need to double check the below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        if nargout == 2 || F.writeRes == 1;
            [SS{ii,1} RES{ii}] = LoopEstimate(y,1,mxm);
            RES{ii} = RES{ii}./repmat(sqrt(SS{ii}/df(ii))',1,size(RES{ii},2));            
            if F.writeRes==1
                fprintf('\n');
                
                tm = [pwd '/ResAll_00'];
                n = num2str(ii);
                tm((end-length(n))+1:end) = n;
                
                warning off;
                delete([tm '.nii']); delete([tm '.mat']);
                warning on;
                
                v = F.v;
                v.fname = [tm '.nii'];
                
                vol = nan(v.dim);
                for jj = 1:size(RES{ii},2);
                    vol(F.vec) = RES{ii}(:,jj);
                    v.n = [jj 1];
                    spm_write_vol(v,vol);
                end
                clear RES;
                h = spm_vol(v.fname);
                %disp([numel(h) df(ii)]);
                [FWHM,VRpv,R] = spm_est_smoothness(h,spm_vol('AllMask.nii'),[numel(h) df(ii)]);
                F.FWHM{ii} = FWHM;
                %F.ReslInfo{ii} = R;
                [tmpM tmpH] = openIMG('RPV.img');
                tmpH.fname = ['RPV' tm(end-2:end) '.nii'];
                spm_write_vol(tmpH,tmpM);
                delete('RPV.*')
            end
        else
            SS{ii,1} = LoopEstimate(y,1,mxm);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    first = a;
end

F.ErrorSS = SS;
F.ErrorDF = df;

warning off
delete('AllMask.nii')
warning on