function [vec L R] = CreateConVec(left,right,SPM,WeightWithinConds,BlockThresh)
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

L.durs = [];
L.counts = [];
R.durs = [];
R.counts = [];

conds = [];
for ii = 1:length(SPM.Sess);
    for jj = 1:length(SPM.Sess(ii).U);
        conds{end+1} = SPM.Sess(ii).U(jj).name{1};
    end
end
conds = unique(conds);

if isstruct(left);
    Lvec = zeros(size(SPM.xX.name));
    Ldur = zeros(size(SPM.xX.name));
    Linds = [];
    Lnames = cell(size(SPM.xX.name));
    c = 0;
    for ii = 1:length(SPM.Sess);
        for jj = 1:length(SPM.Sess(ii).U);
            c = c+1;
            if contains(left.mid,SPM.Sess(ii).U(jj).name);
                tr(ii,jj) = numel(SPM.Sess(ii).U(jj).ons);
                tind = contains(['Sn\(' num2str(ii) '\).*' SPM.Sess(ii).U(jj).name{1}], SPM.xX.name);

                Lvec(tind) = numel(SPM.Sess(ii).U(jj).ons);
                Ldur(tind) = mean(SPM.Sess(ii).U(jj).dur);
                Linds(tind) = tind;
                Lnames(tind) = {SPM.Sess(ii).U(jj).name{1}};
            end
        end
    end
    
    ind = contains([left.pre left.mid left.post],SPM.xX.name);
    Lvec(setdiff(1:numel(SPM.xX.name),ind))=0;
    Lcounts = Lvec(ind);
    
    i1 = ind(find(Ldur(ind)>=0 & Ldur(ind)<BlockThresh));
    i2 = ind(find(Ldur(ind)>=0 & Ldur(ind)>=BlockThresh));

    LvecE = zeros(size(Lvec));
    if WeightWithinConds
        [a,b,c] = unique(Lnames(i1));
        for ii = 1:length(a)
            ti = find(c==ii);
            LvecE(i1(ti)) = (Lvec(i1(ti))./sum(Lvec(i1(ti))))*(1/numel(a));
        end
    else
        LvecE(i1) = Lvec(i1)./sum(Lvec(i1));
        LvecE(isnan(LvecE))=0;
    end
    
    
    lvecT = double(Lvec(i2)>0);
    LvecB = zeros(size(Lvec));
    if WeightWithinConds
        [a,b,c] = unique(Lnames(i2));
        for ii = 1:length(a)
            ti = find(c==ii);
            LvecB(i2(ti)) = (lvecT(ti)./numel(ti)).*(1/numel(a));
        end
    else
        LvecB(i2) = lvecT./numel(lvecT);
    end
    
    if isempty(LvecE) || sum(LvecE)==0;
        Lvec = LvecB;
    elseif isempty(LvecB) || sum(LvecB)==0;
        Lvec = LvecE;
    else
        Lvec = (LvecE+LvecB)/2;
    end
        
    Lnames = Lnames(ind);
    Ldur = Ldur(ind);
    ind = find(Ldur>BlockThresh);
    Lcounts(ind) = inf;

    L.durs = Ldur;
    L.counts = Lcounts;
else
    Lvec = zeros(size(SPM.xX.name));
end

if isstruct(right)
    Rvec = zeros(size(SPM.xX.name));
    Rdur = zeros(size(SPM.xX.name));
    Rinds = [];
    Rnames = cell(size(SPM.xX.name));
    c = 0;
    for ii = 1:length(SPM.Sess);
        for jj = 1:length(SPM.Sess(ii).U);
            c = c+1;
            if contains(right.mid,SPM.Sess(ii).U(jj).name);
                tind = contains(['Sn\(' num2str(ii) '\).*' SPM.Sess(ii).U(jj).name{1}], SPM.xX.name);
                Rvec(tind) = numel(SPM.Sess(ii).U(jj).ons);
                Rdur(tind) = mean(SPM.Sess(ii).U(jj).dur);
                Rinds(tind) = tind;
                Rnames(tind) = {SPM.Sess(ii).U(jj).name{1}};
            end
        end
    end
    
    ind = contains([right.pre right.mid right.post],SPM.xX.name);
    Rvec(setdiff(1:numel(SPM.xX.name),ind))=0;
    Rcounts = Rvec(ind);
    
    i1 = ind(find(Rdur(ind)>=0 & Rdur(ind)<BlockThresh));
    i2 = ind(find(Rdur(ind)>=0 & Rdur(ind)>=BlockThresh));

    RvecE = zeros(size(Rvec));
    if WeightWithinConds
        [a,b,c] = unique(Rnames(i1));
        for ii = 1:length(a)
            ti = find(c==ii);
            RvecE(i1(ti)) = (Rvec(i1(ti))./sum(Rvec(i1(ti))))*(1/numel(a));
        end
    else
        RvecE(i1) = Rvec(i1)./sum(Rvec(i1));
        RvecE(isnan(RvecE))=0;
    end
    
    
    rvecT = double(Rvec(i2)>0);
    RvecB = zeros(size(Rvec));
    if WeightWithinConds
        [a,b,c] = unique(Rnames(i2));
        for ii = 1:length(a)
            ti = find(c==ii);
            RvecB(i2(ti)) = (rvecT(ti)./numel(ti)).*(1/numel(a));
        end
    else
        RvecB(i2) = rvecT./numel(rvecT);
    end
    
    if isempty(RvecE) || sum(RvecE)==0;
        Rvec = RvecB;
    elseif isempty(RvecB) || sum(RvecB)==0;
        Rvec = RvecE;
    else
        Rvec = (RvecE+RvecB)/2;
    end
    
    Rnames = Rnames(ind);
    Rdur = Rdur(ind);
    ind = find(Rdur>BlockThresh);
    Rcounts(ind) = inf;

    R.durs = Rdur;
    R.counts = Rcounts;
else
    Rvec = zeros(size(SPM.xX.name));
end

vec = Lvec-Rvec;


if (sum(abs(vec))-1)<1e-10 || (sum(abs(vec))-2)<1e-10
    
else
    error('This Contrast was not generated correctly');
end
