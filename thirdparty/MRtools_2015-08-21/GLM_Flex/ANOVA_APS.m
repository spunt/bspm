function OUTPUT = ANOVA_APS(dat,mod,posthocs,runModel,dm,asEffects)
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

%%% SEB stuff.
% % X = mod.X;
% % rr = eye(size(X,1))-(X*pinv(X));
% % MSE = LoopEstimate(y,1,rr)./mod.RFMs(1).EDF;
% % b = pinv(X)*y;
% % cv = MSE*inv(X'*X);
% % seb = sqrt(diag(cv));
% % t = b./seb;

OUTPUT.data = dat;
iscontinuous = [];

isconstant = 0;
supress = 0;

if nargin<3
    posthocs = [];
else
    if any(size(posthocs)==1) && any(size(posthocs)>2)
       posthocs = posthocs(:); 
    end
end

if nargin<4 || isempty(runModel)
    runModel = 1;
end

if runModel==3
    runModel=1;
    supress = 1;
end

if nargin<5 || isempty(dm)
    dm = 1;
    % disp('All covariates will be demeaned');
end

if nargin<6 || isempty(asEffects) || runModel==1
    asEffects = 1;
end

%%% Split the model into sets of additive terms
ind = find(mod=='~');
if numel(ind)>1
    error('Only 1 ~ should be present in model specification');
end
if ~isempty(ind)
    varname = strtrim(mod(1:ind(1)-1));
    Nobs = size(dat.(varname),1);
else
    runModel = 0;
    ind = 0;
end

mod2 = mod(ind(1)+1:end);

mod3 = regexp(mod2,'+','split');
for ii = 1:numel(mod3);
    mod3{ii} = strtrim(mod3{ii});
end

%%% Get a list of the unique terms in the model : Not inclusing error specification
terms = [];
for ii = 1:numel(mod3)
    if ~isempty(strfind(mod3{ii},'random'))
        continue
    end
    
    if ~isempty(regexp(mod3{ii},'*')) && ~isempty(regexp(mod3{ii},':'))
        error('You cannot specify both * and : within the same term');
    end
    
    tmp = regexp(mod3{ii},'[*:]','split');
    for jj = 1:numel(tmp);
        terms{end+1} = strtrim(tmp{jj});
    end
end
terms= unique(terms,'stable');
if exist('varname','var')==0
    Nobs = size(dat.(terms{1}),1);
end

%%% now create the dummy coded main effect terms
Effects = [];
for ii = 1:numel(terms);
    if isnumeric(dat.(terms{ii}))
        iscontinuous(ii) = 1;
        if dm == 1
            Nobs = numel(dat.(terms{ii}));
            Effects{ii} = demean(dat.(terms{ii})); % + ((rand(Nobs,1)-.5)*1e-16); % Real values covariates should never have a value of 0 will mess up contrast later on.
        else
            Effects{ii} = dat.(terms{ii});
        end
    else
        iscontinuous(ii) = 0;
        Effects{ii} = makedummy(dat.(terms{ii}),asEffects);
        if size(Effects{ii},2)==1 && all(Effects{ii}==1)
            isconstant = 1;
        end
    end
end

%%% Now create the interaction terms
for ii = 1:numel(mod3);
    if ~isempty(strfind(lower(mod3{ii}),'random'));
        continue
    end
    tt = regexp(mod3{ii},'[*:]' ,'split');
    tt = strtrim(tt);
    
    for jj = 2:numel(tt); 
        sets = combnk(1:numel(tt),jj);
        for kk = 1:size(sets,1)
            which = sets(kk,:);
            
            if any(iscontinuous(which))
                iscontinuous(end+1) = 1;
            else
                iscontinuous(end+1) = 0;
            end
            
            [i1 i2] = matchAll(tt(which),terms);
            Effects{end+1} = FactorCross(Effects(i2));

            name = [];
            for ll = 1:numel(which);
                name = [name tt{which(ll)} ':'];
            end
            terms{end+1} = name(1:end-1);
        end
    end
end


%%% Now it's time to setup the error terms.
ErrorTerms = []; eTerms = [];
if isempty(contains('random',{mod}))
    eTerms{1} = 'Observations';
    ErrorTerms{1} = makedummy(1:Nobs,asEffects);
    random = 'Observations';
else
    
    tmp = regexp(mod2,'random','split');
    if numel(tmp)>2
        error('Only 1 random effect may be specified.  If you need more look to an LME model');
    end
    err = strtrim(tmp{2});
    err = strtrim(regexprep(err,'[()]',''));
    
    tmp = regexp(err,'\|','split');
    random = strtrim(tmp{1});
    nesting = strtrim(tmp{2});
    np = regexp(nesting,'*','split');
    
    
    eTerms{end+1} = tmp{1};
    ErrorTerms{end+1} = makedummy(dat.(tmp{1}),asEffects);
    
    for ii = 1:numel(np)
        sets = combnk(1:numel(np),ii);
        for jj = 1:size(sets,1)
            [a b] = matchAll(np(sets(jj,:)), terms);
            ErrorTerms{end+1} = FactorCross([ErrorTerms{1} Effects(b)]);
            
            %%% Random effects with : specification?
            name = [random ':'];
            for kk = 1:size(sets(jj,:),2)
                name = [name np{sets(jj,kk)} ':'];
            end
            
            eTerms{end+1} = name(1:end-1);
        end
    end
end

%%% If there is a random effect, make sure that there are no missing
%%% observations.
if numel(eTerms)>1
    RF = dat.(eTerms{1});
    S1 = regexp(eTerms{end},[eTerms{1} ':'],'split');
    S2 = regexp(S1{end},':','split');
    ps = {};
    %keyboard;
    for ii = 1:numel(S2);
        if isnumeric(dat.(S2{ii}))
            continue
        end
        ps{ii} = makedummy(dat.(S2{ii}),0);
    end
    
    if ~isempty(ps);
        
        combo = FactorCross(ps);
        
        missing = {}; %drop = [];
        list = unique(RF);
        for ii = 1:numel(list);
            try
                i2 = find(nominal(RF)==list(ii));
            catch
                i2 = strmatch(list{ii},RF,'exact');
            end
            
            if ~all(sum(combo(i2,:)))
                missing{end+1,1} = list{ii};
            end
        end
        
        if ~isempty(missing)
            disp(missing)
            error('At least one observation is missing for each level of the random factor listed above');
        end
    else
        disp('Assuming that everything is ok');
    end
end

%%% Which terms should go into the model?
ModelTerms = []; whichTerms = [];
for ii = 1:numel(mod3)
    if ~isempty(regexp(mod3{ii},'random'));
        continue
    end
    if isempty(regexp(mod3{ii},'*'))
        ModelTerms{end+1} = mod3{ii};
        whichTerms(end+1) = strmatch(ModelTerms{end},terms,'exact');
    else
        tt = regexp(mod3{ii},'*','split');
        for jj = 1:numel(tt);
            sets = combnk(1:numel(tt),jj);
            for kk = 1:size(sets,1)
                which = sets(kk,:);
                name = [];
                for ll = 1:numel(which);
                    name = [name tt{which(ll)} ':'];
                end
                ModelTerms{end+1} = regexprep(name(1:end-1),' ','');
                tmp = strmatch(ModelTerms{end},terms,'exact');
                whichTerms(end+1) = tmp(1);
            end
        end
    end
end
[ModelTerms idx] = unique(ModelTerms,'stable');
whichTerms = whichTerms(idx);
Effects = Effects(whichTerms);
iscontinuous = iscontinuous(whichTerms);

%%% Now partition the effects and put them with the proper error term
ti = 1:numel(ModelTerms);
DesignParts = [];
for ii = numel(eTerms):-1:1
    
    tmp = eTerms{ii};
    
    tmp = regexp(tmp,[random ':'],'split');
    
    if numel(tmp)==1
        if ~isempty(ti)
            DesignParts{ii,1} = ti;
            continue
        else
            continue
        end
    end
    
    tmp = tmp{2};
        
    i1 = contains(tmp,ModelTerms(ti));
%     %%% New and Experimental
%     for jj = 1:numel(Effects)
%         if iscontinuous(jj);
%             i1 = [i1(:); jj];
%         end
%     end
%     i1 = unique(i1,'stable');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DesignParts{ii,1} = ti(i1);
    
    ti = setdiff(ti,ti(i1));
end

%%%  Now Run/Setup the Model
RFMs = [];
if nargin>2 && ~isempty(posthocs)
    spacer = size(char([ModelTerms 'Error' posthocs(:,1)']),2)+4;
else
    spacer = size(char([ModelTerms 'Error' ]),2)+4;
end
   
if isconstant
    C = [];
else
    C = ones(Nobs,1);
end

effectResults = [];

for ii = 1:size(DesignParts,1);
    if isempty(DesignParts{ii});
        continue
    end
    
    if runModel && ~supress
        fprintf(['\n%-' num2str(spacer) 's%s\t%s\t%s\t%s\t%s\n'],'','df','Sum Sq', 'Mean Sq','F-value', 'P-value');
    end
    %%% Compute the error    
    tx1 = [C Effects{DesignParts{ii}} ErrorTerms{ii}];
    tx2 = [C Effects{DesignParts{ii}}];

    r1 = eye(size(tx1,1))-(tx1*pinv(tx1));
    r2 = eye(size(tx2,1))-(tx2*pinv(tx2));
    
    edf = ResidualDFs(tx1)-ResidualDFs(tx2);
    
    RFMs(ii).EDF = edf;
    RFMs(ii).name = eTerms{ii};
    RFMs(ii).tx1 = tx1;
    RFMs(ii).tx2 = tx2;
    
    if runModel
        RSS = LoopEstimate(dat.(varname),1,r2-r1);
    end

    All = DesignParts{ii};

    for jj = 1:numel(DesignParts{ii})
        in = DesignParts{ii}(jj);
        out = setdiff(DesignParts{ii},in);

        tx1 = [C Effects{All}];
        tx2 = [C Effects{out}];
        dif = [Effects{in}];
        %if ~isempty(contains('aschultz',{UserTime})); keyboard; end

            
        r1 = eye(size(tx1,1))-(tx1*pinv(tx1));
        r2 = eye(size(tx2,1))-(tx2*pinv(tx2));
        if isempty(r2)
            r2 = eye(Nobs);
        end
        
        df = ResidualDFs(tx1)-ResidualDFs(tx2);
        
        RFMs(ii).Effect(jj).df = df;
        RFMs(ii).Effect(jj).name = ModelTerms{DesignParts{ii}(jj)};
        RFMs(ii).Effect(jj).tx1 = tx1;
        RFMs(ii).Effect(jj).tx2 = tx2;
        RFMs(ii).Effect(jj).dif = dif;
            
        if runModel 
            RFMs(ii).MSE =RSS/edf;
             
            SS = LoopEstimate(dat.(varname),1,r2-r1);
            F = (SS/df)/(RSS/edf);
            p = 1-spm_Fcdf(F,df,edf);
            
            effectResults{end+1,1} = {ModelTerms{in},df,SS,SS/df,F,p};
            
            if ~supress
                fprintf(['%-' num2str(spacer) 's%0.0f\t%6.3f\t%6.3f\t%6.3f\t%1.6f\n'],ModelTerms{in},df,SS,SS/df,F,p);
            end
        end
    end
    if runModel && ~supress
        fprintf(['%-' num2str(spacer) 's%0.0f\t%6.3f\t%6.3f\n'],'Error',edf,RSS,RSS/edf);
    end
end


if nargin>2 && ~isempty(posthocs)
    if runModel && ~supress
        fprintf('\nPost Hocs:');
        fprintf(['\n%-' num2str(spacer) 's%s\t%s\t%s\t%s\t%s\n'],'','df','Sum Sq', 'Mean Sq','F-value', 'P-value');
    end
    for ii = 1:size(posthocs,1)
        [ttx1 ttx2 ET df track] = posthoc(posthocs{ii,1},[C Effects{:}],dat,eTerms);
        
        tx1 = RFMs(ET).tx1;
        tx2 = RFMs(ET).tx2;
        
        r1 = eye(size(tx1,1))-(tx1*pinv(tx1));
        r2 = eye(size(tx2,1))-(tx2*pinv(tx2));
        
        edf = ResidualDFs(tx1)-ResidualDFs(tx2);
        
        if runModel
            RSS = LoopEstimate(dat.(varname),1,r2-r1);
        end
        
        
        tx1 = ttx1;
        tx2 = ttx2;
        r1 = eye(size(tx1,1))-(tx1*pinv(tx1));
        r2 = eye(size(tx2,1))-(tx2*pinv(tx2));
        if isempty(r2)
            r2 = eye(Nobs);
        end
        
        df = ResidualDFs(tx1)-ResidualDFs(tx2);
        
        if size(posthocs,2)>1 && ~isempty(posthocs{ii,2});
            name = [posthocs{ii,2} ];
        else
            name = [posthocs{ii,1} ];
        end
        
        RFMs(ET).Effect(end+1).df = df;
        RFMs(ET).Effect(end).name = name;
        RFMs(ET).Effect(end).tx1 = tx1;
        RFMs(ET).Effect(end).tx2 = tx2;
        RFMs(ET).Effect(end).dif = [];
        RFMs(ET).Effect(end).pieces = track;
        if runModel
            SS = LoopEstimate(dat.(varname),1,r2-r1);
            F = (SS/df)/(RSS/edf);
            p = 1-spm_Fcdf(F,df,edf);
            
            effectResults{end+1,1} = {name,df,SS,SS/df,F,p};
            if ~supress
                fprintf(['%-' num2str(spacer) 's%0.0f\t%6.3f\t%6.3f\t%6.3f\t%1.6f\n'],name,df,SS,SS/df,F,p);
            end
        end
        if runModel && ~supress
            fprintf(['%-' num2str(spacer) 's%0.0f\t%6.3f\t%6.3f\n'],['ET #' num2str(ET)],edf,RSS,RSS/edf);
        end
    end
end

OUTPUT.EffTerms = ModelTerms;
OUTPUT.ErrTerms = eTerms;
OUTPUT.EffParts = Effects;
OUTPUT.ErrParts = ErrorTerms;
OUTPUT.X = [C Effects{:} ErrorTerms{1:end-1}];
OUTPUT.EffMod = [C Effects{:}];
OUTPUT.ErrMod = [C ErrorTerms{:}];
OUTPUT.DesignParts = DesignParts;
OUTPUT.RFMs = RFMs;
OUTPUT.const = isconstant;
OUTPUT.Results = effectResults;
% warning off
% OUTPUT.Beta = OUTPUT.X\dat.(varname);
% warning on
try
OUTPUT.Y = dat.(varname);
end

if runModel && ~supress
    TSS =  SumOfSquares(dat.(varname));
    pred = PredictedData(dat.(varname),OUTPUT.X);
    ESS =  SumOfSquares(pred);
    RSS = TSS-ESS;
    %if ~isempty(contains('aschultz',{UserTime})); keyboard; end
    mdf = ResidualDFs(OUTPUT.X)-1;
    rdf = size(pred,1)-mdf-1; 
    
    R2 = ESS/TSS;
    aR2 = R2 - ( (1-R2) * ( (mdf-1)/rdf) );
    F = (ESS/mdf)/(RSS/rdf);
    %F = (R2/mdf)/((1-R2)/rdf);
    
    if isinf(F)
        p=0;
    else
        p = 1-spm_Fcdf(F,mdf,rdf);
    end
    
    OUTPUT.ModelSummary.Y = dat.(varname);
    OUTPUT.ModelSummary.TSS = TSS;
    OUTPUT.ModelSummary.ESS = ESS;
    OUTPUT.ModelSummary.RSS = RSS;
    OUTPUT.ModelSummary.R2 =  R2;
    OUTPUT.ModelSummary.aR2 = aR2;
    OUTPUT.ModelSummary.df =  [mdf edf];
    OUTPUT.ModelSummary.F =  F;
    OUTPUT.ModelSummary.p =  p;
    

    s1 = ['Model R^2 = ' sprintf('%0.3f',R2) '; Adjusted R^2 = ' sprintf('%0.3f',aR2)];
    s2 = ['Model Stats: F(' num2str(mdf) ',' num2str(rdf) ') = ' sprintf('%0.3f',F) '; p = ' sprintf('%0.5d',p)];
    fprintf('\n%s\n%s\n\n',s1,s2)
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tx1 tx2 ET df track] = posthoc(contrast,X,data,ErrTerms)
    cond4 = [];
    extra = [];

    flag = [];
    step1 = regexprep(contrast,' ', '');
    step2 = regexp(step1,'&','split');
    for mm = 1:numel(step2);    
        step3 = regexp(step2{mm},'#','split');

        cond3 = [];
        track1 = {}; track2 = {};
        for jj = 1:numel(step3)
            cond2 = [];

            step5 = regexp(step3{jj},'\|','split');
            cond = [];
            for kk = 1:numel(step5)
                step6 = regexp(step5{kk},'\$','split');
                track1(end+1,1:2) = step6;
                track2{end+1,1} = step5{kk};

                var = data.(step6{1});
                if isnumeric(var)
                    cond = [cond demean(var)];
                    flag(kk) = 0;
                else
                    flag(kk) = 1;
                    if iscell(var)
                        tmp = zeros(numel(var),1);
                        i1 = strmatch(step6{2},var,'exact');
                        tmp(i1)=1;
                    else
                        tmp = data.(step6{1})==step6{2};
                    end
                    cond = [cond tmp];
                end
            end

            if any(flag)
                cond2 = [cond2 prod(cond,2)];
            else
                cond2 = [cond2 cond];
            end

            if all(flag)
                cond3 = [cond3 sum(cond2,2)>0]; 
            else
                cond3 = [cond3 cond2];
            end
        end
        cond4 = [cond4 cond3];

        if ~isempty(contains('#',step2(mm)))
            extra = [extra mean(cond3,2)];
        end
    end
    con = cond4;


    tx1 = X;
    tx2 = X;
    
    tmp = [];
    tr = [];
    
    for ii = 1:size(con,2);
        ind = find(con(:,ii)~=0);
        tmp{ii} = ind;
        in = [find(mean(X(ind,:)==repmat(con(ind,ii),1,size(X,2)))==1) find(mean(-X(ind,:)==repmat(con(ind,ii),1,size(X,2)))==1)];
        tr = [tr in];
        tx2(ind,in)=0;
    end
    
    tr = unique(tr);
    if isempty(tr)
        error('Invalid Contrast')
    end

    tx2 = [tx2 extra];
    
    %%%
    [a i] = unique(track2);
    track = track1(i,:);

    list = unique(track(:,1));
    wh = [];
    for ii = 1:numel(list)
        cc = numel(strmatch(list{ii},track1(:,1),'exact'));
        if cc>1
            if ~isempty(contains(list{ii},ErrTerms))
                wh{end+1} = contains(list{ii},ErrTerms);
            end
        end
    end

    if isempty(wh) || isempty(wh{1})
        ET = 1;
    else
        ET = intersections(wh);
        ET = ET(1);
    end
    df = ResidualDFs(tx1)-ResidualDFs(tx2);
end