function [F p perms] = matGLM(MAT,mod,effect,errorTerm, runPerms)

F=NaN; p=NaN; perms = NaN;

if nargin<5
    runPerms=0;
end
% errorTerm = 1;
% effect = ?;
% mod = ANOVA_APS(td,'y~T807_middletemporal_bh_FS+PIB',3);
% keyboard;

[RM r1 r2] = makeSSmat(mod.RFMs(errorTerm).tx1, mod.RFMs(errorTerm).tx2);
ESS = LoopEstimate(MAT,1,RM);
EDF = mod.RFMs.EDF;


if isempty(mod.RFMs(errorTerm).Effect(effect).tx2)
    [EM r1 r2] = makeSSmat(mod.RFMs(errorTerm).Effect(effect).tx1,  zeros(size(MAT,1),1));
else
    [EM r1 r2] = makeSSmat(mod.RFMs(errorTerm).Effect(effect).tx1, mod.RFMs(errorTerm).Effect(effect).tx2);
end

SS = LoopEstimate(MAT,1,EM);
DF = mod.RFMs(1).Effect(1).df;

F = (SS./DF)./(ESS./EDF);
p = 1-cdf('f',F,DF,EDF);

if runPerms>0
    perms = zeros(runPerms,size(MAT,2));
    for ii = 1:runPerms
        Yp=MAT(randperm(size(MAT,1)),:);
        
        RSS = LoopEstimate(Yp,1,RM);
        ESS = LoopEstimate(Yp,1,EM);
        FF = (ESS/DF)./(RSS/EDF);
        
        perms(ii,:) = 1-cdf('f',FF,DF,EDF);
    end
end